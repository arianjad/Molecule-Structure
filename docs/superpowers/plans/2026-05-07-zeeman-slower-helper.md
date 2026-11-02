# RaF Zeeman-Slower Helper Implementation Plan (consolidated)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers-extended-cc:subagent-driven-development (recommended) or superpowers-extended-cc:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a notebook helper to compute and plot transition energies, spontaneous-emission branching ratios, and laser-driven excitation strengths vs B field for selected RaF X→A level pairs (Zeeman-slower study).

**Architecture:** Compute/plot split (Approach B from the design spec). One `compute_transitions_vs_B` function loops over the precomputed `ZeemanMap` eigensystems and returns a data dict; three independent `plot_*` functions consume the dict and emit one matplotlib figure each. Lives in a fresh notebook; does not touch `Source Code/` or other notebooks.

**Tech Stack:** Python 3, NumPy, matplotlib, the local `MoleculeLevels` library (`Source Code/Energy_Levels.py`), Jupyter (conda env `Structure`).

**Design spec:** `~/.claude/plans/enchanted-sprouting-scroll.md` (approved). Refer to it for the full data dict schema, polarization mapping table, and physics conventions.

**Verification convention:** Per `Jupyter Notebooks/CLAUDE.md`, notebook-only changes verify via "Restart & Run All" on the affected notebook. Each task's verification cell is a Python assertion run in the notebook; if it raises, the task isn't done. There is no pytest infrastructure in this repo. The verify command is:

```bash
conda run -n Structure jupyter execute "Jupyter Notebooks/RaX/RaF_Zeeman_Slower.ipynb"
```

**Commit convention:** Per `~/.claude/rules/code-style.md` and observed `git log`: commit early and often, concise subject, no `Co-Authored-By` trailers (the repo doesn't use them).

**Current state (2026-05-07, end of session):** All three tasks complete. Three commits on `fix-for-distrib`:
- `1b83e0f` Add RaF Zeeman slower scaffold + compute_transitions_vs_B
- `a322ce8` Add plot_transition_energies, plot_branching_ratios, plot_excitation_strength
- `31f7dce` Optimize compute fn (3× faster, no print spam); add Task C example + cross-check; coarse Bz=21 for fast verify

**Notable deviations from the original plan, decided during execution:**
- **Bz=21 instead of 401** in cell 4. Coarse for fast verify (~15s vs >25min before optimization). Bump to 401+ for smooth demo plots; the helpers don't care.
- **Cycling-line query** (`X(N=1, parity='-') → A(J=1/2, parity='+')`) instead of the original spec's `{'N':[0]}` + parity_e='+'. The original query produces all-zero TDMs by parity selection rule (N=0 has parity +, A(+) requires parity flip). Caught when smoke-testing during Task C.
- **`compute_transitions_vs_B` no longer calls `Calculate_TDM_evecs`**; it inlines the case conversion (hoisted outside the per-`p` loop) and calls `e.library.TDM_p_builders[...]` directly. Equivalent math (verified bit-exact on a small grid), but 3× fewer convert_evecs calls and stdout-suppressed (the Source Code function prints by default and `Calculate_TDM_evecs` doesn't propagate `verbose=False`).
- **Sanity check (cell 5b)** got an extra assertion `_data_test['TDM_pol_sq'].max() > 1e-3` to catch the all-noise case (the closure `BR.sum + leakage = 1` is tautological by construction and would have passed even with all-zero TDMs).
- **`plot_excitation_strength`** had an unconditional `ax.legend()` that triggered a UserWarning when a ground state had no allowed pairs; gated behind `if partners`.

Verify command: `conda run -n Structure jupyter nbconvert --to notebook --execute --inplace "Jupyter Notebooks/RaX/RaF_Zeeman_Slower.ipynb"` → exit 0; sanity prints from cells 5/7/9/11/13.

---

## File Structure

**Create / modify**:
- `Jupyter Notebooks/RaX/RaF_Zeeman_Slower.ipynb` — the only file. Final cell layout:
  1. Path boilerplate (`config_path.add_to_sys_path()`)  *(already present)*
  2. Imports  *(already present)*
  3. State init for `g` (RaF X) and `e` (RaF A); `g.eigensystem(0,1e-3)`, `e.eigensystem(0,1e-3)`  *(already present)*
  4. B-field grid + `g.ZeemanMap(Bz)`, `e.ZeemanMap(Bz)` + shape assertions  *(already present)*
  5. `compute_transitions_vs_B` definition + sanity-check assertion  *(Task A)*
  6. `plot_transition_energies` definition + sanity-check assertion  *(Task B)*
  7. `plot_branching_ratios` definition + sanity-check assertion  *(Task B)*
  8. `plot_excitation_strength` definition + sanity-check assertion  *(Task B)*
  9. Example call: `compute_transitions_vs_B(...)` + render three figures  *(Task C)*
  10. Cross-check vs `branching_ratios` reference (print, do not assert)  *(Task C)*

---

## Task A: Notebook scaffold + `compute_transitions_vs_B`

**Goal:** Verified, committed notebook scaffold (cells 1–4) plus the data-producing function `compute_transitions_vs_B` (cell 5) that the three plot functions will consume.

**Files:**
- Modify: `Jupyter Notebooks/RaX/RaF_Zeeman_Slower.ipynb` (currently exists at cells 1–4 only).

**Acceptance Criteria:**
- [ ] Cells 1–5 plus the sanity-check cell execute on a fresh kernel under `conda run -n Structure jupyter execute` with exit 0 and no tracebacks.
- [ ] After cell 4: `g.evals_B.shape == (len(Bz), n_g_levels)`, `e.evals_B.shape == (len(Bz), n_e_levels)`, `g.evecs_B.ndim == 3`, `e.evecs_B.ndim == 3`.
- [ ] After the sanity-check cell: `compute_transitions_vs_B(...)` returns dict with keys `Bz, g_idx, e_idx, dE, BR, TDM_pol_sq, leakage, keep_pair, g_label, e_label, pair_label, meta`.
- [ ] `dE.shape == (N_B, n_e_sel, n_g_sel)`, `BR.shape == (N_B, n_e_sel, n_g_sel)`, `TDM_pol_sq.shape == (N_B, n_e_sel, n_g_sel)`, `leakage.shape == (N_B, n_e_sel)`.
- [ ] `np.all(BR >= 0)` and `np.all(BR <= 1 + 1e-9)`.
- [ ] `np.allclose(BR.sum(axis=2) + leakage, 1.0, atol=1e-6)` per excited state at every B (closure).

**Verify:** `conda run -n Structure jupyter execute "Jupyter Notebooks/RaX/RaF_Zeeman_Slower.ipynb"` → exit 0; sanity-check cell prints `compute_transitions_vs_B: shapes OK, BR closure OK`.

**Steps:**

- [ ] **A.1: Confirm cells 1–4 are present and unmodified.** They were written previously and match the spec verbatim. If anything looks off, restore them to:

Cell 1:
```python
from config_path import add_to_sys_path
add_to_sys_path()
```

Cell 2:
```python
import numpy as np
import matplotlib.pyplot as plt
from Energy_Levels import MoleculeLevels, Calculate_TDM_evecs, branching_ratios
```

Cell 3:
```python
N_g = np.arange(0, 4)
N_e = np.arange(1, 5)

g = MoleculeLevels.initialize_state(
    molecule_name='RaF', elec_state='X', vib_state=0,
    N_list=N_g, fermion_or_boson='boson',
    M_sublevels='all', M_list=[1/2],
    I_nuclei=[0, 1/2], isotope=226, round=8,
    params=None, P_values=[1/2],
)
e = MoleculeLevels.initialize_state(
    molecule_name='RaF', elec_state='A', vib_state=0,
    N_list=N_e, fermion_or_boson='boson',
    M_sublevels='all', M_list=[1/2],
    I_nuclei=[0, 1/2], isotope=226, round=8,
    params=None, P_values=[1/2],
)
g.eigensystem(0, 1e-3)
e.eigensystem(0, 1e-3)
```

Cell 4:
```python
Bz = np.linspace(1e-6, 100, 401)   # Gauss; coarser than RaF X-A's 4000 to keep first run fast
g.ZeemanMap(Bz, plot=False)
e.ZeemanMap(Bz, plot=False)

assert g.evals_B.shape[0] == len(Bz) and g.evals_B.ndim == 2
assert e.evals_B.shape[0] == len(Bz) and e.evals_B.ndim == 2
assert g.evecs_B.ndim == 3 and e.evecs_B.ndim == 3
```

- [ ] **A.2: Add cell 5 with the compute function and helpers.**

```python
def _polarization_to_p_weights(polarization):
    """Map polarization label → dict {p: weight} for coherent superposition.

    Returns weights for p in {-1, 0, +1}. Used by compute_transitions_vs_B
    to build the per-polarization TDM from per-p TDMs.
    """
    if isinstance(polarization, dict):
        out = {-1: 0.0+0j, 0: 0.0+0j, +1: 0.0+0j}
        for k, v in polarization.items():
            kk = int(k) if isinstance(k, str) else k
            out[kk] = complex(v)
        return out
    s = polarization.lower()
    if s in ('sigma+', 'sigma_plus', 's+'):  return {-1: 0+0j, 0: 0+0j, +1: 1+0j}
    if s in ('sigma-', 'sigma_minus', 's-'): return {-1: 1+0j, 0: 0+0j, +1: 0+0j}
    if s in ('pi', 'z'):                     return {-1: 0+0j, 0: 1+0j, +1: 0+0j}
    if s == 'x':
        return {-1: 1/np.sqrt(2)+0j, 0: 0+0j, +1: -1/np.sqrt(2)+0j}
    if s == 'y':
        return {-1: 1j/np.sqrt(2),    0: 0+0j, +1: 1j/np.sqrt(2)}
    raise ValueError(f"Unknown polarization {polarization!r}")


def _qn_label(state, idx):
    """Build a short label like 'N=0,J=1/2,F=1,M=+1/2' from dominant-basis QNs."""
    qn = state.q_numbers
    parts = []
    for key in ('N', 'J', 'F', 'M'):
        if key in qn:
            v = qn[key][idx]
            if abs(v - round(v)) < 1e-6:
                parts.append(f"{key}={int(round(v))}")
            else:
                parts.append(f"{key}={int(round(2*v))}/2")
    if hasattr(state, 'parities'):
        parts.append('+' if state.parities[idx] > 0 else '-')
    return ','.join(parts)


def compute_transitions_vs_B(g, e, ground_query, excited_query,
                             polarization='sigma+',
                             q_electronic=(-1, +1),
                             parity_g=None, parity_e=None,
                             tdm_thresh=1e-6):
    """Compute transition energies, spontaneous-emission BR, and excitation
    strengths vs B for selected ground/excited level pairs.

    g, e must have ZeemanMap(Bz) already run. q_electronic is the molecule-frame
    component set fixed by the electronic transition character (default (-1,+1)
    for X-A perpendicular). polarization is the lab-frame light polarization.

    Returns dict — see plan / design spec.
    """
    Bz = g.Bz
    assert np.array_equal(Bz, e.Bz), "g and e must share the same Bz grid"
    N_B = len(Bz)

    g_idx = g.select_q(ground_query, parity=parity_g)
    e_idx = e.select_q(excited_query, parity=parity_e)
    if len(g_idx) == 0 or len(e_idx) == 0:
        raise ValueError(f"Empty selection: g_idx={g_idx}, e_idx={e_idx}")

    n_g_full = g.evecs_B.shape[1]
    n_e_full = e.evecs_B.shape[1]
    p_list = (-1, 0, +1)
    TDM_full = np.zeros((N_B, 3, n_e_full, n_g_full), dtype=complex)
    for ip, p in enumerate(p_list):
        for i in range(N_B):
            TDM_full[i, ip] = Calculate_TDM_evecs(
                p, g.evecs_B[i], g, e.evecs_B[i], e, q=list(q_electronic),
            )

    TDM_sel = TDM_full[:, :, e_idx[:, None], g_idx[None, :]]   # (N_B, 3, n_e_sel, n_g_sel)

    # Spontaneous-emission BR: sum |TDM_p|^2 over photon polarization p
    num   = (np.abs(TDM_sel)**2).sum(axis=1)                       # (N_B, n_e_sel, n_g_sel)
    denom = (np.abs(TDM_full[:, :, e_idx, :])**2).sum(axis=(1, 3)) # (N_B, n_e_sel)
    BR = num / np.where(denom[:, :, None] > 0, denom[:, :, None], 1.0)
    leakage = 1.0 - BR.sum(axis=2)

    weights = _polarization_to_p_weights(polarization)
    TDM_pol = sum(weights[p] * TDM_sel[:, ip] for ip, p in enumerate(p_list))
    TDM_pol_sq = np.abs(TDM_pol)**2

    dE = e.evals_B[:, e_idx, None] - g.evals_B[:, None, g_idx]

    keep_pair = TDM_pol_sq.max(axis=0) >= tdm_thresh   # (n_e_sel, n_g_sel)

    g_label = {gi: _qn_label(g, gi) for gi in g_idx}
    e_label = {ej: _qn_label(e, ej) for ej in e_idx}
    pair_label = {(ej, gi): f"g({g_label[gi]}) → e({e_label[ej]})"
                  for ej in e_idx for gi in g_idx}

    return {
        'Bz': Bz,
        'g_idx': np.asarray(g_idx),
        'e_idx': np.asarray(e_idx),
        'dE': dE,
        'BR': BR,
        'TDM_pol_sq': TDM_pol_sq,
        'leakage': leakage,
        'keep_pair': keep_pair,
        'g_label': g_label,
        'e_label': e_label,
        'pair_label': pair_label,
        'meta': dict(polarization=polarization, q_electronic=tuple(q_electronic),
                     parity_g=parity_g, parity_e=parity_e, tdm_thresh=tdm_thresh),
    }
```

- [ ] **A.3: Add a sanity-check cell (cell 5b) right after the function definition.**

```python
_data_test = compute_transitions_vs_B(
    g, e, {'N': [0]}, {'J': 0.5},
    polarization='sigma+', q_electronic=(-1, +1), parity_e='+',
)
N_B = len(Bz)
n_e_sel = len(_data_test['e_idx'])
n_g_sel = len(_data_test['g_idx'])
assert _data_test['dE'].shape         == (N_B, n_e_sel, n_g_sel)
assert _data_test['BR'].shape         == (N_B, n_e_sel, n_g_sel)
assert _data_test['TDM_pol_sq'].shape == (N_B, n_e_sel, n_g_sel)
assert _data_test['leakage'].shape    == (N_B, n_e_sel)
assert (_data_test['BR'] >= 0).all() and (_data_test['BR'] <= 1 + 1e-9).all()
closure = _data_test['BR'].sum(axis=2) + _data_test['leakage']
assert np.allclose(closure, 1.0, atol=1e-6), f"closure max err {np.max(np.abs(closure-1))}"
print("compute_transitions_vs_B: shapes OK, BR closure OK")
```

- [ ] **A.4: Run the verify command.** Expect exit 0 and the printed line.

- [ ] **A.5: Commit.**

```bash
git add "Jupyter Notebooks/RaX/RaF_Zeeman_Slower.ipynb"
git commit -m "Add RaF Zeeman slower scaffold + compute_transitions_vs_B"
```

---

## Task B: Three plot functions

**Goal:** All three plot functions, each with a sanity-check cell.

**Files:**
- Modify: `Jupyter Notebooks/RaX/RaF_Zeeman_Slower.ipynb` (add cells 6–8 + their sanity checks).

**Acceptance Criteria:**
- [ ] `plot_transition_energies(_data_test)` returns `(fig, panels)`; `len(panels) == len(e_idx)`; each panel's `y.shape[0] == len(Bz)`.
- [ ] `plot_branching_ratios(_data_test)` returns `(fig, panels)`; each panel includes a `leakage` trace (1-D, length `len(Bz)`).
- [ ] `plot_excitation_strength(_data_test)` returns `(fig, panels)`; `len(panels) == len(g_idx)`.
- [ ] All three sanity-check cells call `plt.close(fig)` to keep `jupyter execute` clean.
- [ ] `conda run -n Structure jupyter execute ...` exits 0; sanity prints appear for each function.

**Verify:** `conda run -n Structure jupyter execute "Jupyter Notebooks/RaX/RaF_Zeeman_Slower.ipynb"` → exit 0; sees three printed lines `plot_transition_energies: OK`, `plot_branching_ratios: OK`, `plot_excitation_strength: OK`.

**Steps:**

- [ ] **B.1: Add cell 6 — `plot_transition_energies` with shared `_grid_shape` helper.**

```python
def _grid_shape(n):
    cols = int(np.ceil(np.sqrt(n)))
    rows = int(np.ceil(n / cols))
    return rows, cols


def plot_transition_energies(data, GHz=False, kG=False, figsize_per=(5, 4)):
    """One panel per excited level. Lines = ground partners. y = absolute ΔE."""
    Bz = data['Bz']; e_idx = data['e_idx']; g_idx = data['g_idx']
    dE = data['dE']; keep = data['keep_pair']
    pair_label = data['pair_label']; e_label = data['e_label']

    x = Bz * (1e-3 if kG else 1.0)
    yscale = 1e-3 if GHz else 1.0
    xlabel = 'B (kGauss)' if kG else 'B (Gauss)'
    ylabel = 'ΔE (GHz)' if GHz else 'ΔE (MHz)'

    nrows, ncols = _grid_shape(len(e_idx))
    fig, axes = plt.subplots(nrows, ncols, squeeze=False,
                             figsize=(ncols*figsize_per[0], nrows*figsize_per[1]))
    panels = []
    for k, ej in enumerate(e_idx):
        ax = axes[k // ncols, k % ncols]
        partners = [gi for gi in g_idx if keep[k, list(g_idx).index(gi)]]
        if not partners:
            ax.set_title(f"e({e_label[ej]})  (no allowed pairs)")
            ax.axis('off')
            panels.append({'x': x, 'y': np.empty((len(x), 0)),
                           'e_idx': ej, 'partner_idx': [], 'label': e_label[ej], 'ax': ax})
            continue
        partner_cols = [list(g_idx).index(gi) for gi in partners]
        y = yscale * dE[:, k, partner_cols]
        for col, gi in zip(range(y.shape[1]), partners):
            ax.plot(x, y[:, col], lw=1, label=pair_label[(ej, gi)])
        ax.set_title(f"e({e_label[ej]})", fontsize=10)
        ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
        ax.legend(fontsize=7, loc='best')
        panels.append({'x': x, 'y': y, 'e_idx': ej,
                       'partner_idx': partners, 'label': e_label[ej], 'ax': ax})
    for k in range(len(e_idx), nrows*ncols):
        axes[k // ncols, k % ncols].axis('off')
    fig.suptitle(f"Transition energies — pol={data['meta']['polarization']}, q={data['meta']['q_electronic']}",
                 fontsize=12)
    fig.tight_layout()
    return fig, panels
```

Sanity-check cell (cell 6b):
```python
_fig, _panels = plot_transition_energies(_data_test)
assert len(_panels) == len(_data_test['e_idx'])
for p in _panels:
    assert p['y'].shape[0] == len(_data_test['Bz'])
plt.close(_fig)
print("plot_transition_energies: OK")
```

- [ ] **B.2: Add cell 7 — `plot_branching_ratios`.**

```python
def plot_branching_ratios(data, kG=False, figsize_per=(5, 4)):
    """One panel per excited level. Lines = BR(e -> g) per selected ground.
    Dashed line per panel = leakage (1 - sum_g BR over selected g)."""
    Bz = data['Bz']; e_idx = data['e_idx']; g_idx = data['g_idx']
    BR = data['BR']; leakage = data['leakage']; keep = data['keep_pair']
    pair_label = data['pair_label']; e_label = data['e_label']

    x = Bz * (1e-3 if kG else 1.0)
    xlabel = 'B (kGauss)' if kG else 'B (Gauss)'

    nrows, ncols = _grid_shape(len(e_idx))
    fig, axes = plt.subplots(nrows, ncols, squeeze=False,
                             figsize=(ncols*figsize_per[0], nrows*figsize_per[1]))
    panels = []
    for k, ej in enumerate(e_idx):
        ax = axes[k // ncols, k % ncols]
        partners = [gi for gi in g_idx if keep[k, list(g_idx).index(gi)]]
        partner_cols = [list(g_idx).index(gi) for gi in partners]
        if partners:
            y = BR[:, k, partner_cols]
            for col, gi in zip(range(y.shape[1]), partners):
                ax.plot(x, y[:, col], lw=1, label=pair_label[(ej, gi)])
        ax.plot(x, leakage[:, k], lw=1, ls='--', color='k', label='leakage')
        ax.set_title(f"e({e_label[ej]})", fontsize=10)
        ax.set_xlabel(xlabel); ax.set_ylabel('Branching ratio')
        ax.set_ylim(0, 1.05)
        ax.legend(fontsize=7, loc='best')
        panels.append({'x': x,
                       'y': BR[:, k, partner_cols] if partners else np.empty((len(x), 0)),
                       'leakage': leakage[:, k],
                       'e_idx': ej, 'partner_idx': partners,
                       'label': e_label[ej], 'ax': ax})
    for k in range(len(e_idx), nrows*ncols):
        axes[k // ncols, k % ncols].axis('off')
    fig.suptitle(f"Spontaneous-emission BR — q={data['meta']['q_electronic']}", fontsize=12)
    fig.tight_layout()
    return fig, panels
```

Sanity-check cell (cell 7b):
```python
_fig, _panels = plot_branching_ratios(_data_test)
assert len(_panels) == len(_data_test['e_idx'])
for p in _panels:
    assert p['leakage'].shape[0] == len(_data_test['Bz'])
plt.close(_fig)
print("plot_branching_ratios: OK")
```

- [ ] **B.3: Add cell 8 — `plot_excitation_strength`.**

```python
def plot_excitation_strength(data, kG=False, figsize_per=(5, 4)):
    """One panel per ground level. Lines = |TDM_polarization|^2 to excited partners."""
    Bz = data['Bz']; g_idx = data['g_idx']; e_idx = data['e_idx']
    TDM2 = data['TDM_pol_sq']; keep = data['keep_pair']
    pair_label = data['pair_label']; g_label = data['g_label']

    x = Bz * (1e-3 if kG else 1.0)
    xlabel = 'B (kGauss)' if kG else 'B (Gauss)'
    pol = data['meta']['polarization']

    nrows, ncols = _grid_shape(len(g_idx))
    fig, axes = plt.subplots(nrows, ncols, squeeze=False,
                             figsize=(ncols*figsize_per[0], nrows*figsize_per[1]))
    panels = []
    for k, gi in enumerate(g_idx):
        ax = axes[k // ncols, k % ncols]
        partners = [ej for ej in e_idx if keep[list(e_idx).index(ej), k]]
        partner_rows = [list(e_idx).index(ej) for ej in partners]
        if partners:
            y = TDM2[:, partner_rows, k]
            for col, ej in zip(range(y.shape[1]), partners):
                ax.plot(x, y[:, col], lw=1, label=pair_label[(ej, gi)])
        else:
            ax.text(0.5, 0.5, 'no allowed pairs', ha='center', va='center',
                    transform=ax.transAxes)
            y = np.empty((len(x), 0))
        ax.set_title(f"g({g_label[gi]})", fontsize=10)
        ax.set_xlabel(xlabel); ax.set_ylabel(r'$|\mathrm{TDM}_{\mathrm{pol}}|^2$')
        ax.legend(fontsize=7, loc='best')
        panels.append({'x': x, 'y': y, 'g_idx': gi,
                       'partner_idx': partners, 'label': g_label[gi], 'ax': ax})
    for k in range(len(g_idx), nrows*ncols):
        axes[k // ncols, k % ncols].axis('off')
    fig.suptitle(f"Excitation strength — pol={pol}, q={data['meta']['q_electronic']}",
                 fontsize=12)
    fig.tight_layout()
    return fig, panels
```

Sanity-check cell (cell 8b):
```python
_fig, _panels = plot_excitation_strength(_data_test)
assert len(_panels) == len(_data_test['g_idx'])
plt.close(_fig)
print("plot_excitation_strength: OK")
```

- [ ] **B.4: Run the verify command.** Expect exit 0; three sanity prints appear.

- [ ] **B.5: Commit.**

```bash
git add "Jupyter Notebooks/RaX/RaF_Zeeman_Slower.ipynb"
git commit -m "Add plot_transition_energies, plot_branching_ratios, plot_excitation_strength"
```

---

## Task C: Example call + cross-check

**Goal:** A user-facing example that calls `compute_transitions_vs_B` and renders all three figures, plus a cross-check cell printing (not asserting) the new BR alongside the scalar-Bz reference (`branching_ratios` from `Source Code/Energy_Levels.py:1236`).

**Files:**
- Modify: `Jupyter Notebooks/RaX/RaF_Zeeman_Slower.ipynb` (add cells 9–10).

**Acceptance Criteria:**
- [ ] Cell 9 produces three inline figures.
- [ ] Cell 10 prints `BR cross-check max err at B=... G: ...` plus the first row of both BRs and the leakage. No hard assertion.
- [ ] End-to-end `jupyter execute` exits 0.
- [ ] (Manual) Visual inspection in JupyterLab: energy panels show smooth ΔE vs B; BR panels' traces ≤1 with visible leakage trace; excitation panels show one curve per ground state with its excited partners.

**Verify:** `conda run -n Structure jupyter execute "Jupyter Notebooks/RaX/RaF_Zeeman_Slower.ipynb"` → exit 0. Then open in JupyterLab and eyeball the three figures.

**Steps:**

- [ ] **C.1: Add cell 9 — example call.**

```python
data = compute_transitions_vs_B(
    g, e,
    ground_query={'N': [0]},
    excited_query={'J': 0.5},
    polarization='sigma+',
    q_electronic=(-1, +1),
    parity_e='+',
    tdm_thresh=1e-6,
)
fig_E,   _ = plot_transition_energies(data)
fig_BR,  _ = plot_branching_ratios(data)
fig_abs, _ = plot_excitation_strength(data)
plt.show()
```

- [ ] **C.2: Add cell 10 — cross-check.**

```python
# Cross-check: at one B value, our normalized BR should be comparable to a
# renormalized branching_ratios() reference. Strict equality is NOT guaranteed
# because branching_ratios uses TDM_builders (not TDM_p_builders) and re-
# diagonalizes internally; print and inspect.
i_check = len(Bz) // 2
Bz_val  = float(Bz[i_check])
BR_ref_full = branching_ratios(g, e, 0, Bz_val)   # shape (n_g_full, n_e_full)

e_idx = data['e_idx']; g_idx = data['g_idx']
BR_ref_sel = BR_ref_full[g_idx[:, None], e_idx[None, :]].T   # (n_e_sel, n_g_sel)
denom_ref = BR_ref_full[:, e_idx].sum(axis=0)                # (n_e_sel,)
BR_ref_norm = BR_ref_sel / np.where(denom_ref[:, None] > 0, denom_ref[:, None], 1.0)

err = np.max(np.abs(data['BR'][i_check] - BR_ref_norm))
print(f"BR cross-check max err at B={Bz_val:.4g} G: {err:.2e}")
print(f"  data BR[0]:    {np.array2string(data['BR'][i_check, 0], precision=4)}")
print(f"  ref  BR[0]:    {np.array2string(BR_ref_norm[0],         precision=4)}")
print(f"  data leakage:  {data['leakage'][i_check]}")
```

- [ ] **C.3: Run the verify command.** Expect exit 0; prints appear.

- [ ] **C.4: Open `Jupyter Notebooks/RaX/RaF_Zeeman_Slower.ipynb` in JupyterLab and visually inspect the three figures.** If something looks off (sharp spikes that aren't avoided crossings, missing traces, unexpected axes), STOP and ask the user before claiming done.

- [ ] **C.5: Commit.**

```bash
git add "Jupyter Notebooks/RaX/RaF_Zeeman_Slower.ipynb"
git commit -m "Add example call + BR cross-check vs scalar-B reference"
```

---

## After all tasks complete

Ask the user which physics-verification checks they want to run beyond the BR closure already built into Task A. Per `~/.claude/rules/physics-verification.md`, candidates:

- ΔE at B≈0 vs published RaF X→A line center for selected (J, F).
- dν/dB slopes at low B vs `g_e − g_g` × Bohr magneton; cross-reference `plot_zeeman_shifts(..., derivative=1, g_factor=True)` already in `RaF X-A.ipynb`.
- Try a non-σ⁺ polarization (`'pi'`, `'x'`) and confirm the ΔM filter activates the right pairs.
- Inspect the cross-check delta from Task C — if `BR_ref_norm` and `data['BR']` agree to within ~few %, the BR math is consistent with the scalar reference; if they differ by more than that, dig into the q-convention difference.
