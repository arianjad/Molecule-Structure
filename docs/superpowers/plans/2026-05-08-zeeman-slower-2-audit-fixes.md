# `zeeman_slower_2.ipynb` Audit-Fix Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers-extended-cc:subagent-driven-development (recommended) or superpowers-extended-cc:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Apply the six top fixes from the 2026-05-07 audit of `Jupyter Notebooks/RaX/zeeman_slower_2.ipynb` — performance, physics-correctness diagnostics, and methodological clarity. No new physics, no algorithmic changes beyond what each fix specifies.

**Architecture:** All edits land in a single notebook (`zeeman_slower_2.ipynb`). Per-cell edits are isolated; only one task changes the helper-function cell that downstream cells depend on (Task 1). Verify each task with a `/tmp/` smoke script that asserts the targeted behavior, then a final end-to-end notebook execute confirms integration.

**Tech Stack:** Python (conda env `Structure`), numpy, matplotlib, the `Energy_Levels` module from `Source Code/`, jupyter for notebook editing/execution.

**Conventions reminder (from prior helper plan):**
- Every cell-edit script in `/tmp/` MUST first re-strip notebook outputs (`cell["outputs"]=[]; cell["execution_count"]=None`) so commits show only source diffs.
- Use `python3 - <<'PY' ... PY` pattern with REPL `sh()` to patch the notebook JSON; do NOT use `nbconvert --to notebook --execute --inplace` for edits.
- Never have two notebook executions running concurrently — kill prior `jupyter execute` / `nbconvert` processes before starting a new run (`pkill -f 'jupyter.*zeeman_slower_2'`).
- Execute via `conda run -n Structure jupyter nbconvert --to notebook --execute --inplace "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"`.

---

## File map

| File | Role | Touched by |
|---|---|---|
| `Jupyter Notebooks/RaX/zeeman_slower_2.ipynb` | The target — all source edits land here | All tasks |
| `Source Code/molecule_parameters.py` | Read-only — source of `Gamma_SR`, `bF` for Task 4 | Task 4 (read only) |
| `/tmp/<task>_patch.py`, `/tmp/<task>_verify.py` | Throwaway scripts, one pair per task | All tasks |

**No new files. No edits to `Source Code/`.**

---

## Cell index reference (current state, 31 cells)

| Cell | Type | Content |
|---|---|---|
| 0 | md | Header / conventions |
| 1 | code | `from config_path import add_to_sys_path; add_to_sys_path()` |
| 2 | code | Constants — has `a_recoil` bug, redundant `branching_ratios` import |
| 3 | md | "Step 1 — build state objects" |
| 4 | code | `MoleculeLevels.initialize_state` for `g`, `e` |
| 5 | code | `Bz = np.linspace(1e-3, 1000.0, 101)`, `ZeemanMap` |
| 6 | md | "Step 2 — full transition map" |
| 7 | code | **`build_transition_map`** — print-spam loop |
| 8 | code | `TDM2_ref` normalization — has misleading "case-b stretched" label |
| 9 | md | "Ground-state character labels" |
| 10 | code | `ground_character` |
| 11 | md | "Subset utility" |
| 12 | code | `subset_TM` + call |
| 13 | md | "Step 3 — diagnostic plots" |
| 14 | md | "Plot A — transition frequencies" |
| 15 | code | **`plot_A_transition_frequencies`** — has hardcoded `B_SR_PB`/`B_HF_PB` |
| 16 | md | "Plot B — slope vs B, Doppler-trackable bands" |
| 17 | code | **`find_doppler_trackable_bands`** — slope-only clustering; `plot_B_slopes` |
| 18 | md | "Plot C — required B(z)" |
| 19 | code | `slowing_profile`, `required_Bz`, `a_max_est` |
| 20 | code | `plot_C_required_Bz` + call |
| 21 | md | "Plot D — branching-ratio map" |
| 22 | code | `compute_branching` + call |
| 23 | code | `plot_D_branching` + call |
| 24 | md | "Plot E — dark-state count" |
| 25 | code | `dark_state_count`, `plot_E_dark_states` + call |
| 26 | md | **"Step 4 — rate-equation force"** — needs caveat |
| 27 | code | **`steady_state_rate_eq`** — needs leakage diagnostic; smoke test |
| 28 | code | `plot_force_vs_v` + call |
| 29 | code | `F_max_vs_B` + call |
| 30 | md | "Notes & next steps" |

---

## Self-contained `/tmp/` patch idiom

Every code-cell edit follows this template. Save as `/tmp/<task_name>_patch.py`:

```python
import json, pathlib
p = pathlib.Path("Jupyter Notebooks/RaX/zeeman_slower_2.ipynb")
nb = json.loads(p.read_text())

def src_to_list(s):
    lines = s.split("\n")
    return [l + "\n" for l in lines[:-1]] + ([lines[-1]] if lines[-1] else [])

# Strip outputs from every code cell so the diff shows only source changes
for c in nb["cells"]:
    if c["cell_type"] == "code":
        c["outputs"] = []
        c["execution_count"] = None

# === per-task replacement(s) here ===
old_src = "".join(nb["cells"][CELL_IDX]["source"])
new_src = old_src.replace("OLD_SNIPPET", "NEW_SNIPPET")
assert new_src != old_src, "patch did not match — check OLD_SNIPPET"
nb["cells"][CELL_IDX]["source"] = src_to_list(new_src)

p.write_text(json.dumps(nb, indent=1) + "\n")
print("patched cell", CELL_IDX)
print("--- new cell preview ---")
print("".join(nb["cells"][CELL_IDX]["source"])[:600])
```

For tasks that REPLACE an entire cell, use `nb["cells"][CELL_IDX]["source"] = src_to_list(NEW_FULL_SOURCE)` instead of `.replace()`.

---

## Task 1: Optimize `build_transition_map` — hoist case-conversion, suppress prints

**Goal:** Eliminate the ~600 "Successfully converted eigenvectors" prints during the transition-map build, and cut runtime ~3× by hoisting the case-(b)→case-(a) conversion outside the per-(B, p) loop.

**Files:**
- Modify: `Jupyter Notebooks/RaX/zeeman_slower_2.ipynb` cell 7 (`build_transition_map`)

**Acceptance Criteria:**
- [ ] Notebook execution emits zero "Successfully converted eigenvectors" lines from cell 7.
- [ ] New `build_transition_map` produces `TDM_full` bit-identical to the old version on `Bz=np.linspace(1e-3, 100, 11)` (max abs diff = 0).
- [ ] Cell 7 returns same dict keys, same shapes, same dtypes.
- [ ] Public signature unchanged (same args, same defaults).

**Verify:** Run `/tmp/task1_verify.py` (defined below) → prints `BIT-EXACT MATCH; 0 prints captured`.

**Steps:**

- [ ] **Step 1: Write the failing equivalence test**

Save as `/tmp/task1_verify.py`:

```python
import sys, io, contextlib, pathlib, importlib.util
sys.path.insert(0, "Source Code")
from config_path import add_to_sys_path; add_to_sys_path()
import numpy as np
from Energy_Levels import MoleculeLevels, Calculate_TDM_evecs

# Build small grid for fast equivalence
g = MoleculeLevels.initialize_state(
    molecule_name='RaF', elec_state='X', vib_state=0,
    N_list=[1], fermion_or_boson='boson',
    M_sublevels='all', I_nuclei=[0, 1/2], isotope=226,
    round=8, P_values=[1/2])
e = MoleculeLevels.initialize_state(
    molecule_name='RaF', elec_state='A', vib_state=0,
    N_list=[1, 2], fermion_or_boson='boson',
    M_sublevels='all', I_nuclei=[0, 1/2], isotope=226,
    round=8, P_values=[1/2])
g.eigensystem(0, 1e-3); e.eigensystem(0, 1e-3)
Bz = np.linspace(1e-3, 100., 11)
g.ZeemanMap(Bz, plot=False); e.ZeemanMap(Bz, plot=False)

# OLD reference: per-(i,p) Calculate_TDM_evecs (same as cell 7 before patch)
def build_old(g, e):
    Bz = g.Bz; N_B = len(Bz)
    p_lab = (-1, 0, +1); q_electronic = (-1, +1)
    TDM_full = np.zeros((N_B, 3, e.size, g.size), dtype=complex)
    with contextlib.redirect_stdout(io.StringIO()):
        for ip, p in enumerate(p_lab):
            for i in range(N_B):
                TDM_full[i, ip] = Calculate_TDM_evecs(
                    p, g.evecs_B[i], g, e.evecs_B[i], e,
                    q=list(q_electronic))
    return TDM_full

# NEW: import the patched build_transition_map from the notebook by exec
import json
nb = json.loads(pathlib.Path("Jupyter Notebooks/RaX/zeeman_slower_2.ipynb").read_text())
src_imports = "".join(nb["cells"][2]["source"])  # gets numpy etc.
src_fn      = "".join(nb["cells"][7]["source"])
ns = {}
buf = io.StringIO()
with contextlib.redirect_stdout(buf):
    exec(src_imports, ns)
    exec(src_fn, ns)
    TM = ns["build_transition_map"](g, e,
                                    ground_query={'N': [1]},
                                    excited_query={'J': 0.5})
captured = buf.getvalue()
n_prints = captured.count("Successfully converted eigenvectors")

TDM_old = build_old(g, e)
TDM_new = TM["TDM_full"]
diff = np.abs(TDM_new - TDM_old).max()

print(f"max |TDM_new - TDM_old| = {diff:.3e}")
print(f"prints captured during build_transition_map = {n_prints}")
assert diff == 0.0, f"NOT bit-exact: {diff}"
assert n_prints == 0, f"still spamming: {n_prints} prints"
print("BIT-EXACT MATCH; 0 prints captured")
```

Run BEFORE patching: `conda run -n Structure python /tmp/task1_verify.py`
Expected output: assertion failure on `n_prints` (>0), since old code spams.

- [ ] **Step 2: Patch cell 7**

Save as `/tmp/task1_patch.py`:

```python
import json, pathlib
p = pathlib.Path("Jupyter Notebooks/RaX/zeeman_slower_2.ipynb")
nb = json.loads(p.read_text())
def src_to_list(s):
    lines = s.split("\n")
    return [l + "\n" for l in lines[:-1]] + ([lines[-1]] if lines[-1] else [])
for c in nb["cells"]:
    if c["cell_type"] == "code":
        c["outputs"] = []; c["execution_count"] = None

nb["cells"][7]["source"] = src_to_list(
    "import io as _io, contextlib as _ctx\n"
    "\n"
    "def _hoist_to_aBJ(state):\n"
    "    \"\"\"Convert state.evecs_B (all B-points) to aBJ once. Suppresses\n"
    "    convert_evecs's verbose=True default.  Returns (evecs_aBJ, q_numbers_aBJ).\"\"\"\n"
    "    if 'a' in state.hunds_case:\n"
    "        return state.evecs_B, state.q_numbers\n"
    "    with _ctx.redirect_stdout(_io.StringIO()):\n"
    "        probe = state.convert_evecs('aBJ', evecs=state.evecs_B[0], Normalize=False)\n"
    "        out = np.empty((len(state.Bz),) + probe.shape, dtype=probe.dtype)\n"
    "        out[0] = probe\n"
    "        for i in range(1, len(state.Bz)):\n"
    "            out[i] = state.convert_evecs('aBJ', evecs=state.evecs_B[i], Normalize=False)\n"
    "    return out, state.alt_q_numbers['aBJ']\n"
    "\n"
    "\n"
    "def build_transition_map(g, e, ground_query=None, excited_query=None,\n"
    "                          parity_g=None, parity_e=None,\n"
    "                          q_electronic=(-1, +1)):\n"
    "    '''Return dict of arrays for every (i, j, p) over Bz.\n"
    "\n"
    "    Shapes: TDM[N_B, 3, n_e, n_g], dE[N_B, n_e, n_g], slope[N_B, n_e, n_g].\n"
    "    p index: [0,1,2] -> [sigma-, pi, sigma+]  (lab spherical p = -1, 0, +1).\n"
    "    Case conversion is hoisted outside the (B, p) loop; the inner work is a\n"
    "    pure batched matmul against TDM_p_builders[iso_state].\n"
    "    '''\n"
    "    Bz = g.Bz\n"
    "    assert np.array_equal(Bz, e.Bz)\n"
    "    N_B = len(Bz)\n"
    "\n"
    "    g_idx = (np.arange(g.size) if ground_query  is None\n"
    "             else g.select_q(ground_query,  parity=parity_g))\n"
    "    e_idx = (np.arange(e.size) if excited_query is None\n"
    "             else e.select_q(excited_query, parity=parity_e))\n"
    "\n"
    "    G_a, G_qn = _hoist_to_aBJ(g)\n"
    "    E_a, E_qn = _hoist_to_aBJ(e)\n"
    "\n"
    "    p_lab = (-1, 0, +1)\n"
    "    TDM_full = np.zeros((N_B, 3, e.size, g.size), dtype=complex)\n"
    "    for ip, p in enumerate(p_lab):\n"
    "        TDM_matrix = e.library.TDM_p_builders[e.iso_state](\n"
    "            p, list(q_electronic), G_qn, E_qn)\n"
    "        TDM_full[:, ip] = E_a @ TDM_matrix @ G_a.transpose(0, 2, 1)\n"
    "\n"
    "    TDM = TDM_full[:, :, e_idx[:, None], g_idx[None, :]]   # (N_B, 3, n_e, n_g)\n"
    "    dE  = e.evals_B[:, e_idx, None] - g.evals_B[:, None, g_idx]  # MHz\n"
    "    slope = np.gradient(dE, Bz, axis=0)                          # MHz / G\n"
    "\n"
    "    return {\n"
    "        'Bz': Bz, 'g_idx': g_idx, 'e_idx': e_idx,\n"
    "        'TDM': TDM, 'TDM_full': TDM_full,\n"
    "        'dE': dE, 'slope': slope,\n"
    "        'p_lab': p_lab, 'q_electronic': tuple(q_electronic),\n"
    "    }\n"
    "\n"
    "# ALL of N=1 ground x A J'=1/2 (both parities of the Lambda-doublet).\n"
    "TM_full = build_transition_map(g, e,\n"
    "                               ground_query={'N': [1]},\n"
    "                               excited_query={'J': 0.5})\n"
    "print({k: (v.shape if hasattr(v,'shape') else type(v).__name__) for k,v in TM_full.items()})"
)
p.write_text(json.dumps(nb, indent=1) + "\n")
print("Patched cell 7")
print("--- preview ---")
print("".join(nb["cells"][7]["source"])[:800])
```

Run: `python3 /tmp/task1_patch.py`
Expected: prints "Patched cell 7" + preview showing `_hoist_to_aBJ` and the new `build_transition_map`.

- [ ] **Step 3: Run verify; assert it now passes**

Run: `conda run -n Structure python /tmp/task1_verify.py`
Expected:
```
max |TDM_new - TDM_old| = 0.000e+00
prints captured during build_transition_map = 0
BIT-EXACT MATCH; 0 prints captured
```

- [ ] **Step 4: Commit**

```bash
git add "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"
git commit -m "zeeman_slower_2: hoist case-conversion in build_transition_map (3x faster, no print spam)"
```

---

## Task 2: Add leakage diagnostic to `steady_state_rate_eq`

**Goal:** Print/warn the BR leakage at the queried B inside `steady_state_rate_eq` so the user sees when the closed-cycle assumption is over-stating the force. No change to the math — diagnostic only.

**Files:**
- Modify: `Jupyter Notebooks/RaX/zeeman_slower_2.ipynb` cell 27 (`steady_state_rate_eq` + smoke test below it)

**Acceptance Criteria:**
- [ ] `steady_state_rate_eq` accepts a new optional kwarg `verbose_leakage=True` (default `True`).
- [ ] When `verbose_leakage=True` and leakage at the queried B exceeds 1% (`leakage_thresh=0.01`), the function prints a one-line warning: `"⚠ closed-cycle assumption: subset leakage = X% at B = Y G"`.
- [ ] When leakage ≤ threshold, no warning printed.
- [ ] When `verbose_leakage=False`, never printed regardless of leakage.
- [ ] Returned dict gains a new key `'leakage'` containing the per-excited-level leakage vector (shape `(n_e_sub,)`) at the queried B.
- [ ] Existing return keys / shapes unchanged.
- [ ] Smoke-test cell 27 still passes (`Σ N = 1`, finite force).

**Verify:** Run `/tmp/task2_verify.py` → smoke test produces the expected printed warning at a known-leaky B and stays silent at a known-clean B.

**Steps:**

- [ ] **Step 1: Write the failing test**

Save as `/tmp/task2_verify.py`:

```python
import sys, io, contextlib, json, pathlib
sys.path.insert(0, "Source Code")
from config_path import add_to_sys_path; add_to_sys_path()
import numpy as np
from Energy_Levels import MoleculeLevels

# Re-create state + transition map identical to notebook
g = MoleculeLevels.initialize_state(
    molecule_name='RaF', elec_state='X', vib_state=0,
    N_list=[1], fermion_or_boson='boson',
    M_sublevels='all', I_nuclei=[0, 1/2], isotope=226,
    round=8, P_values=[1/2])
e = MoleculeLevels.initialize_state(
    molecule_name='RaF', elec_state='A', vib_state=0,
    N_list=[1, 2], fermion_or_boson='boson',
    M_sublevels='all', I_nuclei=[0, 1/2], isotope=226,
    round=8, P_values=[1/2])
g.eigensystem(0, 1e-3); e.eigensystem(0, 1e-3)
Bz = np.linspace(1e-3, 1000., 21)
g.ZeemanMap(Bz, plot=False); e.ZeemanMap(Bz, plot=False)

nb = json.loads(pathlib.Path("Jupyter Notebooks/RaX/zeeman_slower_2.ipynb").read_text())
ns = {}
with contextlib.redirect_stdout(io.StringIO()):
    exec("".join(nb["cells"][2]["source"]), ns)   # constants
    exec("".join(nb["cells"][7]["source"]), {**ns, "g": g, "e": e}, ns)  # build_transition_map
    exec("".join(nb["cells"][8]["source"]), ns)   # TDM2 normalization
    exec("".join(nb["cells"][12]["source"]), {**ns, "g": g, "e": e}, ns)  # subset_TM
    exec("".join(nb["cells"][22]["source"]), ns)   # compute_branching
    exec("".join(nb["cells"][27]["source"]), ns)   # steady_state_rate_eq + smoke

# Check the new return key + warning behavior
TM = ns["TM"]
ssr = ns["steady_state_rate_eq"]

# Find a B where leakage is small
leakage = ns["leakage"]   # from compute_branching (N_B, n_e_sub)
iB_clean = int(np.argmin(leakage.mean(axis=1)))
iB_leaky = int(np.argmax(leakage.mean(axis=1)))
B_clean = TM["Bz"][iB_clean]; B_leaky = TM["Bz"][iB_leaky]

# Clean B: should NOT print warning if leakage < 1%
buf = io.StringIO()
with contextlib.redirect_stdout(buf):
    res_clean = ssr(TM, B_clean, 0., polarization='sigma+', s_per_p=2., remix=1e-3)
clean_out = buf.getvalue()
print(f"At B={B_clean:.1f} G: leakage(mean)={leakage[iB_clean].mean():.3e}; "
      f"warning printed = {'closed-cycle' in clean_out}")

# Leaky B (or fall back to forcing one): should print warning if > 1%
buf = io.StringIO()
with contextlib.redirect_stdout(buf):
    res_leaky = ssr(TM, B_leaky, 0., polarization='sigma+', s_per_p=2., remix=1e-3)
leaky_out = buf.getvalue()
print(f"At B={B_leaky:.1f} G: leakage(mean)={leakage[iB_leaky].mean():.3e}; "
      f"warning printed = {'closed-cycle' in leaky_out}")

# verbose_leakage=False suppresses
buf = io.StringIO()
with contextlib.redirect_stdout(buf):
    ssr(TM, B_leaky, 0., polarization='sigma+', s_per_p=2., remix=1e-3,
        verbose_leakage=False)
silent_out = buf.getvalue()
print(f"verbose_leakage=False: no warning = {'closed-cycle' not in silent_out}")

# 'leakage' returned in dict
assert 'leakage' in res_clean, "missing 'leakage' key in returned dict"
assert res_clean['leakage'].shape == (TM['TDM2_norm'].shape[2],), "wrong leakage shape"
print(f"returned leakage shape OK: {res_clean['leakage'].shape}")
print("Σ N = 1 check:", np.isclose(res_clean['N'].sum(), 1.0))
```

Run BEFORE patching: `conda run -n Structure python /tmp/task2_verify.py`
Expected: `KeyError: 'leakage'` or similar — the new key/behavior doesn't exist yet.

- [ ] **Step 2: Patch cell 27**

Save as `/tmp/task2_patch.py` — use the boilerplate template above with this targeted replacement on cell 27. The patch adds two pieces: (a) accept and use `verbose_leakage` and `leakage_thresh` kwargs, (b) compute and store/print leakage at the queried B.

Replace the function signature line:

```
old:
def steady_state_rate_eq(TM, Bz_query, v, polarization='sigma+', s_per_p=1.0,
                          detuning_MHz=0.0, ref_transition=None,
                          Gamma_2pi_MHz=Gamma_2pi_MHz, remix=1e-3):

new:
def steady_state_rate_eq(TM, Bz_query, v, polarization='sigma+', s_per_p=1.0,
                          detuning_MHz=0.0, ref_transition=None,
                          Gamma_2pi_MHz=Gamma_2pi_MHz, remix=1e-3,
                          verbose_leakage=True, leakage_thresh=0.01):
```

Replace the BR computation block (right before `Gamma_total_s = 2*np.pi*Gamma_2pi_MHz*1e6`) with a leakage-aware version. The current block:

```
old:
    # spontaneous-emission BR within subset (closed-cycle assumption)
    weights_e2g = TDM2[iB].sum(axis=0)                          # (n_e, n_g)
    Z = weights_e2g.sum(axis=1, keepdims=True)
    BR_e2g = weights_e2g / np.where(Z > 0, Z, 1.0)
    Gamma_total_s = 2*np.pi*Gamma_2pi_MHz*1e6                   # 1/s

new:
    # spontaneous-emission BR within subset (closed-cycle assumption)
    weights_e2g = TDM2[iB].sum(axis=0)                          # (n_e, n_g) — subset
    weights_full_e2g = (np.abs(TM['TDM_full'][iB])**2).sum(axis=0)[
        TM['e_idx'][:, None], np.arange(TM['TDM_full'].shape[3])[None, :]
    ]                                                            # (n_e_sub, n_g_full)
    Z_subset = weights_e2g.sum(axis=1)                          # (n_e,)
    Z_full   = weights_full_e2g.sum(axis=1)                     # (n_e,)
    leakage_at_B = 1.0 - Z_subset / np.where(Z_full > 0, Z_full, 1.0)  # (n_e,)
    if verbose_leakage and leakage_at_B.max() > leakage_thresh:
        print(f"⚠ closed-cycle assumption: subset leakage = "
              f"{100*leakage_at_B.max():.2f}% at B = {Bz[iB]:.1f} G")
    BR_e2g = weights_e2g / np.where(Z_subset[:, None] > 0,
                                    Z_subset[:, None], 1.0)
    Gamma_total_s = 2*np.pi*Gamma_2pi_MHz*1e6                   # 1/s
```

Replace the return dict to include `'leakage'`:

```
old:
    return {
        'N': N, 'Nl': Nl, 'Nu': Nu,
        'R_per_pol': R_p_force, 'BR_e2g': BR_e2g,
        'F': F, 'scattering_rate': scattering_rate,
        'Bz_used': Bz[iB], 'omega_L_MHz': omega_L_MHz,
        'iB': iB, 'ref': (je0, ig0),
    }

new:
    return {
        'N': N, 'Nl': Nl, 'Nu': Nu,
        'R_per_pol': R_p_force, 'BR_e2g': BR_e2g,
        'F': F, 'scattering_rate': scattering_rate,
        'Bz_used': Bz[iB], 'omega_L_MHz': omega_L_MHz,
        'iB': iB, 'ref': (je0, ig0),
        'leakage': leakage_at_B,
    }
```

Patch script (full):

```python
import json, pathlib
p = pathlib.Path("Jupyter Notebooks/RaX/zeeman_slower_2.ipynb")
nb = json.loads(p.read_text())
def src_to_list(s):
    lines = s.split("\n")
    return [l + "\n" for l in lines[:-1]] + ([lines[-1]] if lines[-1] else [])
for c in nb["cells"]:
    if c["cell_type"] == "code":
        c["outputs"] = []; c["execution_count"] = None

src = "".join(nb["cells"][27]["source"])

# Replacement 1: signature
old_sig = ("def steady_state_rate_eq(TM, Bz_query, v, polarization='sigma+', s_per_p=1.0,\n"
           "                          detuning_MHz=0.0, ref_transition=None,\n"
           "                          Gamma_2pi_MHz=Gamma_2pi_MHz, remix=1e-3):\n")
new_sig = ("def steady_state_rate_eq(TM, Bz_query, v, polarization='sigma+', s_per_p=1.0,\n"
           "                          detuning_MHz=0.0, ref_transition=None,\n"
           "                          Gamma_2pi_MHz=Gamma_2pi_MHz, remix=1e-3,\n"
           "                          verbose_leakage=True, leakage_thresh=0.01):\n")
assert old_sig in src, "signature replacement failed — check old_sig"
src = src.replace(old_sig, new_sig)

# Replacement 2: BR block
old_br = ("    # spontaneous-emission BR within subset (closed-cycle assumption)\n"
          "    weights_e2g = TDM2[iB].sum(axis=0)                          # (n_e, n_g)\n"
          "    Z = weights_e2g.sum(axis=1, keepdims=True)\n"
          "    BR_e2g = weights_e2g / np.where(Z > 0, Z, 1.0)\n"
          "    Gamma_total_s = 2*np.pi*Gamma_2pi_MHz*1e6                   # 1/s\n")
new_br = ("    # spontaneous-emission BR within subset (closed-cycle assumption)\n"
          "    weights_e2g = TDM2[iB].sum(axis=0)                          # (n_e, n_g) — subset\n"
          "    e_idx_full = TM['e_idx']\n"
          "    weights_full_e2g = (np.abs(TM['TDM_full'][iB])**2).sum(axis=0)[e_idx_full, :]  # (n_e, n_g_full)\n"
          "    Z_subset = weights_e2g.sum(axis=1)                          # (n_e,)\n"
          "    Z_full   = weights_full_e2g.sum(axis=1)                     # (n_e,)\n"
          "    leakage_at_B = 1.0 - Z_subset / np.where(Z_full > 0, Z_full, 1.0)\n"
          "    if verbose_leakage and leakage_at_B.max() > leakage_thresh:\n"
          "        print(f\"⚠ closed-cycle assumption: subset leakage = \"\n"
          "              f\"{100*leakage_at_B.max():.2f}% at B = {Bz[iB]:.1f} G\")\n"
          "    BR_e2g = weights_e2g / np.where(Z_subset[:, None] > 0, Z_subset[:, None], 1.0)\n"
          "    Gamma_total_s = 2*np.pi*Gamma_2pi_MHz*1e6                   # 1/s\n")
assert old_br in src, "BR-block replacement failed"
src = src.replace(old_br, new_br)

# Replacement 3: return dict
old_ret = ("    return {\n"
           "        'N': N, 'Nl': Nl, 'Nu': Nu,\n"
           "        'R_per_pol': R_p_force, 'BR_e2g': BR_e2g,\n"
           "        'F': F, 'scattering_rate': scattering_rate,\n"
           "        'Bz_used': Bz[iB], 'omega_L_MHz': omega_L_MHz,\n"
           "        'iB': iB, 'ref': (je0, ig0),\n"
           "    }\n")
new_ret = ("    return {\n"
           "        'N': N, 'Nl': Nl, 'Nu': Nu,\n"
           "        'R_per_pol': R_p_force, 'BR_e2g': BR_e2g,\n"
           "        'F': F, 'scattering_rate': scattering_rate,\n"
           "        'Bz_used': Bz[iB], 'omega_L_MHz': omega_L_MHz,\n"
           "        'iB': iB, 'ref': (je0, ig0),\n"
           "        'leakage': leakage_at_B,\n"
           "    }\n")
assert old_ret in src, "return-dict replacement failed"
src = src.replace(old_ret, new_ret)

nb["cells"][27]["source"] = src_to_list(src)
p.write_text(json.dumps(nb, indent=1) + "\n")
print("Patched cell 27 (3 replacements)")
```

Run: `python3 /tmp/task2_patch.py`
Expected: "Patched cell 27 (3 replacements)" — no AssertionError.

- [ ] **Step 3: Re-run verify, assert it now passes**

Run: `conda run -n Structure python /tmp/task2_verify.py`
Expected: prints leakage values for clean and leaky B, asserts the dict key is present and Σ N = 1.

- [ ] **Step 4: Commit**

```bash
git add "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"
git commit -m "zeeman_slower_2: leakage diagnostic in steady_state_rate_eq (warn if subset closure broken)"
```

---

## Task 3: Extend Doppler-trackable cluster criterion (slope AND frequency offset)

**Goal:** Make `find_doppler_trackable_bands` cluster on slope-vs-B AND mean frequency offset, so two transitions far apart in ν don't get falsely merged. Update the markdown above the cell so the new criterion matches the wording.

**Files:**
- Modify: `Jupyter Notebooks/RaX/zeeman_slower_2.ipynb` cell 16 (markdown) and cell 17 (`find_doppler_trackable_bands`)

**Acceptance Criteria:**
- [ ] `find_doppler_trackable_bands` accepts new kwarg `freq_tol_MHz=200.` (default).
- [ ] Two transitions with similar slope (within `tol`) are NOT merged if `|<ν_i> − <ν_j>|_avg_over_B > freq_tol_MHz`.
- [ ] Cluster count when `freq_tol_MHz=np.inf` is identical to the pre-patch behavior (sanity check — old behavior recoverable).
- [ ] Markdown cell 16 updated to reflect the joint slope+frequency criterion (one short sentence; preserve the "Doppler-trackable bands" heading).

**Verify:** Run `/tmp/task3_verify.py` → asserts (a) `freq_tol_MHz=np.inf` reproduces old cluster count, (b) `freq_tol_MHz=200` produces ≥ as many clusters (splits, never merges).

**Steps:**

- [ ] **Step 1: Write the failing test**

Save as `/tmp/task3_verify.py`:

```python
import sys, io, contextlib, json, pathlib
sys.path.insert(0, "Source Code")
from config_path import add_to_sys_path; add_to_sys_path()
import numpy as np
from Energy_Levels import MoleculeLevels

g = MoleculeLevels.initialize_state(
    molecule_name='RaF', elec_state='X', vib_state=0,
    N_list=[1], fermion_or_boson='boson',
    M_sublevels='all', I_nuclei=[0, 1/2], isotope=226,
    round=8, P_values=[1/2])
e = MoleculeLevels.initialize_state(
    molecule_name='RaF', elec_state='A', vib_state=0,
    N_list=[1, 2], fermion_or_boson='boson',
    M_sublevels='all', I_nuclei=[0, 1/2], isotope=226,
    round=8, P_values=[1/2])
g.eigensystem(0, 1e-3); e.eigensystem(0, 1e-3)
Bz = np.linspace(1e-3, 1000., 21)
g.ZeemanMap(Bz, plot=False); e.ZeemanMap(Bz, plot=False)

nb = json.loads(pathlib.Path("Jupyter Notebooks/RaX/zeeman_slower_2.ipynb").read_text())
ns = {}
with contextlib.redirect_stdout(io.StringIO()):
    exec("".join(nb["cells"][2]["source"]), ns)
    exec("".join(nb["cells"][7]["source"]), {**ns, "g": g, "e": e}, ns)
    exec("".join(nb["cells"][8]["source"]), ns)
    exec("".join(nb["cells"][12]["source"]), {**ns, "g": g, "e": e}, ns)
    exec("".join(nb["cells"][17]["source"]), ns)

TM = ns["TM"]; find = ns["find_doppler_trackable_bands"]
clusters_inf, _ = find(TM, tol=0.2, min_TDM2=0.05, freq_tol_MHz=np.inf)
clusters_200, _ = find(TM, tol=0.2, min_TDM2=0.05, freq_tol_MHz=200.)
print(f"clusters with freq_tol=inf: {len(clusters_inf)}")
print(f"clusters with freq_tol=200 MHz: {len(clusters_200)}")
assert len(clusters_200) >= len(clusters_inf), "freq_tol should split, not merge"
# Also: with infinite tol, clusters should match old behavior.  We rebuild old:
def find_old(TM, tol=0.2, min_TDM2=0.05, min_overlap=0.4):
    sl = TM["slope"]; TDM2 = TM["TDM2_norm"]
    N_B, _, n_e, n_g = TDM2.shape
    bright = (TDM2.max(axis=(0,1)) > min_TDM2)
    j_e, j_g = np.where(bright)
    keys=[]; curves=[]; weights=[]
    for je, ig in zip(j_e, j_g):
        for ip in range(3):
            w = TDM2[:, ip, je, ig]
            if w.max() < min_TDM2: continue
            keys.append((int(je), int(ig), int(ip)))
            curves.append(sl[:, je, ig]); weights.append(w)
    curves=np.array(curves); weights=np.array(weights)
    clusters=[]; assigned=np.zeros(len(curves), dtype=bool)
    for i in range(len(curves)):
        if assigned[i]: continue
        members=[i]; assigned[i]=True
        for j in range(i+1, len(curves)):
            if assigned[j]: continue
            mask=(weights[i]>min_TDM2)&(weights[j]>min_TDM2)
            if mask.sum() < min_overlap*N_B: continue
            denom=np.maximum(np.abs(curves[i][mask]), np.abs(curves[j][mask]))
            denom=np.where(denom<1e-9,1e-9,denom)
            rel=np.abs(curves[i][mask]-curves[j][mask])/denom
            if np.median(rel) <= tol:
                members.append(j); assigned[j]=True
        clusters.append([keys[m] for m in members])
    return clusters
clusters_old = find_old(TM, tol=0.2, min_TDM2=0.05)
assert len(clusters_inf) == len(clusters_old), \
    f"freq_tol=inf must reproduce old behavior: new={len(clusters_inf)}, old={len(clusters_old)}"
print("OK: freq_tol=inf reproduces old; finite freq_tol can only split")
```

Run BEFORE patching: assertion fail on `freq_tol_MHz` keyword (TypeError).

- [ ] **Step 2: Patch cell 17 (function only — leave `plot_B_slopes` and the call alone)**

Save as `/tmp/task3_patch.py`. Replace the `find_doppler_trackable_bands` function and update markdown cell 16.

```python
import json, pathlib
p = pathlib.Path("Jupyter Notebooks/RaX/zeeman_slower_2.ipynb")
nb = json.loads(p.read_text())
def src_to_list(s):
    lines = s.split("\n")
    return [l + "\n" for l in lines[:-1]] + ([lines[-1]] if lines[-1] else [])
for c in nb["cells"]:
    if c["cell_type"] == "code":
        c["outputs"] = []; c["execution_count"] = None

src17 = "".join(nb["cells"][17]["source"])
old_fn = ('def find_doppler_trackable_bands(TM, tol=0.2, min_TDM2=0.05, min_overlap=0.4):\n'
          '    """Cluster (j,i,p) transitions by similar slope-vs-B curves.\n'
          '\n'
          '    Returns list of clusters: each is a list of (je,ig,ip) tuples sharing slope\n'
          '    within tol over at least min_overlap of the B grid (where they are bright).\n'
          '    """\n')
new_fn = ('def find_doppler_trackable_bands(TM, tol=0.2, min_TDM2=0.05,\n'
          '                                   min_overlap=0.4, freq_tol_MHz=200.):\n'
          '    """Cluster (j,i,p) transitions sharing slope-vs-B AND mean frequency.\n'
          '\n'
          '    Two transitions cluster only if (a) median relative slope difference is\n'
          '    within `tol` over at least `min_overlap` of the bright-overlap B-window,\n'
          '    AND (b) their mean frequency over that window is within `freq_tol_MHz`.\n'
          '    Set freq_tol_MHz=np.inf to recover the slope-only behavior.\n'
          '    """\n')
assert old_fn in src17, "function-header replacement failed"
src17 = src17.replace(old_fn, new_fn)

# Insert the freq-coincidence test after the existing rel-slope test.
old_loop = ('            denom = np.maximum(np.abs(curves[i][mask]), np.abs(curves[j][mask]))\n'
            '            denom = np.where(denom < 1e-9, 1e-9, denom)\n'
            '            rel   = np.abs(curves[i][mask] - curves[j][mask]) / denom\n'
            '            if np.median(rel) <= tol:\n'
            '                members.append(j); assigned[j] = True\n')
new_loop = ('            denom = np.maximum(np.abs(curves[i][mask]), np.abs(curves[j][mask]))\n'
            '            denom = np.where(denom < 1e-9, 1e-9, denom)\n'
            '            rel   = np.abs(curves[i][mask] - curves[j][mask]) / denom\n'
            '            if np.median(rel) > tol: continue\n'
            '            # frequency-coincidence test\n'
            '            je_i, ig_i, _ = keys[i]; je_j, ig_j, _ = keys[j]\n'
            '            df = np.abs(TM["dE"][mask, je_i, ig_i].mean()\n'
            '                      - TM["dE"][mask, je_j, ig_j].mean())\n'
            '            if df > freq_tol_MHz: continue\n'
            '            members.append(j); assigned[j] = True\n')
assert old_loop in src17, "inner-loop replacement failed"
src17 = src17.replace(old_loop, new_loop)

nb["cells"][17]["source"] = src_to_list(src17)

# Markdown cell 16
nb["cells"][16]["source"] = src_to_list(
    "### Plot B — slope $d\\nu_{ij}/dB$ vs $B$, with Doppler-trackable bands\n"
    "\n"
    "Cluster `(j,i,p)` transitions that share BOTH `dν/dB` (within `tol`) AND mean\n"
    "frequency (within `freq_tol_MHz`) over the B-window where both are bright. The\n"
    "joint criterion is what \"Doppler-trackable\" actually requires — slope agreement\n"
    "alone is necessary but not sufficient (two ~zero-offset transitions at very\n"
    "different frequencies need different lasers, even if their B-tracking is\n"
    "identical). Zero-slope transitions are flagged as untrackable for any laser."
)

p.write_text(json.dumps(nb, indent=1) + "\n")
print("Patched cells 16 (md) and 17 (fn)")
```

Run: `python3 /tmp/task3_patch.py`
Expected: "Patched cells 16 (md) and 17 (fn)".

- [ ] **Step 3: Re-run verify, assert pass**

Run: `conda run -n Structure python /tmp/task3_verify.py`
Expected: prints both cluster counts, asserts `freq_tol=200` ≥ `freq_tol=inf`, and `freq_tol=inf` reproduces hand-rolled old behavior.

- [ ] **Step 4: Commit**

```bash
git add "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"
git commit -m "zeeman_slower_2: cluster by slope AND frequency in find_doppler_trackable_bands"
```

---

## Task 4: PB shading from molecule_parameters constants

**Goal:** Replace hardcoded `B_SR_PB=(50, 200)` and `B_HF_PB=(500, 2500)` defaults in `plot_A_transition_frequencies` with values derived from `Gamma_SR` and `bF` in `Source Code/molecule_parameters.py` (read at notebook load time, isotope-aware).

**Files:**
- Modify: `Jupyter Notebooks/RaX/zeeman_slower_2.ipynb` cell 15 (`plot_A_transition_frequencies` defaults + a new cell-2-style read of params)

**Acceptance Criteria:**
- [ ] A new line in cell 15 (or in a tiny prelude inside the function) reads `Gamma_SR` and `bF` from `molecule_parameters.molecules[name][stat][state]`.
- [ ] Default `B_SR_PB = (0.3 * Gamma_SR/mu_B, 3.0 * Gamma_SR/mu_B)`, default `B_HF_PB = (0.3 * bF/mu_B, 3.0 * bF/mu_B)`. (The 0.3-3 factors capture the PB transition width — onset at ratio ~1, decoupled by ~3.)
- [ ] User-supplied `B_SR_PB` / `B_HF_PB` arguments still take precedence (defaults are only the fallback).
- [ ] For ²²⁶RaF X²Σ N=1, the new defaults compute to `B_SR_PB ≈ (38, 376) G` and `B_HF_PB ≈ (21, 207) G` — sanity-printed by the verify script.

**Verify:** Run `/tmp/task4_verify.py` → asserts numerical defaults match the analytical formula for ²²⁶RaF.

**Steps:**

- [ ] **Step 1: Write the failing test**

Save as `/tmp/task4_verify.py`:

```python
import sys, io, contextlib, json, pathlib
sys.path.insert(0, "Source Code")
from config_path import add_to_sys_path; add_to_sys_path()
import numpy as np
import molecule_parameters as mp

mu_B = 1.399624494
Gamma_SR = mp.molecules['RaF']['boson']['X0']['Gamma_SR']   # 175.38
bF       = mp.molecules['RaF']['boson']['X0']['bF']         # 96.3
print(f"Gamma_SR = {Gamma_SR}, bF = {bF}")
B_SR_expect = (0.3*Gamma_SR/mu_B, 3.0*Gamma_SR/mu_B)
B_HF_expect = (0.3*bF/mu_B,       3.0*bF/mu_B)
print(f"expected B_SR_PB = ({B_SR_expect[0]:.1f}, {B_SR_expect[1]:.1f})")
print(f"expected B_HF_PB = ({B_HF_expect[0]:.1f}, {B_HF_expect[1]:.1f})")

nb = json.loads(pathlib.Path("Jupyter Notebooks/RaX/zeeman_slower_2.ipynb").read_text())
ns = {}
with contextlib.redirect_stdout(io.StringIO()):
    exec("".join(nb["cells"][2]["source"]), ns)
    exec("".join(nb["cells"][15]["source"]), ns)

import inspect
sig = inspect.signature(ns["plot_A_transition_frequencies"])
B_SR_default = sig.parameters["B_SR_PB"].default
B_HF_default = sig.parameters["B_HF_PB"].default
print(f"actual   B_SR_PB = {B_SR_default}")
print(f"actual   B_HF_PB = {B_HF_default}")
assert np.allclose(B_SR_default, B_SR_expect, rtol=1e-3), "B_SR_PB default off"
assert np.allclose(B_HF_default, B_HF_expect, rtol=1e-3), "B_HF_PB default off"
print("OK: PB shading defaults derived from molecule_parameters")
```

Run BEFORE patching: assertion failure (defaults still 50/200, 500/2500).

- [ ] **Step 2: Patch cell 15**

Save as `/tmp/task4_patch.py`:

```python
import json, pathlib
p = pathlib.Path("Jupyter Notebooks/RaX/zeeman_slower_2.ipynb")
nb = json.loads(p.read_text())
def src_to_list(s):
    lines = s.split("\n")
    return [l + "\n" for l in lines[:-1]] + ([lines[-1]] if lines[-1] else [])
for c in nb["cells"]:
    if c["cell_type"] == "code":
        c["outputs"] = []; c["execution_count"] = None

src = "".join(nb["cells"][15]["source"])
old_sig = ("def plot_A_transition_frequencies(TM, g_char, g_pos,\n"
           "                                  B_SR_PB=(50, 200), B_HF_PB=(500, 2500),\n"
           "                                  lw_min=0.2, lw_max=2.5,\n"
           "                                  fig_kwargs=None):\n")
new_sig = ("import molecule_parameters as _mp\n"
           "_X0 = _mp.molecules['RaF']['boson']['X0']\n"
           "_B_SR_default = (0.3 * _X0['Gamma_SR'] / mu_B, 3.0 * _X0['Gamma_SR'] / mu_B)\n"
           "_B_HF_default = (0.3 * _X0['bF']       / mu_B, 3.0 * _X0['bF']       / mu_B)\n"
           "\n"
           "def plot_A_transition_frequencies(TM, g_char, g_pos,\n"
           "                                  B_SR_PB=_B_SR_default, B_HF_PB=_B_HF_default,\n"
           "                                  lw_min=0.2, lw_max=2.5,\n"
           "                                  fig_kwargs=None):\n")
assert old_sig in src, "signature replacement failed"
src = src.replace(old_sig, new_sig)
nb["cells"][15]["source"] = src_to_list(src)

p.write_text(json.dumps(nb, indent=1) + "\n")
print("Patched cell 15")
```

Run: `python3 /tmp/task4_patch.py`
Expected: "Patched cell 15".

- [ ] **Step 3: Re-run verify, assert pass**

Run: `conda run -n Structure python /tmp/task4_verify.py`
Expected: prints expected vs actual defaults, both match within 1e-3, prints "OK: PB shading defaults derived...".

- [ ] **Step 4: Commit**

```bash
git add "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"
git commit -m "zeeman_slower_2: derive PB shading defaults from molecule_parameters (Gamma_SR, bF)"
```

---

## Task 5: Cleanup pass — `a_recoil`, unused import, reference label, rate-eq caveat

**Goal:** Four small textual edits with no behavior change.

**Files:**
- Modify: `Jupyter Notebooks/RaX/zeeman_slower_2.ipynb` cells 0 (header conventions list, optional), 2, 8, 26

**Acceptance Criteria:**
- [ ] Cell 2: `a_recoil = …` line removed (variable is unused; misleadingly equal to `v_recoil`).
- [ ] Cell 2: `from Energy_Levels import …, branching_ratios` becomes `from Energy_Levels import MoleculeLevels, Calculate_TDM_evecs` (drop the unused `branching_ratios`).
- [ ] Cell 8: comment / print on `TDM2_ref` no longer claims "case-(b) stretched ME"; reworded as "max coupling at lowest sampled B (~ stretched-state ME for case-b ground)".
- [ ] Cell 26 markdown: append a one-paragraph caveat block: "**Caveat: laser auto-retunes per query.** `steady_state_rate_eq` re-anchors `omega_L_MHz` at the queried `B` so the reference transition is on resonance. F(v, B) plots therefore represent the optimum-tuned-laser envelope at each B, NOT a fixed-frequency laser sweep."

**Verify:** Run `/tmp/task5_verify.py` → grep-style asserts for each change.

**Steps:**

- [ ] **Step 1: Write the failing test**

Save as `/tmp/task5_verify.py`:

```python
import json, pathlib
nb = json.loads(pathlib.Path("Jupyter Notebooks/RaX/zeeman_slower_2.ipynb").read_text())
src2  = "".join(nb["cells"][2]["source"])
src8  = "".join(nb["cells"][8]["source"])
src26 = "".join(nb["cells"][26]["source"])

assert "a_recoil" not in src2, "cell 2 still defines a_recoil"
assert "branching_ratios" not in src2, "cell 2 still imports branching_ratios"
assert "Calculate_TDM_evecs" in src2, "cell 2 lost Calculate_TDM_evecs import"
assert "case-(b) stretched ME" not in src8 and "case-b stretched" not in src8, \
    "cell 8 still claims 'case-(b) stretched ME'"
assert "max coupling at lowest sampled B" in src8, \
    "cell 8 missing reworded reference label"
assert "Caveat: laser auto-retunes per query" in src26, \
    "cell 26 missing rate-eq caveat block"
print("OK: cleanup pass — all 4 assertions passed")
```

Run BEFORE patching: 4 assertion failures.

- [ ] **Step 2: Patch cells 2, 8, 26**

Save as `/tmp/task5_patch.py`:

```python
import json, pathlib
p = pathlib.Path("Jupyter Notebooks/RaX/zeeman_slower_2.ipynb")
nb = json.loads(p.read_text())
def src_to_list(s):
    lines = s.split("\n")
    return [l + "\n" for l in lines[:-1]] + ([lines[-1]] if lines[-1] else [])
for c in nb["cells"]:
    if c["cell_type"] == "code":
        c["outputs"] = []; c["execution_count"] = None

# --- cell 2 ---
src2 = "".join(nb["cells"][2]["source"])
old_imp = "from Energy_Levels import MoleculeLevels, Calculate_TDM_evecs, branching_ratios"
new_imp = "from Energy_Levels import MoleculeLevels, Calculate_TDM_evecs"
assert old_imp in src2, "import line not found in cell 2"
src2 = src2.replace(old_imp, new_imp)

old_recoil = ("v_recoil = hbar_J * k_XA / mass_RaF        # m/s per scattered photon\n"
              "a_recoil = hbar_J * k_XA / mass_RaF        # m/s^2 / scattering rate factor\n")
new_recoil = ("v_recoil = hbar_J * k_XA / mass_RaF        # m/s per scattered photon\n")
assert old_recoil in src2, "v_recoil/a_recoil block not found in cell 2"
src2 = src2.replace(old_recoil, new_recoil)
nb["cells"][2]["source"] = src_to_list(src2)

# --- cell 8 ---
src8 = "".join(nb["cells"][8]["source"])
old_norm = ("# case-(b) stretched-state normalization: max |TDM|^2 over (i,j,p) at B=0.\n"
            "TDM2     = np.abs(TM_full['TDM'])**2\n"
            "TDM2_ref = float(TDM2[0].max())\n")
new_norm = ("# Reference normalization: max coupling at lowest sampled B (~ stretched-state ME\n"
            "# for case-(b) ground; A^2Pi_{1/2} is case-(a)/aBJ so the equality is approximate).\n"
            "TDM2     = np.abs(TM_full['TDM'])**2\n"
            "TDM2_ref = float(TDM2[0].max())\n")
assert old_norm in src8, "TDM2_ref block not found in cell 8"
src8 = src8.replace(old_norm, new_norm)
nb["cells"][8]["source"] = src_to_list(src8)

# --- cell 26 ---
src26 = "".join(nb["cells"][26]["source"])
caveat = ("\n\n**Caveat: laser auto-retunes per query.** `steady_state_rate_eq` re-anchors\n"
          "`omega_L_MHz` at the queried `B` so the reference transition is on resonance.\n"
          "F(v, B) plots therefore represent the optimum-tuned-laser envelope at each B,\n"
          "NOT a fixed-frequency laser sweep.")
if caveat.strip() not in src26:
    src26 = src26.rstrip() + caveat + "\n"
nb["cells"][26]["source"] = src_to_list(src26)

p.write_text(json.dumps(nb, indent=1) + "\n")
print("Patched cells 2, 8, 26")
```

Run: `python3 /tmp/task5_patch.py`
Expected: "Patched cells 2, 8, 26".

- [ ] **Step 3: Re-run verify, assert pass**

Run: `python3 /tmp/task5_verify.py`
Expected: "OK: cleanup pass — all 4 assertions passed".

- [ ] **Step 4: Commit**

```bash
git add "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"
git commit -m "zeeman_slower_2: cleanup — drop a_recoil + unused import, rework TDM2_ref label, add rate-eq caveat"
```

---

## Task 6: End-to-end notebook execute + integration verify

**Goal:** Execute the full notebook on a fresh kernel and confirm (a) it runs without errors, (b) the print-spam is gone, (c) all sanity outputs (Σ N = 1, leakage warnings) appear as expected, (d) all five plots A–E + two force plots render.

**Files:**
- Re-execute: `Jupyter Notebooks/RaX/zeeman_slower_2.ipynb` (no source change — outputs only)

**Acceptance Criteria:**
- [ ] `jupyter nbconvert --to notebook --execute --inplace` exits 0.
- [ ] No "Successfully converted eigenvectors" lines anywhere in cell outputs.
- [ ] Cell 27 smoke test prints `Σ N ≈ 1` (within 1e-6).
- [ ] At least one `steady_state_rate_eq` call in cells 28/29 prints the leakage warning if any queried B has leakage > 1% (informational, not failure).
- [ ] Total runtime < 180 sec (down from the previous spam-loaded baseline).

**Verify:** Run end-to-end execute below; then run `/tmp/task6_verify.py` to scan outputs.

**Steps:**

- [ ] **Step 1: Kill any concurrent jupyter on this notebook**

```bash
pkill -f 'jupyter.*zeeman_slower_2' || true
sleep 2
ps aux | grep -i 'jupyter.*zeeman_slower_2' | grep -v grep
```
Expected: no matching processes.

- [ ] **Step 2: Execute end-to-end**

```bash
time conda run -n Structure jupyter nbconvert --to notebook --execute --inplace \
    "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb" 2>&1 | tee /tmp/task6_exec.log
```
Expected:
- Exit code 0.
- `real` time well under 180s (the helper notebook with Bz=21 ran in ~15s; this one has Bz=101 + rate-eq scans, expect 60–120s).
- No `Traceback`.

- [ ] **Step 3: Scan executed notebook for spam-free + sanity prints**

Save as `/tmp/task6_verify.py`:

```python
import json, pathlib, re
nb = json.loads(pathlib.Path("Jupyter Notebooks/RaX/zeeman_slower_2.ipynb").read_text())
all_out = []
for c in nb["cells"]:
    if c["cell_type"] != "code": continue
    for o in c.get("outputs", []):
        if o.get("output_type") == "stream":
            t = o.get("text"); 
            if isinstance(t, list): t = "".join(t)
            all_out.append(t)
all_text = "\n".join(all_out)
n_spam = all_text.count("Successfully converted eigenvectors")
print(f"'Successfully converted eigenvectors' lines: {n_spam}")
assert n_spam == 0, f"print spam still present ({n_spam} lines)"

# Sigma N == 1
m = re.search(r"sum_N\s*=\s*([0-9.eE+-]+)", all_text)
assert m, "smoke test 'sum_N = ...' not found"
sumN = float(m.group(1))
print(f"smoke-test sum_N = {sumN}")
assert abs(sumN - 1.0) < 1e-6, f"Σ N off: {sumN}"

# Leakage diagnostic shows up at least once if any leakage > 1% in the run
n_leak_warn = all_text.count("closed-cycle assumption: subset leakage")
print(f"leakage warnings printed: {n_leak_warn} (informational)")
print("OK: end-to-end clean")
```

Run: `python3 /tmp/task6_verify.py`
Expected: "OK: end-to-end clean".

- [ ] **Step 4: Commit (outputs)**

```bash
git add "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"
git commit -m "zeeman_slower_2: re-execute end-to-end (clean run, no print spam)"
```

---

## Self-review (run by author before handoff)

1. **Spec coverage** — six top fixes from the audit:
   - B2 (print spam) → Task 1 ✓
   - B4 (leakage diagnostic) → Task 2 ✓
   - B3 (cluster criterion) → Task 3 ✓
   - B5 (PB shading) → Task 4 ✓
   - B1 (a_recoil) → Task 5 ✓
   - C1 (caveat markdown) → Task 5 ✓
   - Plus C2 (label) and C4 (unused import) opportunistically folded into Task 5.
   - C3 (Bz[0]=1e-3 endpoint slope) explicitly NOT addressed — minor cosmetic.
2. **Placeholder scan** — every code/patch step contains complete code. No "TODO", "TBD", "implement later".
3. **Type consistency** — `_hoist_to_aBJ` returns `(evecs, qn)` consistently in Task 1; `leakage_at_B` shape `(n_e_sub,)` is the same in patch + verify in Task 2; `freq_tol_MHz` named consistently across Task 3 patch/verify; `B_SR_PB`/`B_HF_PB` defaults match formula in Task 4 patch/verify.

---

## Notable deviations from the audit

- The audit text claimed `γ_SR` is "sub-MHz" and `b_F` is "tens of MHz". The actual values in `molecule_parameters.py` for ²²⁶RaF X0 are `Gamma_SR = 175.38 MHz` and `bF = 96.3 MHz`. PB onsets are therefore ~125 G (SR) and ~69 G (HF). Task 4 uses the correct constants. The original `B_SR_PB=(50, 200)` is in the right ballpark; `B_HF_PB=(500, 2500)` is still way too high.
- The `_hoist_to_aBJ` helper is identical to the one already in `RaF_Zeeman_Slower.ipynb` (the helper notebook from the prior session). DRY note: a future task could lift it to a small shared module under `Jupyter Notebooks/RaX/_shared.py`, but YAGNI for two callers.
