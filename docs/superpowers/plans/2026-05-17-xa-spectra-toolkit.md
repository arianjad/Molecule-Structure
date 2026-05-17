# X–A Spectra + Averaged-BR Toolkit Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers-extended-cc:subagent-driven-development (recommended) or superpowers-extended-cc:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** One reusable, physics-correct X–A spectrum/branching module, driven off the N²→R² converter for BaF A0.

**Architecture:** New molecule-agnostic `Jupyter Notebooks/RaX/xa_spectra.py` that calls the existing `Energy_Levels.branching_ratios` / `Calculate_TDMs` (never reimplements matrix elements). BaF A0 `Origin` switches from a hand G-row to the converter (`formalism:'N2'`). Consumers (`BaF_spectrum_plot.py`, the cooling-notebook generator) import the module; the two hand-built notebooks get a snippet, not edits.

**Tech Stack:** Python, conda env `Structure`, NumPy, pandas, matplotlib, SciPy (`scipy.special.voigt_profile`). Spec: `docs/superpowers/specs/2026-05-16-xa-spectra-toolkit-design.md`.

**Conventions:** All commands run from repo root. Use `conda run -n Structure`. Never `cd` inside a shell call (Google-Drive cwd trap). Chain `git add`/`commit`/`log` in one shell invocation (GDrive index-flush trap). Commit only the named files — never `git add -A` (uncommitted user WIP in `RaF X-A.ipynb` / `BaF X-A.ipynb`). macOS git: `/Library/Developer/CommandLineTools/usr/bin/git`.

---

### Task 1: Item 3 — converter-driven BaF A0/X0 Origin

**Goal:** Replace the hand G-row in BaF A0 with `formalism:'N2'` so the converter produces the R² backend; tag X0 for a consistent convention declaration.

**Files:**
- Modify: `Source Code/molecule_parameters.py` (X0 lines 77–84, A0 lines 85–101; converter call site at line 426 is unchanged)
- Test: `Jupyter Notebooks/RaX/baf_xa_validate.py` (existing hard-assert gate, unmodified — used as the regression test)

**Acceptance Criteria:**
- [ ] `get_molecule_params('BaF','A',0,'boson')` returns `Origin` within 1e-6 cm⁻¹ of `11946.109609 + 6347.847/29979.2458 - 0.006007/29979.2458` and `formalism == 'R2'`
- [ ] A0 dict source has `'Origin': 11946.109609`, `'formalism': 'N2'`, `'Lambda': 1`; the `+ 6347.847/c` hand term is gone
- [ ] X0 dict has `'formalism': 'N2'`, `'Lambda': 0`
- [ ] `baf_xa_validate.py` prints `VALIDATION PASSED`, exit 0
- [ ] Tutorial gate exits 0 with no `--allow-errors`

**Verify:** `conda run -n Structure python "Jupyter Notebooks/RaX/baf_xa_validate.py"; echo EXIT=$?` → ends `VALIDATION PASSED`, `EXIT=0`

**Steps:**

- [ ] **Step 1: Write the failing test** — create `Jupyter Notebooks/RaX/_t1_origin.py`:

```python
import sys, os
sys.path.insert(0, os.path.join('Source Code'))
from molecule_parameters import get_molecule_params as g, params_general
c = params_general['c']
a = g('BaF', 'A', 0, 'boson')
expect = 11946.109609 + 6347.847 / c - 0.006007 / c   # T0,0 + Be_A/c - D/c (converter)
assert abs(a['Origin'] - expect) < 1e-6, (a['Origin'], expect)
assert a.get('formalism') == 'R2', a.get('formalism')   # tag follows data: N2 -> R2
assert a.get('Lambda') == 1, a.get('Lambda')
print("T1 OK", a['Origin'])
```

- [ ] **Step 2: Run it, verify it FAILS**

Run: `conda run -n Structure python "Jupyter Notebooks/RaX/_t1_origin.py"`
Expected: `AssertionError` (current hand path gives `formalism` absent and `Origin == 11946.109609 + 6347.847/c`, no `- D/c`).

- [ ] **Step 3: Edit `Source Code/molecule_parameters.py` A0 Origin line (line 100).** Replace:

```python
    'Origin': 11946.109609 + 6347.847/c,  # T0,0[cm^-1] (arXiv:2511.06986 Tbl III, "This work"; pgopher-DEFAULT N^2 fit) + Be_A[MHz]/c R^2 G-row. Code rot op is R^2-form B*(N^2-Lam^2); paper T0,0 is N^2-convention (pgopher default; paper p.8: 0.21 cm^-1~B_A offset vs Steimle). Be_A=6347.847 (THIS state's own B) NOT Be_X=6473.9588 (Be_X was the +126 MHz bug class the converter/RaF-migration removed). NOT formalism:'N2' (would mis-apply -2Lam^2D to R^2-form Be/p+2q, spec §3.4). Abs line confirmed vs Table II 348666424.4 MHz in Task 2.
    }
```

with:

```python
    'Origin': 11946.109609,  # PHYSICAL N^2 band origin T0,0 [cm^-1] (arXiv:2511.06986 Tbl III, "This work"; pgopher-DEFAULT N^2 fit). The +Be_A/c R^2 G-row is now applied by the converter (formalism:'N2' below), NOT by hand. Verified deltas vs the old hand path: Origin -0.006, Be -0.012, p+2q -0.007, ASO +0.93 MHz (on ~1.9e7) -- all << 1-2 MHz validation tols. Spec docs/superpowers/specs/2026-05-16-xa-spectra-toolkit-design.md §2; thesis App A.3.1.
    'formalism': 'N2',        # paper fit is pgopher-default N^2; converter -> engine R^2 (adds Be*Lam^2/c G-row + Table 7.2 terms)
    'Lambda': 1,              # A^2Pi: |Lambda|=1 (REQUIRED whenever 'formalism' set)
    }
```

- [ ] **Step 4: Edit the A0 `'Be'` comment (line 86)** so it no longer claims the hand convention. Replace the inline comment after `'Be': 6347.847,` (keep the value) with:

```python
    'Be': 6347.847,         # Steimle2011 [24], 0.2117414 cm^-1, held fixed in fit (paper N^2 value). Converter applies the N^2->R^2 B-row (-2*Lam^2*D = -0.012 MHz, immaterial) since this dict is now formalism:'N2'.
```

- [ ] **Step 5: Tag X0 (insert two lines after line 77 `molecules['BaF']['boson']['X0'] = {`)**, before `'Be': 6473.9588,`:

```python
molecules['BaF']['boson']['X0'] = {
    'formalism': 'N2',      # convention declaration (Sigma: convention-independent; converter no-ops at Lambda=0)
    'Lambda': 0,            # X^2Sigma+: |Lambda|=0
    'Be': 6473.9588,        # Ryzlewicz1980 [22]; 138Ba19F; 0.21594802 cm^-1
```

- [ ] **Step 6: Run the failing test, verify it PASSES**

Run: `conda run -n Structure python "Jupyter Notebooks/RaX/_t1_origin.py"`
Expected: `T1 OK 11946.32134...`

- [ ] **Step 7: Run the regression + tutorial gates**

Run: `conda run -n Structure python "Jupyter Notebooks/RaX/baf_xa_validate.py"; echo EXIT=$?`
Expected: `x_split ≈ 66.250`, `a_split ≈ 21.928`, Table II `+0.6/+0.8/+1.4`, parity closure `~1e-34`, `VALIDATION PASSED`, `EXIT=0`.

Run: `conda run -n Structure jupyter execute "Jupyter Notebooks/RaF_Calcs_Tutorial.ipynb"; echo EXIT=$?`
Expected: `EXIT=0`, no `--allow-errors`.

- [ ] **Step 8: Remove the scratch test and commit**

```bash
rm -f "Jupyter Notebooks/RaX/_t1_origin.py"
/Library/Developer/CommandLineTools/usr/bin/git add "Source Code/molecule_parameters.py" && /Library/Developer/CommandLineTools/usr/bin/git commit -m "fix(molecule_parameters): BaF A0 Origin via converter (formalism:'N2', physical T0,0); tag X0" && /Library/Developer/CommandLineTools/usr/bin/git log -1 --stat
```

---

### Task 2: `xa_spectra.py` core — `line_list` (no-M + M-resolved)

**Goal:** Create the module with `line_list`, returning a tidy DataFrame of (freq, strength, labels) using the thesis line-strength definition; no Boltzmann, no `J_adjust`.

**Files:**
- Create: `Jupyter Notebooks/RaX/xa_spectra.py`
- Create: `Jupyter Notebooks/RaX/test_xa_spectra.py`

**Acceptance Criteria:**
- [ ] no-M `line_list` strength for each (g,e) row equals the `branching_ratios` matrix entry exactly (Physics A)
- [ ] BaF cooling-line (−)/N=0 frequencies from `line_list` reproduce the three Table II values within 2 MHz (regression vs `baf_xa_validate.py`)
- [ ] Empty selection returns an empty DataFrame with the documented columns and emits a warning (no crash)

**Verify:** `conda run -n Structure python "Jupyter Notebooks/RaX/test_xa_spectra.py"; echo EXIT=$?` → `T2 OK ...`, `EXIT=0`

**Steps:**

- [ ] **Step 1: Write the failing test** — create `Jupyter Notebooks/RaX/test_xa_spectra.py`:

```python
"""Hard-assert tests for xa_spectra (exit 0 iff all pass)."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'Source Code'))
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from Energy_Levels import MoleculeLevels, branching_ratios
import xa_spectra as xs

c = 29979.2458

def build(elec, N_list):
    return MoleculeLevels.initialize_state(
        molecule_name='BaF', elec_state=elec, vib_state=0,
        N_list=np.array(N_list), fermion_or_boson='boson',
        M_sublevels='none', I_nuclei=[0, 1/2], isotope=138,
        round=8, params=None, P_values=[1/2])

g = build('X', [0, 1]); e = build('A', [1, 2])
g.eigensystem(0, 0); e.eigensystem(0, 0)

# Physics A: no-M strength == branching_ratios entry exactly
BR = branching_ratios(g, e, 0, 0)                       # (nG, nE)
gx0 = g.select_q({'N': 0, 'J': 0.5})
exm = e.select_q({'J': 0.5}, parity='-')
ll = xs.line_list(g, e, gx0, exm, origin=e.parameters['Origin'])
for _, r in ll.iterrows():
    assert abs(r['strength'] - BR[int(r['g_idx']), int(r['e_idx'])]) < 1e-12, r
# Regression: the 3 Table II (-)/N=0 components within 2 MHz
TBL2 = {(1, 0): 348666402.6, (1, 1): 348666424.4, (0, 1): 348666490.0}
hits = 0
for _, r in ll.iterrows():
    key = (int(r['g_F']), int(r['e_F']))
    if key in TBL2:
        assert abs(r['freq'] - TBL2[key]) <= 2.0, (key, r['freq'], TBL2[key])
        hits += 1
assert hits == 3, hits
# Empty selection -> empty df + warning, no crash
import warnings
with warnings.catch_warnings(record=True) as w:
    warnings.simplefilter("always")
    empty = xs.line_list(g, e, np.array([], dtype=int), exm)
    assert len(empty) == 0 and list(empty.columns) == xs.LINE_COLUMNS
    assert any(issubclass(x.category, UserWarning) for x in w)
print("T2 OK", len(ll), "rows; hits", hits)
```

- [ ] **Step 2: Run it, verify it FAILS**

Run: `conda run -n Structure python "Jupyter Notebooks/RaX/test_xa_spectra.py"`
Expected: `ModuleNotFoundError: No module named 'xa_spectra'`.

- [ ] **Step 3: Create `Jupyter Notebooks/RaX/xa_spectra.py`** with full content:

```python
"""Molecule-agnostic X-A spectrum / branching toolkit.

Generalizes the copy-pasted simulate_spectra / simulate_spectra_noM /
plot_gaussian_spectrum helpers (RaF X-A.ipynb, BaF X-A.ipynb,
BaF_spectrum_plot.py) into one physics-correct module.

Physics (thesis §3.2.4.1, pp.107-115; spec
docs/superpowers/specs/2026-05-16-xa-spectra-toolkit-design.md §3):
  line strength S = sum_{M'',M',p} |<g,i|T^1_p(d)|e,j>|^2   (thesis p.109)
  averaged BR  r~ = |<g,J''||d||e,J'>|^2 / [(2J'+1) D_ge^2]  (thesis Eq.3.4)
omega_0^3 is dropped (cancels in any branching ratio; <1e-5 across one band).
NO Boltzmann population weight, NO ad-hoc J_adjust.

Callers put "Source Code" on sys.path first (notebooks: config_path).
This module never reimplements a matrix element -- it calls
Energy_Levels.branching_ratios / Calculate_TDMs.
"""
import warnings
import numpy as np
import pandas as pd
from Energy_Levels import branching_ratios, Calculate_TDMs

C_CM = 29979.2458                       # MHz per cm^-1 (== molecule_parameters.c)
_PARITY = {1: '+', -1: '-'}
_POL = {'all': [-1, 0, 1], 'z': [0], '+': [1], '-': [-1], 'x': [-1, 1]}
LINE_COLUMNS = ['freq', 'strength',
                'g_N', 'g_J', 'g_F', 'g_par', 'e_J', 'e_F', 'e_par',
                'g_idx', 'e_idx']


def _has_M(state):
    """True iff the state was built with resolved M sublevels."""
    M = getattr(state, 'q_numbers', {}).get('M', None)
    if M is None:
        return False
    return np.unique(np.asarray(M, dtype=float)).size > 1


def _dom(state, i, key):
    """Dominant basis label of quantum `key` for eigenstate i."""
    qd = int(np.argmax(state.evecs0[i] ** 2))
    return state.q_numbers[key][qd]


def _strength_matrix(Ground, Excited, Ez, Bz, pol):
    """(nG, nE) line-strength matrix.

    no-M  : branching_ratios entry = |<g,J''|d|e,J'>|^2.
    M-res : S[g,e] = sum_{p in pol} |Calculate_TDMs(p)|^2  (thesis p.109 S,
            the M'' / M' sum is realized by clustering the degenerate
            M-components downstream in broaden()).
    """
    if not _has_M(Ground):
        return branching_ratios(Ground, Excited, Ez, Bz)        # (nG, nE)
    S = None
    for pp in _POL[pol]:
        Tp = Calculate_TDMs(pp, Ground, Excited, Ez, Bz, q=[-1, 0, 1])  # (nE,nG)
        S = np.abs(Tp) ** 2 if S is None else S + np.abs(Tp) ** 2
    return S.T                                                  # -> (nG, nE)


def line_list(Ground, Excited, g_idx, e_idx, *,
              origin=None, units='MHz', field=(0.0, 0.0),
              pol='all', weight='dipole2', thresh=1e-9):
    """Tidy DataFrame of X-A transitions between selected level sets.

    Parameters
    ----------
    Ground, Excited : MoleculeLevels (eigensystem run or run here via `field`)
    g_idx, e_idx     : index arrays (e.g. from state.select_q)
    origin           : band origin in cm^-1; default Excited.parameters['Origin']
    units            : 'MHz' (default) or 'cm-1' for the freq column
    field            : (Ez, Bz)
    pol              : 'all'|'z'|'+'|'-'|'x'  (sums |TDM|^2 over chosen p)
    weight           : 'dipole2' (bare |d|^2, in-band) or 'rate' (x nu^3)
    thresh           : drop rows with strength < thresh (absolute)

    Returns a DataFrame with columns LINE_COLUMNS. Empty selection ->
    empty DataFrame (same columns) + UserWarning.
    """
    Ez, Bz = field
    Ground.eigensystem(Ez, Bz)
    Excited.eigensystem(Ez, Bz)
    if origin is None:
        origin = Excited.parameters['Origin']
    g_idx = np.asarray(g_idx, dtype=int)
    e_idx = np.asarray(e_idx, dtype=int)
    S = _strength_matrix(Ground, Excited, Ez, Bz, pol)          # (nG, nE)

    rows = []
    for ig in g_idx:
        Eg = Ground.evals0[ig]
        for ie in e_idx:
            strength = float(S[ig, ie])
            if strength < thresh:
                continue
            nu_MHz = (Excited.evals0[ie] - Eg) + C_CM * origin
            if weight == 'rate':
                strength *= nu_MHz ** 3                          # relative; const dropped
            elif weight != 'dipole2':
                raise ValueError("weight must be 'dipole2' or 'rate'")
            freq = nu_MHz if units == 'MHz' else nu_MHz / C_CM
            rows.append({
                'freq': freq, 'strength': strength,
                'g_N': _dom(Ground, ig, 'N'), 'g_J': _dom(Ground, ig, 'J'),
                'g_F': _dom(Ground, ig, 'F'), 'g_par': _PARITY[Ground.parities[ig]],
                'e_J': _dom(Excited, ie, 'J'), 'e_F': _dom(Excited, ie, 'F'),
                'e_par': _PARITY[Excited.parities[ie]],
                'g_idx': int(ig), 'e_idx': int(ie)})
    if not rows:
        warnings.warn("line_list: no transitions above threshold for the "
                      "given selection; returning empty DataFrame.")
        return pd.DataFrame(columns=LINE_COLUMNS)
    return pd.DataFrame(rows, columns=LINE_COLUMNS).sort_values(
        'freq').reset_index(drop=True)
```

- [ ] **Step 4: Run the test, verify it PASSES**

Run: `conda run -n Structure python "Jupyter Notebooks/RaX/test_xa_spectra.py"`
Expected: `T2 OK 4 rows; hits 3` (4 = the 2×2 F_X×F_A grid; 3 match Table II, F0→F0 is E1-forbidden but still listed).

- [ ] **Step 5: Commit**

```bash
/Library/Developer/CommandLineTools/usr/bin/git add "Jupyter Notebooks/RaX/xa_spectra.py" "Jupyter Notebooks/RaX/test_xa_spectra.py" && /Library/Developer/CommandLineTools/usr/bin/git commit -m "feat(xa_spectra): line_list (no-M + M-resolved line strength); Physics-A + Table-II regression" && /Library/Developer/CommandLineTools/usr/bin/git log -1 --stat
```

---

### Task 3: `averaged_branching` (thesis Eq. 3.4)

**Goal:** Add the excited-state-averaged branching ratio: average the per-excited-sublevel normalized BR over every excited F and M_F with equal weight; sums to 1 over the ground manifold.

**Files:**
- Modify: `Jupyter Notebooks/RaX/xa_spectra.py` (append `averaged_branching`)
- Modify: `Jupyter Notebooks/RaX/test_xa_spectra.py` (append T3 block)

**Acceptance Criteria:**
- [ ] `avg_df['BR'].sum()` == 1 within 1e-9 for BaF A²Π₁/₂(J=½,+)→X N=1
- [ ] The average equals the equal-weight mean of the per-sublevel columns (`per_sublevel_df`)
- [ ] Σ over X N=0 ≤ 1e-3 (parity closure preserved through the averaging)

**Verify:** `conda run -n Structure python "Jupyter Notebooks/RaX/test_xa_spectra.py"; echo EXIT=$?` → `T3 OK ...`, `EXIT=0`

**Steps:**

- [ ] **Step 1: Append the failing test to `test_xa_spectra.py`** (before the final `print`):

```python
# ---- T3: averaged_branching == thesis Eq.3.4 generalization ----
exp = e.select_q({'J': 0.5}, parity='+')                # cooling excited state
avg_df, per_df = xs.averaged_branching(g, e, exp)
assert abs(avg_df['BR'].sum() - 1.0) < 1e-9, avg_df['BR'].sum()
# equal-weight mean of per-sublevel columns reproduces the average
mean_cols = per_df[[col for col in per_df.columns if col != 'ground']].mean(axis=1)
assert np.allclose(avg_df['BR'].values, mean_cols.values, atol=1e-12)
# parity closure: X N=0 gets ~0
n0_mask = avg_df['g_N'].astype(float) == 0
assert avg_df.loc[n0_mask, 'BR'].sum() <= 1e-3, avg_df.loc[n0_mask, 'BR'].sum()
print("T3 OK  sum(BR)=%.6f  N0=%.1e" % (avg_df['BR'].sum(),
      avg_df.loc[n0_mask, 'BR'].sum()))
```

- [ ] **Step 2: Run T3, verify it FAILS**

Run: `conda run -n Structure python "Jupyter Notebooks/RaX/test_xa_spectra.py"`
Expected: `AttributeError: module 'xa_spectra' has no attribute 'averaged_branching'`.

- [ ] **Step 3: Append `averaged_branching` to `xa_spectra.py`:**

```python
def averaged_branching(Ground, Excited, e_select, *,
                        field=(0.0, 0.0), normalize=True):
    """Branching ratio averaged over every excited sublevel in `e_select`.

    Thesis Eq.3.4 generalized to hyperfine: for each excited eigenstate
    (every F and, in an M basis, every M_F) take the BR distribution over
    the ground manifold, normalize it, then average with EQUAL weight per
    excited sublevel (the (2J'+1)^-1 / N_excited orientation average).

    Returns
    -------
    avg_df : DataFrame [g_N, g_J, g_F, g_par, BR, g_idx]  (sum BR == 1)
    per_df : DataFrame ['ground' label + one column per excited sublevel]
    """
    Ez, Bz = field
    BR = branching_ratios(Ground, Excited, Ez, Bz)              # (nG, nE)
    e_select = np.asarray(e_select, dtype=int)
    cols, labels = [], []
    for ie in e_select:
        col = BR[:, ie].astype(float)
        tot = col.sum()
        if tot <= 0:
            raise ValueError(f"excited idx {ie} has zero total decay")
        cols.append(col / tot if normalize else col)
        labels.append("e[J=%s,F=%s,%s]" % (
            _dom(Excited, ie, 'J'), _dom(Excited, ie, 'F'),
            _PARITY[Excited.parities[ie]]))
    P = np.vstack(cols).T                                       # (nG, nSel)
    avg = P.mean(axis=1)                                        # equal weight

    g_all = np.arange(Ground.size)
    avg_df = pd.DataFrame({
        'g_N': [_dom(Ground, i, 'N') for i in g_all],
        'g_J': [_dom(Ground, i, 'J') for i in g_all],
        'g_F': [_dom(Ground, i, 'F') for i in g_all],
        'g_par': [_PARITY[Ground.parities[i]] for i in g_all],
        'BR': avg, 'g_idx': g_all})
    per_df = pd.DataFrame({'ground': [
        "X N=%s J=%s F=%s%s" % (_dom(Ground, i, 'N'), _dom(Ground, i, 'J'),
                                _dom(Ground, i, 'F'),
                                _PARITY[Ground.parities[i]]) for i in g_all]})
    for lab, col in zip(labels, cols):
        per_df[lab] = col
    return avg_df, per_df
```

- [ ] **Step 4: Run T3, verify it PASSES**

Run: `conda run -n Structure python "Jupyter Notebooks/RaX/test_xa_spectra.py"`
Expected: `T3 OK  sum(BR)=1.000000  N0=~1e-34`.

- [ ] **Step 5: Commit**

```bash
/Library/Developer/CommandLineTools/usr/bin/git add "Jupyter Notebooks/RaX/xa_spectra.py" "Jupyter Notebooks/RaX/test_xa_spectra.py" && /Library/Developer/CommandLineTools/usr/bin/git commit -m "feat(xa_spectra): averaged_branching (thesis Eq.3.4, equal-weight over excited F,M_F)" && /Library/Developer/CommandLineTools/usr/bin/git log -1 --stat
```

---

### Task 4: `broaden` + `doppler_fwhm` + `plot_spectrum`

**Goal:** Add lineshape convolution (Gaussian/Lorentzian/Voigt, area- or peak-normalized), the Doppler-width helper (thesis Eq. 3.15), and a one-call plot with optional experimental overlay.

**Files:**
- Modify: `Jupyter Notebooks/RaX/xa_spectra.py` (append three functions)
- Modify: `Jupyter Notebooks/RaX/test_xa_spectra.py` (append T4 block)

**Acceptance Criteria:**
- [ ] `doppler_fwhm` matches `2√(2 ln 2)·(ν₀/c)·√(k_B T/M)` to 1e-6 relative
- [ ] Area-normalized Gaussian kernel integrates to the summed line strength (±1%); peak-normalized single line peaks at its strength
- [ ] `broaden` clustering merges co-frequent lines (mean freq, summed strength)
- [ ] `plot_spectrum` returns a matplotlib Axes with stick + broadened + experimental artists, no exception (smoke)

**Verify:** `conda run -n Structure python "Jupyter Notebooks/RaX/test_xa_spectra.py"; echo EXIT=$?` → `T4 OK ...`, `EXIT=0`

**Steps:**

- [ ] **Step 1: Confirm the SciPy API is present** (stable since SciPy 1.4; verify in the env):

Run: `conda run -n Structure python -c "from scipy.special import voigt_profile; print('voigt ok')"`
Expected: `voigt ok`. (If absent, the plan's Voigt branch must fall back — but it is present in the `Structure` env's SciPy.)

- [ ] **Step 2: Append the failing test to `test_xa_spectra.py`** (before the final `print`):

```python
# ---- T4: broaden / doppler_fwhm / plot_spectrum ----
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fw = xs.doppler_fwhm(T=4.0, mass_amu=157.0, nu0=348.66e12)
ref = 2*np.sqrt(2*np.log(2)) * (348.66e12/299792458.0) * np.sqrt(
    1.380649e-23*4.0/(157.0*1.66053906660e-27))
assert abs(fw - ref)/ref < 1e-6, (fw, ref)

# area-normalized Gaussian: integral ~ total strength
lines = xs.line_list(g, e, g.select_q({'N': 1}), e.select_q({'J': 0.5}, parity='+'),
                     origin=e.parameters['Origin'])
x, y, (cp, cs) = xs.broaden(lines, shape='gaussian', fwhm=5.0,
                            norm='area', cluster=False)
area = np.trapz(y, x)
assert abs(area - lines['strength'].sum()) / lines['strength'].sum() < 0.01, area
# peak-normalized single line peaks at its strength
one = lines.iloc[[0]]
xx, yy, _ = xs.broaden(one, shape='gaussian', fwhm=5.0, norm='peak', cluster=False)
assert abs(yy.max() - one['strength'].iloc[0]) / one['strength'].iloc[0] < 1e-3
# clustering merges co-frequent lines
import pandas as pd
dup = pd.DataFrame({'freq': [100.0, 100.0, 300.0], 'strength': [1.0, 2.0, 4.0]})
_, _, (mp, ms) = xs.broaden(dup, shape='lorentzian', fwhm=10.0,
                            cluster=True, cluster_thresh=1.0)
assert len(mp) == 2 and abs(ms[0] - 3.0) < 1e-9 and abs(ms[1] - 4.0) < 1e-9, (mp, ms)
# plot_spectrum smoke (sticks + broadened + experimental)
fig, ax = plt.subplots()
expt = dict(freq=np.array([0.0, 50.0]), signal=np.array([1.0, 2.0]),
            err=np.array([0.1, 0.1]), offset=0, tweak=0, yscale=1.0)
ax2 = xs.plot_spectrum(lines=lines, ax=ax, sticks=True,
                       broaden_kw=dict(shape='voigt', fwhm=5.0, lorentz_fwhm=2.0),
                       experimental=expt, label='Sim')
assert ax2 is ax and len(ax.lines) >= 1
plt.close(fig)
print("T4 OK  doppler=%.3f MHz" % (fw/1e6))
```

- [ ] **Step 3: Append the three functions to `xa_spectra.py`:**

```python
_FWHM_SIG = 2.0 * np.sqrt(2.0 * np.log(2.0))     # FWHM = sigma * this (Gaussian)
_KB = 1.380649e-23
_AMU = 1.66053906660e-27
_C_SI = 299792458.0


def doppler_fwhm(T, mass_amu, nu0):
    """Doppler FWHM (thesis Eq.3.15). Returned in nu0's frequency unit."""
    sigma = (nu0 / _C_SI) * np.sqrt(_KB * T / (mass_amu * _AMU))
    return _FWHM_SIG * sigma


def _cluster(pos, s, thr):
    pos = np.asarray(pos, float); s = np.asarray(s, float)
    o = np.argsort(pos); pos, s = pos[o], s[o]
    if len(pos) == 0:
        return pos, s
    cut = np.where(np.diff(pos) >= thr)[0] + 1
    gp, gs = np.split(pos, cut), np.split(s, cut)
    return (np.array([p.mean() for p in gp]),
            np.array([q.sum() for q in gs]))


def broaden(lines, *, shape='gaussian', fwhm=None, lorentz_fwhm=None,
            T=None, mass_amu=None, nu0=None, x=None, norm='area',
            cluster=True, cluster_thresh=None, n_per_fwhm=20, pad=5):
    """Convolve a line list to a curve.

    shape : 'gaussian' | 'lorentzian' | 'voigt'
    fwhm  : explicit FWHM, OR omit and pass (T, mass_amu, nu0) for Doppler.
            For 'voigt', `fwhm` is the Gaussian (Doppler) FWHM and
            `lorentz_fwhm` the Lorentzian FWHM.
    norm  : 'area' (int g dnu = 1) | 'peak' (g_hat(0)=1)   (thesis p.111)
    cluster : merge lines closer than cluster_thresh (default fwhm/2).
    Returns (x, y_total, (component_pos, component_strength)).
    """
    if isinstance(lines, pd.DataFrame):
        pos = lines['freq'].to_numpy(float)
        s = lines['strength'].to_numpy(float)
    else:
        pos, s = np.asarray(lines[0], float), np.asarray(lines[1], float)
    if fwhm is None:
        if None in (T, mass_amu, nu0):
            raise ValueError("pass fwhm, or all of (T, mass_amu, nu0)")
        fwhm = doppler_fwhm(T, mass_amu, nu0)
    if len(pos) == 0:
        xx = np.array([]) if x is None else np.asarray(x, float)
        return xx, np.zeros_like(xx), (pos, s)
    if cluster:
        thr = fwhm / 2.0 if cluster_thresh is None else cluster_thresh
        pos, s = _cluster(pos, s, thr)
    if x is None:
        lo, hi = pos.min() - pad * fwhm, pos.max() + pad * fwhm
        n = max(int(np.ceil((hi - lo) / fwhm * n_per_fwhm)) + 1, 2)
        x = np.linspace(lo, hi, n)
    x = np.asarray(x, float)

    sig = fwhm / _FWHM_SIG
    if shape == 'gaussian':
        def k(d):
            g = np.exp(-0.5 * (d / sig) ** 2)
            return g / (sig * np.sqrt(2 * np.pi)) if norm == 'area' else g
    elif shape == 'lorentzian':
        hw = fwhm / 2.0
        def k(d):
            L = hw ** 2 / (d ** 2 + hw ** 2)             # peak = 1
            return (L / (np.pi * hw)) if norm == 'area' else L
    elif shape == 'voigt':
        from scipy.special import voigt_profile
        lf = fwhm if lorentz_fwhm is None else lorentz_fwhm
        gam = lf / 2.0
        v0 = voigt_profile(0.0, sig, gam)
        def k(d):
            v = voigt_profile(d, sig, gam)               # area-normalized
            return v if norm == 'area' else v / v0
    else:
        raise ValueError("shape must be gaussian|lorentzian|voigt")

    y = np.zeros_like(x)
    for p, a in zip(pos, s):
        y = y + a * k(x - p)
    return x, y, (pos, s)


def plot_spectrum(lines=None, *, Ground=None, Excited=None,
                  g_idx=None, e_idx=None, ax=None, sticks=True,
                  broaden_kw=None, normalize=True, experimental=None,
                  label='Simulation', **line_list_kw):
    """Plot sticks and/or a broadened curve, with an optional experiment.

    experimental : dict(freq, signal, err=None, offset=0, tweak=0, yscale=1)
        generalizes the BaF_spectrum_plot.py offset/tweak/yscale overlay.
    """
    import matplotlib.pyplot as plt
    if lines is None:
        lines = line_list(Ground, Excited, g_idx, e_idx, **line_list_kw)
    if ax is None:
        ax = plt.gca()
    if len(lines):
        smax = lines['strength'].max()
        norm_s = (lambda v: v / smax) if (normalize and smax > 0) else (lambda v: v)
        if broaden_kw:
            x, y, _ = broaden(lines, **broaden_kw)
            if normalize and y.size and y.max() > 0:
                y = y / y.max()
            ax.plot(x, y, label=label)
        if sticks:
            ax.vlines(lines['freq'].to_numpy(),
                      0.0, norm_s(lines['strength'].to_numpy()),
                      color='0.4', alpha=0.7, linewidth=1)
    if experimental is not None:
        ex = experimental
        xd = np.asarray(ex['freq'], float) + ex.get('tweak', 0)
        sg = np.asarray(ex['signal'], float)
        yd = sg / sg.max() * ex.get('yscale', 1.0)
        if ex.get('err') is not None:
            er = np.asarray(ex['err'], float) / sg.max() * ex.get('yscale', 1.0)
            ax.errorbar(xd, yd, yerr=er, linestyle='none', marker='o',
                        label='Data')
        else:
            ax.plot(xd, yd, linestyle='none', marker='o', label='Data')
    ax.legend(loc='best')
    return ax
```

- [ ] **Step 4: Run T4, verify it PASSES**

Run: `conda run -n Structure python "Jupyter Notebooks/RaX/test_xa_spectra.py"`
Expected: `T2 OK ...`, `T3 OK ...`, `T4 OK  doppler=... MHz`, exit 0.

- [ ] **Step 5: Commit**

```bash
/Library/Developer/CommandLineTools/usr/bin/git add "Jupyter Notebooks/RaX/xa_spectra.py" "Jupyter Notebooks/RaX/test_xa_spectra.py" && /Library/Developer/CommandLineTools/usr/bin/git commit -m "feat(xa_spectra): broaden (G/L/Voigt, area|peak), doppler_fwhm (Eq.3.15), plot_spectrum" && /Library/Developer/CommandLineTools/usr/bin/git log -1 --stat
```

---

### Task 5: Wire `BaF_spectrum_plot.py`, the cooling generator, and the snippet

**Goal:** `BaF_spectrum_plot.py` builds its simulation from structure via `xa_spectra` (keeping the measured `fluor_data` overlay); the cooling-notebook generator uses `averaged_branching`; deliver the user snippet. The two hand-built notebooks are not touched.

**Files:**
- Modify: `Jupyter Notebooks/RaX/BaF_spectrum_plot.py` (full rewrite)
- Modify: `Jupyter Notebooks/RaX/make_baf_xa_cooling_nb.py` (add averaged-BR cell)
- Regenerate: `Jupyter Notebooks/RaX/baf_xa_cooling.ipynb`
- Create: `Jupyter Notebooks/RaX/xa_spectra_snippet.md`

**Acceptance Criteria:**
- [ ] `BaF_spectrum_plot.py` runs headless and writes a PDF; no pasted `sim_freq`/`sim_data`; `fluor_data` retained and overlaid
- [ ] `make_baf_xa_cooling_nb.py` regenerates `baf_xa_cooling.ipynb` with an excited-averaged-BR table; `baf_xa_validate.py` still exits 0
- [ ] `xa_spectra_snippet.md` exists with a copy-paste cell for `RaF X-A.ipynb` / `BaF X-A.ipynb`
- [ ] `git status --porcelain "Jupyter Notebooks/RaX/RaF X-A.ipynb" "Jupyter Notebooks/RaX/BaF X-A.ipynb"` shows only the pre-existing ` M` (untouched by this task)

**Verify:** `conda run -n Structure python "Jupyter Notebooks/RaX/BaF_spectrum_plot.py" && conda run -n Structure python "Jupyter Notebooks/RaX/make_baf_xa_cooling_nb.py" && conda run -n Structure python "Jupyter Notebooks/RaX/baf_xa_validate.py"; echo EXIT=$?` → ends `VALIDATION PASSED`, `EXIT=0`

**Steps:**

- [ ] **Step 1: Capture the existing `fluor_data` arrays.** Run and keep the printed literal for Step 2:

Run: `conda run -n Structure python -c "import ast,io,sys; src=open('Jupyter Notebooks/RaX/BaF_spectrum_plot.py').read(); print(src[:src.index('sim_freq')])"`
Expected: prints the `import` lines + the full `fluor_data = [np.array([...]), np.array([...]), np.array([...])]` block. This block is reused verbatim.

- [ ] **Step 2: Rewrite `Jupyter Notebooks/RaX/BaF_spectrum_plot.py`** as (paste the captured `fluor_data = [...]` block where marked — it is unchanged from Step 1):

```python
"""BaF X-A simulated spectrum vs the 08/27/2025 fluorescence scan.

Simulation is generated from molecular structure via xa_spectra
(no pasted arrays). Freq offset is 348.66 THz.
"""
from config_path import add_to_sys_path
add_to_sys_path()                       # walk up to "Source Code"

import numpy as np
import matplotlib.pyplot as plt
from Energy_Levels import MoleculeLevels
import xa_spectra as xs

# --- measured data (unchanged: [freq(MHz)-offset, counts, stdev]) ---
# >>> PASTE the `fluor_data = [np.array([...]), ...]` block from Step 1 here <<<

offset = 348660000                      # 348.66 THz, the x-axis reference
tweak = -19                             # laser-frequency calibration shift
yscale = 1.3

def build(elec, N_list):
    return MoleculeLevels.initialize_state(
        molecule_name='BaF', elec_state=elec, vib_state=0,
        N_list=np.array(N_list), fermion_or_boson='boson',
        M_sublevels='none', I_nuclei=[0, 1/2], isotope=138,
        round=8, params=None, P_values=[1/2])

g = build('X', [1])                      # the N=1 feature in the 08/27 scan
e = build('A', [1])                      # A 2Pi1/2 J=1/2
g.eigensystem(0, 0); e.eigensystem(0, 0)

gidx = g.select_q({'N': 1})
eidx = e.select_q({'J': 0.5}, parity='+')
lines = xs.line_list(g, e, gidx, eidx, origin=e.parameters['Origin'])
lines = lines.assign(freq=lines['freq'] - offset)        # to the plot axis

fig, ax = plt.subplots()
xs.plot_spectrum(
    lines=lines, ax=ax, sticks=True,
    broaden_kw=dict(shape='gaussian', fwhm=40.0, cluster=True),
    normalize=True,
    experimental=dict(freq=fluor_data[0], signal=fluor_data[1],
                      err=fluor_data[2], tweak=tweak, yscale=yscale),
    label='Simulation')
ax.set_xlabel('Laser Frequency (MHz) - 348.66 THz')
ax.set_ylabel('Fluorescence (arb)')
fig.savefig('BaF_X_A_N1_sim_vs_data.pdf')
print("wrote BaF_X_A_N1_sim_vs_data.pdf;", len(lines), "lines")
```

- [ ] **Step 3: Run it, verify the rewrite works**

Run: `conda run -n Structure python "Jupyter Notebooks/RaX/BaF_spectrum_plot.py"`
Expected: `wrote BaF_X_A_N1_sim_vs_data.pdf; N lines` (N ≥ 1), file created, no exception.

- [ ] **Step 4: Add an averaged-BR cell to `make_baf_xa_cooling_nb.py`.** After the `CELL_BR` string definition, add a new markdown + code cell string:

```python
MD_AVG = """\
## 3b.&nbsp; Excited-state-averaged branching &nbsp;<span style="color:gray">(A hyperfine is small — often unresolved)</span>

Because the A²Π₁/₂ hyperfine splitting is only ≈ 22 MHz, the cooling laser
often does not resolve it. The table below averages the branching ratio over
**every** A²Π₁/₂(J=½,+) hyperfine sublevel with equal weight (thesis Eq. 3.4,
hyperfine-generalized). Σ over X N=1 = 1, Σ over X N=0 ≈ 0 (parity closure).
"""

CELL_AVG = """\
import xa_spectra as xs
avg_df, per_df = xs.averaged_branching(g, e, e.select_q({'J': 0.5}, parity='+'))
_n1 = avg_df['g_N'].astype(float) == 1
view = avg_df.loc[_n1, ['g_J', 'g_F', 'BR']].copy()
view['g_J'] = view['g_J'].map({0.5: '1/2', 1.5: '3/2'}).fillna(view['g_J'])
print('sum over N=1 =', round(float(avg_df.loc[_n1, 'BR'].sum()), 6),
      '  sum over N=0 =', '%.1e' % float(avg_df.loc[~_n1, 'BR'].sum()))
(view.set_index(['g_J', 'g_F']).style
 .format({'BR': '{:.4f}'})
 .background_gradient(subset=['BR'], cmap='Blues')
 .set_caption('Excited-averaged cooling BR  A²Π₁/₂(J=½,+) → X N=1'))
"""
```

Then insert these two cells into the `cells = [...]` list immediately after the `nbf.v4.new_code_cell(CELL_BR)` entry:

```python
    nbf.v4.new_markdown_cell(MD_BR),
    nbf.v4.new_code_cell(CELL_BR),
    nbf.v4.new_markdown_cell(MD_AVG),
    nbf.v4.new_code_cell(CELL_AVG),
    nbf.v4.new_code_cell(CELL_PLOT),
```

(The pre-existing list had `CELL_BR` then `CELL_PLOT`; the two new cells go between them. Leave every other cell unchanged.)

- [ ] **Step 5: Regenerate the notebook and re-verify the gate**

Run: `conda run -n Structure python "Jupyter Notebooks/RaX/make_baf_xa_cooling_nb.py"`
Expected: `wrote .../baf_xa_cooling.ipynb (12 cells)` (was 10; +2).

Run: `conda run -n Structure python "Jupyter Notebooks/RaX/baf_xa_validate.py"; echo EXIT=$?`
Expected: `VALIDATION PASSED`, `EXIT=0` (the validator is independent of the notebook; this confirms no Source-Code regression).

- [ ] **Step 6: Create `Jupyter Notebooks/RaX/xa_spectra_snippet.md`:**

```markdown
# xa_spectra — copy-paste cell for `RaF X-A.ipynb` / `BaF X-A.ipynb`

These two notebooks hold your own uncommitted work and were intentionally
not edited. Paste the cell below where the old `simulate_spectra` /
`simulate_spectra_noM` calls were. It replaces the drifting copy-pasted
helpers with the shared module.

```python
import xa_spectra as xs   # config_path.add_to_sys_path() must have run already

# spectrum (auto: no-M -> TDM^2; M-resolved -> line strength S; no Boltzmann)
gidx = g.select_q({'N': 1})
eidx = e.select_q({'J': 0.5}, parity='+')
lines = xs.line_list(g, e, gidx, eidx, origin=e.parameters['Origin'])
xs.plot_spectrum(lines=lines, sticks=True,
                 broaden_kw=dict(shape='voigt', fwhm=40.0, lorentz_fwhm=10.0))

# excited-state-averaged branching (A hyperfine often unresolved)
avg_df, per_df = xs.averaged_branching(g, e, eidx)
display(avg_df)
```

Old → new map:
`simulate_spectra(...)`           → `line_list` + `plot_spectrum`
`simulate_spectra_noM(...)`       → `line_list` (no-M auto-detected)
`plot_gaussian_spectrum(...)`     → `broaden(shape='gaussian')` / `plot_spectrum`
manual offset/tweak/yscale block  → `plot_spectrum(experimental=dict(...))`
```

- [ ] **Step 7: Confirm the hand-built notebooks were not touched**

Run: `/Library/Developer/CommandLineTools/usr/bin/git status --porcelain "Jupyter Notebooks/RaX/RaF X-A.ipynb" "Jupyter Notebooks/RaX/BaF X-A.ipynb"`
Expected: exactly the two pre-existing ` M` lines (the user's WIP), nothing staged.

- [ ] **Step 8: Commit (named files only — never `-A`)**

```bash
/Library/Developer/CommandLineTools/usr/bin/git add "Jupyter Notebooks/RaX/BaF_spectrum_plot.py" "Jupyter Notebooks/RaX/make_baf_xa_cooling_nb.py" "Jupyter Notebooks/RaX/baf_xa_cooling.ipynb" "Jupyter Notebooks/RaX/xa_spectra_snippet.md" && /Library/Developer/CommandLineTools/usr/bin/git commit -m "feat(notebooks): wire BaF_spectrum_plot + cooling-nb to xa_spectra; add averaged-BR cell + user snippet" && /Library/Developer/CommandLineTools/usr/bin/git log -1 --stat
```

---

### Task 6: Final verification gate + spec close-out

**Goal:** All gates green together; spec marked implemented.

**Files:**
- Modify: `docs/superpowers/specs/2026-05-16-xa-spectra-toolkit-design.md` (status line)

**Acceptance Criteria:**
- [ ] `test_xa_spectra.py` exits 0 (T2+T3+T4)
- [ ] `baf_xa_validate.py` exits 0
- [ ] Tutorial notebook exits 0, no `--allow-errors`
- [ ] `BaF_spectrum_plot.py` runs and writes its PDF
- [ ] Spec status line reads `implemented & validated 2026-05-17`

**Verify:** the four commands below all succeed; spec status updated.

**Steps:**

- [ ] **Step 1: Run the full gate set**

```bash
conda run -n Structure python "Jupyter Notebooks/RaX/test_xa_spectra.py"; echo T=$?
conda run -n Structure python "Jupyter Notebooks/RaX/baf_xa_validate.py"; echo V=$?
conda run -n Structure python "Jupyter Notebooks/RaX/BaF_spectrum_plot.py"; echo P=$?
conda run -n Structure jupyter execute "Jupyter Notebooks/RaF_Calcs_Tutorial.ipynb"; echo NB=$?
```
Expected: `T=0`, `V=0`, `P=0`, `NB=0`.

- [ ] **Step 2: Update the spec status line.** In `docs/superpowers/specs/2026-05-16-xa-spectra-toolkit-design.md` replace:

```
Status: designed 2026-05-16, awaiting user spec review.
```

with:

```
Status: implemented & validated 2026-05-17 (plan docs/superpowers/plans/2026-05-17-xa-spectra-toolkit.md; gates: test_xa_spectra/baf_xa_validate/tutorial all exit 0).
```

- [ ] **Step 3: Commit the close-out**

```bash
/Library/Developer/CommandLineTools/usr/bin/git add "docs/superpowers/specs/2026-05-16-xa-spectra-toolkit-design.md" && /Library/Developer/CommandLineTools/usr/bin/git commit -m "docs(spec): xa-spectra-toolkit implemented & validated 2026-05-17" && /Library/Developer/CommandLineTools/usr/bin/git log -1 --stat
```

---

## Self-Review

**Spec coverage:** §1 deliverables 1–4 → Tasks 1 / 2–4 / 5 / 5. §2 item-3 → Task 1 (+ verified deltas asserted). §3 line strength S → Task 2 `line_list`; ω³ `weight='rate'` → Task 2; averaged-BR Eq. 3.4 → Task 3; `doppler_fwhm`/Voigt/dual-norm → Task 4. §5 wiring + snippet → Task 5. §6 verification → every task's Verify + Task 6. §7 non-goals respected (no Boltzmann; hand-built notebooks untouched — Task 5 Step 7 asserts it). §8 file list matches the tasks.

**Placeholder scan:** no TBD/TODO; the one "paste the captured block" in Task 5 Step 2 is a verbatim-copy instruction with the exact source command in Step 1 (the array is ~2 KB of measured data, reproduced by command, not invented).

**Type consistency:** `LINE_COLUMNS`, `line_list`, `_has_M`, `_dom`, `_strength_matrix`, `averaged_branching` (returns `avg_df` with `g_N/g_J/g_F/g_par/BR/g_idx`, `per_df` with `ground` + per-sublevel cols), `broaden` (`shape/fwhm/lorentz_fwhm/norm/cluster`), `doppler_fwhm(T, mass_amu, nu0)`, `plot_spectrum` — names identical across Tasks 2–5 and the tests.
