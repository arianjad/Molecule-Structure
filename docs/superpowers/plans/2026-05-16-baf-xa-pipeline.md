# BaF A²Π₁/₂ ← X²Σ⁺ Pipeline Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers-extended-cc:subagent-driven-development (recommended) or superpowers-extended-cc:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add ¹³⁸Ba¹⁹F to the molecule database, reproduce arXiv:2511.06986 Table II frequencies (the (−)/N=0 line, F→F line-by-line) + the parity-closed laser-cooling branching table (A²Π₁/₂ v=0,J=½,+ → X²Σ⁺ v=0,N=1), via the validated RaF-style X–A pipeline. *(Scope corrected 2026-05-16: BR target is the cooling (+)↔N=1 closed system, not the physically-wrong (−)→N=0,1; see spec §1/§6.)*

**Architecture:** ¹³⁸Ba¹⁹F (I=0 + I=1/2) is structurally identical to RaF boson; the YbOH backend routes both through the molecule-agnostic `'174X000'/'174A000'` dispatch, so the *only* code change is two dicts in `molecule_parameters.py`. Everything else is `MoleculeLevels.initialize_state` + `branching_ratios` calls. Conventions/values are fixed by the approved spec (`docs/superpowers/specs/2026-05-16-baf-xa-pipeline-design.md`), triple-verified against B&C Eq.(8.511)/p.845 and Steimle PRA 84 012508.

**Tech Stack:** Python, conda env `Structure`, NumPy, the in-repo `Source Code/` library (`Energy_Levels.MoleculeLevels`, `molecule_parameters`), Jupyter.

**Spec reference:** `docs/superpowers/specs/2026-05-16-baf-xa-pipeline-design.md` (§4 = exact dicts, §6 = verification, §3.3 = the one acceptance-gated residual).

**Environment:** All Python runs use `conda run -n Structure` (CLAUDE.md). Repo is on Google Drive — chain dependent git ops in one shell call; use `/Library/Developer/CommandLineTools/usr/bin/git`.

---

### Task 1: Add ¹³⁸Ba¹⁹F boson X0 + A0 to molecule_parameters.py

**Goal:** Database entries so `get_molecule_params('BaF','X'|'A',0,'boson')` returns the spec §4 constants.

**Files:**
- Modify: `Source Code/molecule_parameters.py` (replace the dead stub line `# molecules['BaF']['boson']['X0'] = {`, currently at ~line 70, between the `molecules['RaF']['boson']['A0']` block and the `bF_225, c_225 = abinitio_2_effective_hyperfine(...)` line)

**Context the engineer needs:**
- Module scope already defines `c = params_general['c']  # = 29979.2458` and helper `b_2_bF(b, c) -> b + c/3` (B&C Eq. 8.511 Fermi contact).
- `get_molecule_params(molecule_name, elec_state, vib_state, fermion_or_boson)` builds the key `f'{elec_state}{vib_state}'` → `'X0'` / `'A0'`, indexes `molecules['BaF']['boson'][...]`, and merges `params_general`. `'BaF'` is already in `species_names`, so `molecules['BaF']['boson']` exists (empty) — adding the two dicts is the complete and sufficient change. No `molecule_library_class.py` edit (boson routes through generic `'174X000'/'174A000'`).
- Unit convention (spec §4, user decision): MHz directly from arXiv:2511.06986 Table III's MHz column; cm⁻¹ in the comment. `muE` keeps `<D>*0.503412`.

**Acceptance Criteria:**
- [ ] `get_molecule_params('BaF','X',0,'boson')` returns a dict with `bF` ≈ 66.2503, `c` == 8.224, `Be` == 6473.9588
- [ ] `get_molecule_params('BaF','A',0,'boson')` returns a dict with `Be` == 6347.847 (no −2Λ²D), `p+2q` == **-7719.00** (= -7713.96+2·(-2.52); spec §3.2/§4 wrote -7718.99 — hand-add slip, 0.01 MHz, immaterial), `a` == 29.32, `d` == -3.58
- [ ] `A0['Origin']` == `11946.109609 + 6347.847/c` (N² T₀,₀ + **Be_A** G-row, hand-written; **not** `formalism:'N2'`; **not** Be_X). Absolute-frequency validation deferred to Task 2 (Origin not exercised by this task's verify command)
- [ ] No `KeyError`; existing RaF/YbOH/CaOH entries untouched

**Verify:** `conda run -n Structure python -c "import sys; sys.path.insert(0,'Source Code'); from molecule_parameters import get_molecule_params as g; x=g('BaF','X',0,'boson'); a=g('BaF','A',0,'boson'); print(round(x['bF'],4), x['c'], x['Be']); print(a['Be'], a['p+2q'], a['a'], a['d'])"`
→ expected: `66.2503 8.224 6473.9588` / `6347.847 -7719.0 29.32 -3.58`

**Steps:**

- [ ] **Step 1: Replace the BaF stub with the two dicts**

In `Source Code/molecule_parameters.py`, replace the single line:

```python
# molecules['BaF']['boson']['X0'] = {
```

with (exact, from spec §4 — MHz-direct, comments carry cm⁻¹ + provenance):

```python
molecules['BaF']['boson']['X0'] = {
    'Be': 6473.9588,        # Ryzlewicz1980 [22]; 138Ba19F; 0.21594802 cm^-1
    'Gamma_SR': 80.984,     # Ryzlewicz1980 [22]; 0.0027013 cm^-1
    'bF': b_2_bF(63.509, 8.224),  # =66.2503; Frosch-Foley b,c Ernst1986 [23]; bF=b+c/3 (B&C Eq.8.511 p.605/p.845; Steimle PRA84 012508 TblIII bF(F)=66.25 direct)
    'c': 8.224,             # Ernst1986 [23]; 0.0002743 cm^-1
    'D': 0.0055250,         # Effantin1990 [18]; 1.84294e-7 cm^-1
    'muE': 3.170*0.503412,  # 3.170(3) D, Steimle PRA84 012508 Table V; Debye->MHz/(V/cm)
    }
molecules['BaF']['boson']['A0'] = {
    'Be': 6347.847,         # Steimle2011 [24], 0.2117414 cm^-1, held fixed in fit. Entered direct (NO -2*Lam^2*D): the R^2<->N^2 B-row term is +-2*Lam^2*D = 0.012 MHz << 0.3 MHz exp.unc & all validation tols -> immaterial regardless of B's source convention. The ~6348 MHz Lam^2*B convention shift lives entirely in Origin (G-row), not here. (RaF A0 differs: its source B was N^2 at a level the converter's -2Lam^2D handles.)
    'ASO': 18955512.5,      # This work arXiv:2511.06986; 632.287838 cm^-1
    'h1/2': 0,
    'a': 26.55 - 0.5*(-5.54),  # =29.32; arXiv:2511.06986 Eq.(4)/B&C p.845 h1/2=a-(b+c)/2; Denis2022 [39]. Gated on Table II 21.8 MHz (plan Task 2)
    'bF': 0,
    'c': 0,
    'd': -3.58,             # Denis2022 [39]; pgopher->B+C needs minus (RaF A0 precedent)
    'p+2q': -7713.96 + 2*(-2.52),  # =-7719.00 (spec/§4 wrote -7718.99: hand-add slip, -7713.96-5.04=-7719.00; 0.01 MHz, immaterial); Effantin1990 [18] combined p+2q (op-grounded hamiltonian_builders.py:147; Steimle (p+2q)=-7721 cross-check)
    'q': -2.52,             # Effantin1990 [18]; -8.40e-5 cm^-1
    'p2q_D': -0.00699,      # Effantin1990 [18] p_D; -2.332e-7 cm^-1
    'D': 0.006007,          # Effantin1990 [18]; 2.0036e-7 cm^-1
    'A_D': 0.93,            # Steimle2011 [24]; consumed iff A2Pi builder reads it (plan Task 2 confirms)
    'muE': 1.50*0.503412,   # 1.50(2) D A2Pi1/2, Steimle PRA84 012508 Table V; Debye->MHz/(V/cm)
    'g_S': 2.0023,
    'Origin': 11946.109609 + 6347.847/c,  # T0,0[cm^-1] (arXiv:2511.06986 Tbl III, "This work"; pgopher-DEFAULT N^2 fit) + Be_A[MHz]/c R^2 G-row. Code rot op is R^2-form B*(N^2-Lam^2); paper T0,0 is N^2-convention (pgopher default; paper p.8: 0.21 cm^-1~B_A offset vs Steimle). Be_A=6347.847 (THIS state's own B) NOT Be_X=6473.9588 (Be_X was the +126 MHz bug class the converter/RaF-migration removed). NOT formalism:'N2' (would mis-apply -2Lam^2D to R^2-form Be/p+2q, spec §3.4). Abs line confirmed vs Table II 348666424.4 MHz in Task 2.
    }
```

- [ ] **Step 2: Run the verify command**

Run: `conda run -n Structure python -c "import sys; sys.path.insert(0,'Source Code'); from molecule_parameters import get_molecule_params as g; x=g('BaF','X',0,'boson'); a=g('BaF','A',0,'boson'); print(round(x['bF'],4), x['c'], x['Be']); print(a['Be'], a['p+2q'], a['a'], a['d'])"`
Expected: `66.2503 8.224 6473.9588` then `6347.847 -7719.0 29.32 -3.58`

- [ ] **Step 3: Commit**

```bash
GIT=/Library/Developer/CommandLineTools/usr/bin/git
"$GIT" add "Source Code/molecule_parameters.py" && "$GIT" commit -m "feat(molecule_parameters): add 138Ba19F boson X0/A0 (arXiv:2511.06986; B&C/Steimle-verified conventions)" && "$GIT" log -1 --oneline
```

---

### Task 2: BaF X–A no-M validation script + resolve the §3.3 hyperfine residual

**Goal:** A runnable script asserting the spec §6 checks: (i) X²Σ⁺(N=0,J=½) F-split ≈ 65.6 MHz and A²Π₁/₂(J=½,−) F-split ≈ 21.8 MHz (both ≤1 MHz; the §3.3 `'a'`/`'d'` gate); (ii) the 3 Table II (−)-line components line-by-line, matched F_X→F_A, ≤2 MHz; (iii) the parity-closed cooling-line BR A²Π₁/₂(J=½,+)→X N=1 (Σ→N=1=1.0, Σ→N=0≤1e-3, spreads over the N=1 hf manifold). **Status: implemented + green 2026-05-16** — §3.3 passed first try (21.928 MHz); freqs +0.6/+0.8/+1.4 MHz; parity closure exact (→N=0 ~1e-34).

**Files:**
- Create: `Jupyter Notebooks/RaX/baf_xa_validate.py`
- Possibly modify (only if Step 3 fails): `Source Code/molecule_parameters.py` (`molecules['BaF']['boson']['A0']` `'a'`/`'d'`) and read `Source Code/matrix_elements.py` (the `IzLz` / `T2_2(I,S)` aBJ ²Π₁/₂ element) — spec §3.3.

**Context the engineer needs:**
- `MoleculeLevels.initialize_state(molecule_name, elec_state, vib_state, N_list, fermion_or_boson='boson', M_sublevels='none', I_nuclei=[0,1/2], isotope=138, round=8, params=None, P_values=[1/2])` — mirror the RaF X-A notebook call (`Jupyter Notebooks/RaX/RaF X-A.ipynb` cell 5) with `molecule_name='BaF'`, `isotope=138`.
- `g.eigensystem(0,0)` populates `g.evals0` (MHz) / `g.evecs0`. `g.select_q({'N':0,'J':0.5})` returns index array; `parity='-'` filters by Λ-doublet parity (spec: the observed A²Π₁/₂ component is `(−)`).
- `branching_ratios(Ground, Excited, Ez, Bz)` (from `Energy_Levels`) returns `(G_evecs @ TDM @ E_evecs.T)**2`, shape `(len(Ground), len(Excited))`, unnormalized; for no-M it auto-uses the `_noM` TDM builder. Memory `raf-xa-branching-tdm-validated`: this machinery is correct — do not "fix" it; only the BaF A0 hyperfine mapping is under test.
- The two F-splittings are intra-state and **Origin-independent** (differences of zero-field eigenvalues), so they isolate `bF` (X) and the `'a'`/`'d'` mapping (A) cleanly.
- Discriminator targets (spec §6, from arXiv:2511.06986 Table II): X(N=0,J=1/2) F=0↔1 = **65.6 MHz**; A²Π₁/₂(J=1/2,−) F=0↔1 = **21.8 MHz**. Tolerance ≤1 MHz (paper's stated absolute accuracy; exp. unc. 0.3 MHz).

**Acceptance Criteria:**
- [ ] Script runs under `conda run -n Structure` with no exception
- [ ] `abs(x_split - 65.6) <= 1.0` (X Fermi-contact mapping correct)
- [ ] `abs(a_split - 21.8) <= 1.0` (A²Π₁/₂ hyperfine mapping correct — spec §3.3 gate)
- [ ] `abs(nu_centroid - 348666439.0) <= 5.0` (absolute X→A line confirms the N² T₀,₀ + Be_A G-row Origin — spec §5; a ~6348 MHz miss = inverted N²/R², ~126 MHz = Be_X-not-Be_A)
- [ ] No-M BR matrix is finite, non-negative; each A²Π₁/₂(J=1/2,−) hyperfine level's normalized decay spreads over X N=0 **and** N=1 (both > 0) and sums to 1.0 ± 1e-6
- [ ] If `'a'`/`'d'` were changed, the spec §3.3 STOP-comment in `molecule_parameters.py` is updated to record the resolved mapping + that 21.8 MHz now passes

**Verify:** `conda run -n Structure python "Jupyter Notebooks/RaX/baf_xa_validate.py"` → expected final line `VALIDATION PASSED` and printed `x_split≈65.6`, `a_split≈21.8`.

**Steps:**

- [x] **Steps 1–5: DONE 2026-05-16.** The committed `Jupyter Notebooks/RaX/baf_xa_validate.py` is the source of truth — refocused to the cooling-line scope (spec §1/§6): the two splitting gates (65.6/21.8 ≤1 MHz), the 3 Table II (−)-line components line-by-line (≤2 MHz; got +0.6/+0.8/+1.4), and the parity-closed cooling BR A²Π₁/₂(J=½,+)→X N=1 (Σ→N=1=1.0, Σ→N=0~1e-34, spreads over the N=1 hf manifold). §3.3 passed first try (`'a'`=29.32, `'d'`=−3.58 → 21.928 MHz; no mapping change needed). The original (−)→N=0,1 script + Steps 2–5 below are **superseded historical scaffolding** (kept for trace; do not re-execute).

<details><summary>superseded original Step 1 script + Steps 2–5 (historical)</summary>

Create `Jupyter Notebooks/RaX/baf_xa_validate.py`:

```python
"""BaF X-A no-M validation against arXiv:2511.06986 Table II (spec §6)."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'Source Code'))
import numpy as np
from Energy_Levels import MoleculeLevels, branching_ratios

def build(elec, N_list):
    return MoleculeLevels.initialize_state(
        molecule_name='BaF', elec_state=elec, vib_state=0,
        N_list=np.array(N_list), fermion_or_boson='boson',
        M_sublevels='none', I_nuclei=[0, 1/2], isotope=138,
        round=8, params=None, P_values=[1/2])

g = build('X', [0, 1])     # X2Sig+ v=0, N=0,1
e = build('A', [1, 2])     # A2Pi1/2 v=0 (contains J=1/2,3/2)
g.eigensystem(0, 0)
e.eigensystem(0, 0)

# X(N=0, J=1/2): F=0,1 splitting ~ Fermi contact bF ~ 65.6 MHz
gx = g.select_q({'N': 0, 'J': 0.5})
gE = np.sort(g.evals0[gx])
x_split = float(abs(gE[-1] - gE[0]))

# A2Pi1/2(J=1/2, parity -): F=0,1 splitting ~ 21.8 MHz
ex = e.select_q({'J': 0.5}, parity='-')
eE = np.sort(e.evals0[ex])
a_split = float(abs(eE[-1] - eE[0]))

print(f"X  N=0,J=1/2 F-split = {x_split:.3f} MHz  (target 65.6, tol 1.0)")
print(f"A  2Pi1/2 J=1/2(-) F-split = {a_split:.3f} MHz  (target 21.8, tol 1.0)")

# no-M branching: A2Pi1/2(J=1/2,-) -> X(N=0,1)
BR = branching_ratios(Ground=g, Excited=e, Ez=0, Bz=0)
assert np.all(np.isfinite(BR)) and np.all(BR >= 0), "BR not finite/non-negative"
gN = np.array(g.q_numbers['N'])
for col in ex:
    tot = BR[:, col].sum()
    assert tot > 0, f"A level {col} has zero total decay"
    p = BR[:, col] / tot
    n0 = p[gN[np.argmax(g.evecs0**2, axis=1)] == 0].sum()
    n1 = p[gN[np.argmax(g.evecs0**2, axis=1)] == 1].sum()
    assert abs(p.sum() - 1.0) <= 1e-6, "normalized BR != 1"
    assert n0 > 0 and n1 > 0, f"A level {col}: decay missing an N branch (N=0:{n0}, N=1:{n1})"
    print(f"  A idx {col}: BR to X N=0 = {n0:.4f}, N=1 = {n1:.4f}")

# Absolute-frequency convention confirmation (spec §5; user decision: N² T0,0 + Be_A G-row).
# Common-mode to all hf components, so use the manifold centroid; a wrong N²/R² Origin
# choice misses by Be_A·Λ² ~ 6348 MHz (>> 5 MHz tol); Be_X-vs-Be_A bug ~ 126 MHz (also caught).
c_cm = 29979.2458                                       # == molecule_parameters.c
Ec = e.evals0[ex]/c_cm + e.parameters['Origin']         # A 2Pi1/2(J=1/2,-) levels, cm^-1
Gc = g.evals0[gx]/c_cm                                  # X N=0,J=1/2 levels, cm^-1
nu_centroid = c_cm*(Ec.mean() - Gc.mean())              # MHz, hf-centroid X(N=0)->A(J=1/2,-)
TBL2_CENTROID = (348666402.6 + 348666424.4 + 348666490.0)/3   # = 348666439.0 MHz (Table II)
print(f"abs X(N=0)->A(2Pi1/2,J=1/2,-) centroid = {nu_centroid:.1f} MHz "
      f"(Table II {TBL2_CENTROID:.1f}; tol 5 MHz)")
assert abs(nu_centroid - TBL2_CENTROID) <= 5.0, (
    f"absolute line off by {nu_centroid - TBL2_CENTROID:.1f} MHz — Origin N²/R² "
    f"convention wrong (~6348 MHz miss => use the other Origin form; spec §5)")

assert abs(x_split - 65.6) <= 1.0, f"X bF mapping off: {x_split} (expect 65.6)"
assert abs(a_split - 21.8) <= 1.0, f"A 2Pi1/2 hyperfine mapping off: {a_split} (expect 21.8) — spec §3.3 residual"
print("VALIDATION PASSED")
```

- [ ] **Step 2: Run it — expect X to pass, A possibly failing (the §3.3 gate)**

Run: `conda run -n Structure python "Jupyter Notebooks/RaX/baf_xa_validate.py"`
Expected: `x_split` within 1 MHz of 65.6 (confirms Task 1 `bF`). If `a_split` assertion fails → the spec §3.3 residual is live; go to Step 3. If both pass → skip to Step 4.

- [ ] **Step 3: (only if A assertion failed) Resolve the §3.3 `'a'`/`'d'` mapping against the data**

Inspect the aBJ ²Π₁/₂ hyperfine matrix elements to see what scalar `params['a']` (×`IzLz`) and `params['d']` (×`T2_2(I,S)`) are expected to be (spec §3.3 — the codebase `a·IzLz` normalization question):

Run: `conda run -n Structure python -c "import sys; sys.path.insert(0,'Source Code'); import inspect, matrix_elements as me; print(inspect.getsource(me.IzLz_aBJ)); print(inspect.getsource(me.T2_2_aBJ if hasattr(me,'T2_2_aBJ') else me.T2ISq2_aBJ))"`
(Adjust the element names to the actual aBJ ²Π hyperfine elements wired in `molecule_library_class.collect_all_matrix_elements` for the A state — grep `aBJ_even_A_matrix_elements` in that file.)

Then, treating the **21.8 MHz Table II splitting as the arbiter** (spec: "test it, the data settles it"), adjust ONLY `molecules['BaF']['boson']['A0']['a']` and/or `['d']` in `Source Code/molecule_parameters.py`. Candidate forms to try, in order (each is a one-line change + re-run Step 2):
  1. `'a': 26.55 - 0.5*(-5.54)` (=29.32; B&C h1/2=a−(b+c)/2; current)
  2. `'a': 26.55` with the `(b+c)` carried via `'bF'`/`'c'` on the A-state if the builder consumes them there
  3. flip `'d'` sign (`+3.58`) if the A F-ordering is inverted

Do not tune to arbitrary values — only switch between these physically-defined parametrizations until `abs(a_split-21.8)<=1.0`. Record the winning form in the inline comment, replacing the "Gated on Table II" note with e.g. `# =29.32; reproduces Table II 21.8 MHz (verified plan Task 2)`.

- [ ] **Step 4: Re-run to confirm green**

Run: `conda run -n Structure python "Jupyter Notebooks/RaX/baf_xa_validate.py"`
Expected: ends with `VALIDATION PASSED`; both splittings within 1 MHz.

- [ ] **Step 5: Commit**

```bash
GIT=/Library/Developer/CommandLineTools/usr/bin/git
"$GIT" add "Jupyter Notebooks/RaX/baf_xa_validate.py" "Source Code/molecule_parameters.py" && "$GIT" commit -m "test(BaF): Table II 21.8/65.6 MHz validation; resolve A2Pi1/2 hyperfine mapping" && "$GIT" log -1 --oneline
```

</details>

---

### Task 3: BaF X–A analysis notebook (deliverable, cooling-line)

**Goal:** A notebook producing hyperfine-resolved energies, the **parity-closed cooling-line BR table** (A²Π₁/₂ v=0,J=½,+ → X²Σ⁺ v=0,N=1), and the 3 Table II (−)-line components line-by-line (F_X→F_A). Mirrors the committed `baf_xa_validate.py` structure/scope (not the superseded (−)→N=0,1 form).

**Files:**
- Create: `Jupyter Notebooks/RaX/BaF X-A.ipynb` (zone has its own CLAUDE.md — analysis artifact, not `Source Code/`)
- Create (generator, kept for reproducibility): `Jupyter Notebooks/RaX/make_baf_xa_nb.py`

**Context the engineer needs:**
- Build the notebook programmatically with `nbformat` (already in the `Structure` env) so it is deterministic and re-runnable headless via `jupyter execute`.
- Absolute transition frequency = `c*(E_A_cm - E_X_cm)` style; the notebook prints relative splittings (origin-independent, the real validation) and the absolute line near 348666424 MHz using `e.parameters['Origin']` (carries the spec §5 convention offset — label it as approximate).
- `c = 29979.2458` (MHz per cm⁻¹).

**Acceptance Criteria:**
- [ ] `BaF X-A.ipynb` exists and executes end-to-end on a fresh kernel with no error
- [ ] It prints: X(N=0,J=½) & A²Π₁/₂(J=½,−) F-splittings (65.6/21.8); the 3 Table II (−)-line components computed-vs-reference line-by-line (F1→F0 402.6, F1→F1 424.4, F0→F1 490.0; F0→F0 forbidden); and the parity-closed cooling-line BR table A²Π₁/₂(J=½,+)→X N=1 (showing Σ→N=0≈0)
- [ ] Generator script `make_baf_xa_nb.py` reproduces the notebook (its cells mirror the committed `baf_xa_validate.py` cooling-line logic; the embedded generator below is to be rewritten to that scope during execution)

**Verify:** `conda run -n Structure jupyter execute "Jupyter Notebooks/RaX/BaF X-A.ipynb"` → exit code 0

**Steps:**

- [ ] **Step 1: Write the notebook generator**

Create `Jupyter Notebooks/RaX/make_baf_xa_nb.py`:

```python
"""Generate 'BaF X-A.ipynb' (mirrors RaF X-A.ipynb structure)."""
import nbformat as nbf, os

cells_src = [
("md", "# BaF A²Π₁/₂ ← X²Σ⁺ — hyperfine energies + no-M branching\n"
       "Validation vs arXiv:2511.06986 Table II. Constants: spec 2026-05-16."),
("code",
 "from config_path import add_to_sys_path\n"
 "add_to_sys_path()\n"
 "import numpy as np\n"
 "from Energy_Levels import MoleculeLevels, branching_ratios\n"
 "np.set_printoptions(precision=5, suppress=True)\n"
 "c = 29979.2458"),
("code",
 "def build(elec, N_list):\n"
 "    return MoleculeLevels.initialize_state(molecule_name='BaF', elec_state=elec,\n"
 "        vib_state=0, N_list=np.array(N_list), fermion_or_boson='boson',\n"
 "        M_sublevels='none', I_nuclei=[0,1/2], isotope=138, round=8,\n"
 "        params=None, P_values=[1/2])\n"
 "g = build('X', [0,1]); e = build('A', [1,2])\n"
 "g.eigensystem(0,0); e.eigensystem(0,0)"),
("code",
 "gx = g.select_q({'N':0,'J':0.5}); gE = np.sort(g.evals0[gx])\n"
 "ex = e.select_q({'J':0.5}, parity='-'); eE = np.sort(e.evals0[ex])\n"
 "x_split = abs(gE[-1]-gE[0]); a_split = abs(eE[-1]-eE[0])\n"
 "print(f'X  N=0,J=1/2 F-split  = {x_split:.3f} MHz (Table II: 65.6)')\n"
 "print(f'A  2Pi1/2 J=1/2(-) F-split = {a_split:.3f} MHz (Table II: 21.8)')"),
("code",
 "BR = branching_ratios(Ground=g, Excited=e, Ez=0, Bz=0)\n"
 "gN = np.array(g.q_numbers['N'])[np.argmax(g.evecs0**2, axis=1)]\n"
 "for col in ex:\n"
 "    p = BR[:,col]/BR[:,col].sum()\n"
 "    print(f'A idx {col}: no-M BR  X(N=0)={p[gN==0].sum():.4f}  X(N=1)={p[gN==1].sum():.4f}')"),
("code",
 "ref = {'F1->F0':348666402.6,'F1->F1':348666424.4,'F0->F1':348666490.0}\n"
 "print('arXiv:2511.06986 Table II (MHz):', ref)\n"
 "print('Computed X F-split %.3f MHz vs Table II 65.6; A F-split %.3f MHz vs 21.8' % (x_split, a_split))"),
]
nb = nbf.v4.new_notebook()
nb.cells = [nbf.v4.new_markdown_cell(s) if t=="md" else nbf.v4.new_code_cell(s)
            for t, s in cells_src]
out = os.path.join(os.path.dirname(__file__), "BaF X-A.ipynb")
nbf.write(nb, out)
print("wrote", out)
```

- [ ] **Step 2: Generate the notebook**

Run: `cd "Jupyter Notebooks/RaX" && conda run -n Structure python make_baf_xa_nb.py`
Expected: `wrote .../BaF X-A.ipynb`

- [ ] **Step 3: Execute headless (fresh kernel)**

Run: `conda run -n Structure jupyter execute "Jupyter Notebooks/RaX/BaF X-A.ipynb"`
Expected: exit code 0, no exception.

- [ ] **Step 4: Commit**

```bash
GIT=/Library/Developer/CommandLineTools/usr/bin/git
"$GIT" add "Jupyter Notebooks/RaX/BaF X-A.ipynb" "Jupyter Notebooks/RaX/make_baf_xa_nb.py" && "$GIT" commit -m "feat(notebooks): BaF X-A analysis notebook (energies + no-M BR vs Table II)" && "$GIT" log -1 --oneline
```

---

### Task 4: Verification gate + spec close-out

**Goal:** Confirm no regression to the project verification gate, and close the resolved spec open items.

**Files:**
- Modify: `docs/superpowers/specs/2026-05-16-baf-xa-pipeline-design.md` (§7 Open items: mark #1 resolved with the Task 2 outcome; note #4 confirmed — no dispatch wiring needed; #5 A_D consumption finding)

**Context the engineer needs:**
- CLAUDE.md single completeness gate for any `Source Code/` change: the tutorial notebook runs end-to-end on a fresh kernel. The only `Source Code/` change here is additive DB entries, but the gate must still pass (no `--allow-errors`, per memory `feedback_no_allow_errors`).

**Acceptance Criteria:**
- [ ] `conda run -n Structure jupyter execute "Jupyter Notebooks/RaF_Calcs_Tutorial.ipynb"` exits 0 (no regression; no `--allow-errors`)
- [ ] Spec §7: item #1 (§3.3 `'a'` residual) updated with the Task 2 result (resolved value or still-open with evidence); item #4 marked "confirmed — generic `'174X000'/'174A000'` path, no wiring"; item #5 records whether the A²Π builder consumed `A_D`
- [ ] Spec §5 Origin framing corrected: paper T₀,₀ = pgopher-default N²; `Origin = T₀,₀ + Be_A·Λ²/c` (hand G-row, not `formalism:'N2'`, not Be_X); the old "documented offset / RaF +B/c template" wording removed
- [ ] Spec status header changed from `awaiting user review before writing-plans` to `implemented <date>`

**Verify:** `conda run -n Structure jupyter execute "Jupyter Notebooks/RaF_Calcs_Tutorial.ipynb"; echo $?` → final line `0`

**Steps:**

- [ ] **Step 1: Run the project verification gate**

Run: `conda run -n Structure jupyter execute "Jupyter Notebooks/RaF_Calcs_Tutorial.ipynb"; echo EXIT=$?`
Expected: `EXIT=0`. If non-zero, the DB edit regressed the tutorial — diff `molecule_parameters.py`, fix, re-run (do not pass `--allow-errors`).

- [ ] **Step 2: Update spec open items**

Edit `docs/superpowers/specs/2026-05-16-baf-xa-pipeline-design.md` §7:
- Item 1: replace "*Highest-risk item.*" with the Task 2 verdict, e.g. `RESOLVED — 'a'=29.32 reproduces Table II 21.8 MHz (Task 2, <date>)` (or the corrected form actually used).
- Item 4: `RESOLVED — BaF boson uses the generic '174X000'/'174A000' dispatch; no molecule_library_class.py change needed (confirmed Task 1).`
- Item 5: record whether `A_D` is consumed by the boson A²Π builder (from Task 2 Step 3 inspection); if unused, note it's inert/kept for documentation.
- §5 Origin: replace the "documented absolute offset / RaF-A0-template +B/c" framing with the resolved finding — paper Eq.(3) is R²-*symbolic* but the fit/Table III is **pgopher-default N²** (corroborated by the p.8 0.21 cm⁻¹≈B_A offset vs Steimle [24]); code operator is R²-form, so `Origin = T₀,₀ + Be_A·Λ²/c` (hand-written G-row with **this state's own B**; **not** `formalism:'N2'`, **not** Be_X). Confirmed vs Table II (Task 2, centroid ≤5 MHz). Also nuance §3.4: its "Eq.(3) R²-form" is the *symbolic* operator only (the fit is pgopher-default N²); the no-−2Λ²D **decision stands** because the R²↔N² B-row term is 0.012 MHz (≪ 0.3 MHz exp. unc.), immaterial — record this bound rather than the "cleanly R²-form" rationale.
- Header `**Status:**` → `implemented <YYYY-MM-DD>`.

- [ ] **Step 3: Commit**

```bash
GIT=/Library/Developer/CommandLineTools/usr/bin/git
"$GIT" add "docs/superpowers/specs/2026-05-16-baf-xa-pipeline-design.md" && "$GIT" commit -m "docs(spec): close BaF pipeline open items; record verification-gate pass" && "$GIT" log -1 --oneline
```

---

## Self-Review

**Spec coverage:**
- Spec §4 (DB dicts, MHz-direct, μ from Steimle) → Task 1. ✓
- Spec §6 (Table II 65.6/21.8 MHz discriminators, BR sum=1/structure, no quantitative intensity) → Task 2. ✓
- Spec §3.3 (the one acceptance-gated residual) → Task 2 Step 3, arbitrated by the 21.8 MHz data. ✓
- Spec §3.4 (Be R²-form, no −2Λ²D) → Task 1 dict + comment. ✓
- Spec §5 (Origin convention) → **resolved**: paper T₀,₀ is N²-convention (pgopher-default fit; verified vs paper Eq.3/p.8 + the 0.21 cm⁻¹≈B_A Steimle offset). Origin = `11946.109609 + 6347.847/c` (N² T₀,₀ + **Be_A** G-row by hand; **not** `formalism:'N2'` — would corrupt R²-form Be/p+2q per §3.4). Task 1 dict + Task 2 absolute Table II confirmation (centroid ≤5 MHz). Spec §5's "documented offset / RaF-template +B/c" framing is a misread of the paper-vs-Steimle note — corrected in Task 4. ✓
- Spec §3.6 (A–B perturbation → BR not quantitatively matched) → Task 2 BR acceptance is structure/normalization only. ✓
- CLAUDE.md verification gate → Task 4. ✓
- Deliverable (RaF-style notebook) → Task 3. ✓
- Spec "Out of scope" (v=1, A²Π₃/₂, non-138, field-on) → not in any task. ✓ (correctly excluded)

**Placeholder scan:** No `TBD`/`later`/"handle errors". The only spec `<…>` (Origin) is resolved to a concrete formula in Task 1. Code blocks are complete and runnable. ✓

**Type/name consistency:** `build(elec, N_list)`, `g`/`e`, `select_q`, `branching_ratios(Ground=,Excited=,Ez=,Bz=)`, `evals0`, `q_numbers['N']`, `parity='-'` used identically across Tasks 2 and 3. Discriminator targets (65.6, 21.8) and tolerance (1.0 MHz) identical to spec §6. ✓
