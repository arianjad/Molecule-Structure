# ²²⁵RaF A(0)–X(0) hyperfine level diagram — design spec

**Status:** draft for review, 2026-06-05. **Owner:** A. Jadbabaie.
**Deliverable:** a paper/talk-quality energy-level diagram of the ²²⁵RaF laser-cooling
cycling transition with hyperfine/rotational branching, grouped by F₁.

---

## 1. Purpose

Show, as an energy-level diagram, the ²²⁵RaF cooling cycle

> **X²Σ⁺ (v=0, N=1, G=1, F₁=1)  →  A²Π₁/₂ (v=0, J′=½⁺, F₁=0)**

and the spontaneous-decay branching out of the excited **F₁=0** state, including the
leakage channels. This is the **hyperfine** companion to the existing ²²⁶RaF
**vibrational** cooling diagram (`~/Figures/raf-cooling-diagram`, `plot_v5.py`). The
²²⁵Ra nucleus (I=½, octupole-deformed) carries the extra **G = S + I_Ra** quantum
number that ²²⁶RaF (I=0) does not, so this is a structurally different figure, not a
re-render of the 226 one.

**Decisions locked with user (2026-06-05):**
- Resolution: **group by F₁** (sum branching over F within each F₁ manifold).
- Scope: **cycling manifold (N=1, G=1) + leakage** (N=1 G=0, and N=3).
- Use/style: **paper/talk figure** (tighter than the DAMOP-poster v5 scaling).
- Annotations: **levels + branching only** (no "12+2", "Γ_eff=2/7 Γ", or "1.75×" labels).
- Location: **Approach A** — build in `Molecule-Structure/Jupyter Notebooks/RaX/`.

## 2. Source of truth & exact setup

All numbers are **computed from Molecule-Structure**, reproducing
`RaX/RaF225 X-A.ipynb`. They are **not** read off the attached slide images. The
branching code's correctness is settled in
`RaX/reports/coherent_vs_incoherent_branching.md` (coherent ≡ incoherent for branching;
`M='none'` is the exact closed-form thesis Eq. 3.3, machine-precision validated).

Initialization (verbatim from the notebook):

```python
# Excited A²Π₁/₂ (v=0)
A0 = MoleculeLevels.initialize_state(
    molecule_name='RaF', elec_state='A', vib_state=0,
    N_list=np.arange(1,4), fermion_or_boson='fermion',
    M_sublevels='none', I_nuclei=[1/2, 1/2],   # [metal=Ra, ligand=F]
    isotope=225, round=8, params=None, P_values=[1/2, 3/2])
A0.update_params({'dYb': -A0.parameters['dYb'], 'h1/2Yb': 0}, recompute=True)  # see §3

# Ground X²Σ⁺ (v=0)
X0 = MoleculeLevels.initialize_state(
    molecule_name='RaF', elec_state='X', vib_state=0,
    N_list=np.arange(0,4), fermion_or_boson='fermion',
    M_sublevels='none', I_nuclei=[1/2, 1/2],
    isotope=225, round=8, params=None, P_values=[1/2])

BR0 = branching_ratios(X0, A0, 0, 1e-8)   # BR0[g,e]; columns ≈ 1 over v''=0 manifolds
```

- `BR0[g,e] = (G_evecs @ TDM @ E_evecs.T)²` — rows = ground, cols = excited
  (`Energy_Levels.branching_ratios`).
- Ground labels: `basis='bBS'`, `label_q=['N','G','F1','F']` → **case (b_βS)**:
  `G = S + I_Ra`, `F₁ = N + G`, `F = F₁ + I_F`.
- Excited labels: `label_q=['J','F1','F']` → **case (a)**: `F₁ = J + I_Ra`, `F = F₁ + I_F`.
- Cycling upper state = the single A `J=½, P=½, parity=+, F₁=0` eigenstate.

## 3. ⚠ Open physics point — param-set discrepancy (resolved empirically)

In the notebook, the **A-state energies** (cells 5/7/13 → the attached A image:
F₁=0 at ≈4.25 GHz, F₁=1 at ≈1.8 GHz) use the **cell-3-modified** params
(`dYb→−dYb`, `h1/2Yb→0`). But the **branching** (cells 26–27 → the heatmap) re-creates
`A0` from **database defaults** and does **not** re-apply cell 3. So the slide pairs an
energy diagram and a branching matrix computed under *different* A-state params.

`dYb`, `h1/2Yb` are **unpublished Skripnikov ab-initio** A-state params (no paper to
cite; per `Molecule-Structure/CLAUDE.md`).

**Resolution (no guessing):** compute branching **both** ways and compare the leak
fractions.
- If insensitive (expected — leaks are set by rotational/G mixing, not the small A
  hyperfine params), use the **cell-3-modified** set **consistently** for energies *and*
  branching (so the figure is internally consistent and the A energies match the image),
  and report the with/without comparison.
- If the modification **materially** changes the branching → **STOP and ask** which set
  is canonical before rendering.

## 4. Physics content of the figure

Selection-rule consequence (E1, ΔF₁=0,±1, 0↛0): **A(F₁=0) decays only to ground
F₁=1 manifolds.** Therefore within X(N=1,G=1) the **F₁=0 and F₁=2 levels are dark**
(not populated by spontaneous decay) — the hyperfine cycle is a closed **F₁=1 ↔ F₁=0**
loop. This is *verified* against the computed matrix (those columns must be ≈0), not
assumed.

Branching from the excited F₁=0 state, **grouped by ground (N,G,F₁)** (sum over F),
*conditional on decay to v″=0*:

| Decay target | grouped f (expected) | role |
|---|---|---|
| X(N=1, G=1, F₁=1) | ≈ 0.99995 | main return (closed loop) |
| X(N=1, G=1, F₁=0), F₁=2 | ≈ 0 | **dark** (selection rule) |
| X(N=1, **G=0**, F₁=1) | ≈ 5×10⁻⁵ | leak (repump target) |
| X(**N=3**, F₁=1) | < 10⁻⁸ | leak |

Exact values come from the compute step; the table is the expectation to verify against.
Total branching factorizes as **FCF_vib × (this hyperfine fraction)**; the figure states
the "conditional on v″=0" caveat so it is not misread as an absolute BR.

## 5. Layout (energy *not* to scale; intra-manifold ordering honored)

- **Top — A²Π₁/₂ (v=0, J′=½⁺):** the **F₁=0** cycling upper state, highlighted; F₁=1 of
  the same J doublet shown faintly below it for context (F₁=1 ≈1.8 GHz < F₁=0 ≈4.25 GHz).
  J=3/2 omitted (not the cycling state).
- **Bottom-center — X²Σ⁺ (v=0, N=1, G=1):** three F₁ bars, ordered F₁=1 < F₁=0 < F₁=2
  (per the X image, ≈25/70 · 350 · 355/400 MHz). F₁=1 is the cycling level; F₁=0 and F₁=2
  drawn but tagged *dark*.
- **Bottom-side (faint) — leakage manifolds:** X(N=1, G=0, F₁=1) and X(N=3) as the decay
  targets for the 5×10⁻⁵ and <10⁻⁸ leaks.
- **Transitions:** solid up-arrow = cycling pump (X G=1 F₁=1 → A F₁=0), labeled
  λ≈752.8 nm; wavy down-arrows = spontaneous decay — solid for the dominant return,
  dashed for the two leaks, each labeled `f = …`.
- **Labels:** `|G, F₁⟩` / `|J′, F₁⟩` mathtext matching the slide images; "energy not to
  scale" + "branching conditional on v″=0" notes.

## 6. Files (Approach A)

In `Molecule-Structure/Jupyter Notebooks/RaX/`:
- `compute_225_xa_branching.py` — `Structure` env; reproduces §2 setup, extracts F₁-level
  energies + the F₁=0-excited branching column (grouped by F₁), runs §3 comparison and §4
  verification, **freezes** exact numbers + provenance into a data module.
- `data_225_level_diagram.py` — frozen numbers + `PROVENANCE` (the only thing the renderer
  imports; keeps fast style iteration off the slow `Structure` env).
- `plot_225_level_diagram.py` — any-Python renderer; matplotlib; reads the data module;
  emits `figures/225RaF_XA_level_diagram_v1.{pdf,svg,png}` (RaX naming convention).

## 7. Verification plan

1. **Energies** match the attached images (X ≈25/70/350/355/400 MHz; A ≈0.6/1.45/1.8/4.25 GHz).
2. **Branching**: excited F₁=0 column sums to ≈1 over v″=0 manifolds; F₁=0/F₁=2 (N=1,G=1)
   columns ≈0 (selection rule); G=0 leak ≈5×10⁻⁵; N=3 leak <10⁻⁸ — matching the heatmap
   annotations.
3. **§3 param comparison** reported (with/without cell-3 modification).
4. **Adversarial physics review** (spectroscopy-reviewer: coupling-case labels, selection
   rules, parity, F₁-grouping) + an **independent recompute** of the F₁=0 branching column.
5. **I render the PNG and inspect it myself** before handing back (per the "read your own
   visual output" rule); then iterate styling with the user.

## 8. Out of scope (this pass)

Vibrational branching (covered by the 226 v5 figure); the J=3/2 A levels; M-sublevel /
Zeeman structure; the "12+2 / Γ_eff / 1.75×" annotations.
