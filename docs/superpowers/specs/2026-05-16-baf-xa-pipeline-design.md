# BaF A²Π₁/₂ ← X²Σ⁺ pipeline — design spec

**Date:** 2026-05-16
**Branch:** `fix-for-distrib`
**Status:** design — awaiting user review before `writing-plans`

## 1. Goal & scope

Add `¹³⁸Ba¹⁹F` to `Source Code/molecule_parameters.py` and run the existing
RaF-style X–A workflow to produce, for the laser-cooling / eEDM transition:

- Hyperfine-resolved level energies for **X²Σ⁺(v=0), N=0,1** and
  **A²Π₁/₂(v=0), J=1/2** (both Λ-doublet parities; the spectroscopically
  observed one is the `(−)` component, per paper §IV).
- Zero-field (`M_values='none'`, "no-M") **branching-ratio table** for the
  laser-cooling transition **A²Π₁/₂(v=0, J=1/2, +) → X²Σ⁺(v=0, N=1)** — the
  rotationally-closed system (A `+`-parity decays E1 only to X `−`-parity =
  odd N=1; N=0/N=2 are `+`-parity, forbidden). *(Scope correction 2026-05-16,
  validated: the earlier "(−)→N=0,1" target was physically wrong — the (−)/N=0
  system the paper measured is parity-decoupled from N=1 and is not the cooling
  cycle. Frequencies are still cross-checked on the (−)/N=0 line, which is what
  Table II measured.)*
- A **frequency cross-check** against the measured absolute frequencies in
  arXiv:2511.06986 Table II.

**In scope:** `¹³⁸Ba¹⁹F` only; v=0 only; X N=0,1; A²Π₁/₂ J=1/2; no-M energies +
BRs + frequency validation.

**Out of scope (deferred, not foreclosed):** v=1 (Table IV); A²Π₃/₂;
non-138 isotopologues; M-resolved / field-on; the A²Π₁/₂–B²Σ⁺ perturbation
(§3.6); rigorous high-J validation beyond N=1 (§3.5).

## 2. Source & verified constants

**Paper:** M. C. Mooij *et al.* (NL-eEDM collaboration), *Precision
spectroscopy of the A²Π ← X²Σ⁺ transition in BaF*, arXiv:2511.06986v1
(10 Nov 2025). Local PDF: `~/Downloads/2511.06986v1.pdf`.

**Integrity note.** An earlier in-session "full read" operated on a
`/tmp/baf_paper.html` that does not exist on the reachable filesystem (the
extraction commands returned no captured output; the working files are gone —
a `/tmp` sandbox boundary). That coverage ledger was retracted and the paper
re-read citation-grade from the PDF via `pdf-mcp` (pages 1–10). All values
below are transcribed from **Table III** of the PDF read this session, not
from prior context.

### Table III — X²Σ⁺(v=0) and A²Π(v=0), verbatim (MHz unless noted)

| State | Constant | Value (MHz) | Value (cm⁻¹) | Provenance |
|---|---|---|---|---|
| X²Σ⁺(v=0) | T₀(X) | 0.0 | 0.0 | — |
| | B | 6473.9588 | 0.21594802 | Ryzlewicz [22] |
| | D | 0.0055250 | 1.84294e-7 | Effantin [18] |
| | H | 4.2e-10 | 1.40e-14 | Effantin [18] |
| | γ (spin-rot) | 80.984 | 0.0027013 | Ryzlewicz [22] |
| | γ_D | −0.0584 | −1.95e-6 | Ryzlewicz [22] |
| | γ_H | 0.000112 | 3.7e-9 | Ryzlewicz [22] |
| | b (Frosch-Foley) | 63.509 | 0.0021184 | Ernst [23] |
| | c (Frosch-Foley) | 8.224 | 0.0002743 | Ernst [23] |
| A²Π(v=0) | T₀,₀ band origin | 358135356.32(21) | 11946.109609(7) | **This work** |
| | B | 6347.847 | 0.2117414 | Steimle [24] |
| | D | 0.006007 | 2.0036e-7 | Effantin [18] |
| | A (spin-orbit) | 18955512.5(4) | 632.287838(13) | **This work** |
| | A_D | 0.93 | 3.1e-5 | Steimle [24] |
| | p | −7713.96 | −0.257310 | Effantin [18] |
| | q | −2.52 | −8.40e-5 | Effantin [18] |
| | p_D | −0.00699 | −2.332e-7 | Effantin [18] |
| | a (Frosch-Foley) | 26.55 | 0.0008856 | Denis [39] |
| | b+c (Frosch-Foley) | −5.54 | −0.000185 | Denis [39] |
| | d (Frosch-Foley) | 3.58 | 0.000119 | Denis [39] |

Refs: [18] Effantin 1990 Mol. Phys. 70 735; [22] Ryzlewicz & Törring 1980
Chem. Phys. 51 329; [23] Ernst, Kändler & Törring 1986 J. Chem. Phys. 84
4769; [24] Steimle *et al.* 2011 Phys. Rev. A 84 012508; [39] Denis *et al.*
2022 Phys. Rev. A 105 052811.

Only **T₀,₀ and A** for A²Π(v=0) are this-work fits (all else fixed from
prior refs; paper §V / Table III caption).

### Effective Hamiltonians (paper Eqs. 1–4, verbatim)

- Eq. (1) X²Σ⁺: `H = B·N² + γ_sr·N·S + b·(I·S) + c·I_zS_z`. Paper note:
  *"The bF hyperfine constant [Steimle 24] equals bF = b − c/3."* (This is
  Steimle's **reporting** scalar; see §3.1 — it is not the codebase mapping.)
- Eq. (3) A²Π rotation: **`H_rot = B·R² − D·[R²]²`** — explicitly R²-form,
  "congruent with pgopher"; yields Λ-doubling `(p+2q)·J` in the A²Π₁/₂ ladder,
  `q·J²` in A²Π₃/₂.
- Eq. (4) A²Π₁/₂ hyperfine: `H_hyp(²Π₁/₂) = [a − ½(b+c)]·I_z +
  ½d(e^{+2iφ}I⁻S⁻ + e^{−2iφ}I⁺S⁺)`; Denis [39], "verified consistent with
  pgopher definitions."

## 3. Convention resolution (operator-grounded)

Resolved against the **actual operators** in
`Source Code/hamiltonian_builders.py` (read this session), not parameter-file
pattern-matching. Confidence: **high** (read the source this turn).

### 3.1 X²Σ⁺ Fermi contact — RESOLVED

Builder [hamiltonian_builders.py:45-48]:
`H = Be·N² + Γ_SR·N·S + bF·(I·S) + (c/3)·√6·T²₀(I,S)`.

`'bF'` multiplies `I·S` ⟹ `'bF'` **is the Fermi-contact constant**.
Verified against **Brown & Carrington, Appendix 8.5, p.605 (book p.573),
Eq. (8.511): `bF = b + (1/3)c`** (FF operator 8.507 `aI·L+bI·S+cIzSz+…` →
B&C book form 8.510 with Fermi contact separated; also 8.516 `t0=c/3`,
8.517 `t2=(2/3)d`). This matches `b_2_bF` in `molecule_parameters.py` and the
codebase operator. The paper tabulates the Frosch-Foley **operator** constants
b, c (Eq. 1 `b·I·S + c·I_zS_z`, identical to B&C 8.507). The paper's prose
"bF = b − c/3 [Steimle 24]" is **inconsistent with B&C 8.511** — it is
Steimle's (Ref [24]) own convention, **not used** for the codebase `'bF'`.

→ `'bF'` = `b_2_bF(63.509, 8.224)` = 63.509 + 8.224/3 = **66.2503 MHz**
→ `'c'` = **8.224 MHz**
→ optionally `'b'` = `bF_2_b(66.2503, 8.224)` = 63.509 (round-trip, YbOH-style)

**Triple-confirmed (high confidence):** (i) B&C Eq. (8.511); (ii) Steimle [24]
Table III tabulates the ¹³⁸BaF ¹⁹F Fermi contact *directly* as
`bF(F) = 0.002209862 cm⁻¹ = 66.25 MHz` and `c(F) = 0.000274323 cm⁻¹ =
8.224 MHz` — i.e. Steimle's own bF **is** 66.25, and `b = bF − c/3 = 63.51`
recovers the paper's Frosch-Foley b; (iii) the measured 65.6 MHz Table II
splitting (§6). All three agree on **66.25**; `b − c/3 = 60.77` is excluded by
Steimle's own tabulated value. The arXiv:2511.06986 prose "bF = b − c/3" is
the lone outlier and is a reporting-convention statement, not the operator.

### 3.2 A²Π Λ-doubling — RESOLVED

Builder [hamiltonian_builders.py:145-167]:
`… + (p+2q)·Λ_{p2q} − q·Λ_q … ; p2q_D centrifugal via p2q_D/2·{p2q,N²}`.

`'p+2q'` multiplies the literal combined-(p+2q) operator; `'q'` enters
separately (builder applies the `−q` sign internally).

→ `'p+2q'` = p + 2q = −7713.96 + 2(−2.52) = **−7719.00 MHz** (earlier draft wrote −7718.99 — hand-add slip; −7713.96−5.04 = −7719.00; 0.01 MHz, immaterial)
→ `'q'` = **−2.52 MHz**
→ `'p2q_D'` = p_D = **−0.00699 MHz**

Sanity (paper §IV): A²Π₁/₂ J=1/2 Λ-splitting ≈ |p+2q|·(J+½) ⟹ at J=½,
≈ |−7719.00| ≈ 7719 MHz vs paper's stated **7723 MHz** (≈4 MHz residual from
p_D / higher order). |p| alone (7714) is 9 MHz off — confirms the combined
value belongs in `'p+2q'`. The `'p+2q'` and `'q'` operators are **distinct**
matrix elements (builder L147), not double-counting; for the A²Π₁/₂ J=1/2
scope the separate `'q'` term is a small higher-order / inter-ladder
correction, so a `'q'` sign error cannot corrupt the §6 hyperfine
discriminators (21.8 / 65.6 MHz). RaF A0 set `'q'`=0 (²Π₁/₂-only); BaF keeps
q=−2.52 (YbOH A000 pattern) for completeness. Cross-check: Steimle [24]
Table III gives the combined `(p+2q) = −0.25755 cm⁻¹ ≈ −7721 MHz` for
¹³⁸BaF (Steimle set q≈0), consistent with `'p+2q'` = −7719.00 MHz here
(Effantin p,q split; ~2 MHz fit-to-fit) — independently confirms the
**combined** value belongs in `'p+2q'`.

### 3.3 A²Π₁/₂ hyperfine — ONE RESIDUAL (acceptance-gated)

Boson A²Π builder hyperfine terms [hamiltonian_builders.py:146-151]:
`a·I_zL_z + h1/2·I_z + d·T²₂(I,S)` — **no bF/c term** in this builder.
Validated analog: RaF A0 (`'a'`=A‖/2-type scalar, `'h1/2'`=`'bF'`=`'c'`=0,
`'d'` with pgopher→B&C minus). RaF A0 comment: `h1/2 = a − (bF+2c/3) = A‖/2`.

Paper Eq. (4) ²Π₁/₂ diagonal scalar = `a − ½(b+c)` = 26.55 − ½(−5.54)
= **29.32 MHz**. B&C-corroborated citation-grade — **Brown & Carrington
p.845 (book p.813), verbatim:** *"The axial component of the total magnetic
hyperfine interaction, h3/2, in the ²Π3/2 component is equal to a + (b + c)/2,
where a, b and c are the Frosch and Foley constants. In the ²Π1/2 component
the axial hyperfine constant, h1/2, is equal to a − (b + c)/2."* The same page
also restates *"Note also that bF = b + (c/3)."* (second independent B&C
confirmation of §3.1). Both BaF-paper Eq. (4) forms (²Π₁/₂ `a−½(b+c)`,
²Π₃/₂ `a+½(b+c)`) match B&C exactly. So the *physical target constant* is
29.32 MHz; the residual is only whether the codebase `a·IzLz` matrix element
expects this scalar directly (Table-II-gated, not a physics question).

→ Candidate `'a'` = **29.32 MHz**, `'h1/2'` = 0, `'bF'`/`'c'` (A-state) = 0
→ `'d'` = **−3.58 MHz** (pgopher→B&C sign flip; RaF A0 precedent + Denis [39]
  states its params are pgopher-consistent)

**Residual:** whether the codebase's `a·I_zL_z` ²Π₁/₂ matrix element expects
the Eq.(4) scalar `a−½(b+c)` directly or a differently-normalized parallel
constant is a genuine cross-parametrization question. It is **not asserted** —
it is gated on the Table II empirical discriminator (§4): the A²Π₁/₂
F=0/F=1 splitting is **21.8 MHz** (≫ 0.3 MHz exp. uncertainty). Implementation
adjusts the `'a'` (and `'d'`) mapping to reproduce 21.8 MHz; if the RaF-pattern
`a−½(b+c)` does not, the matrix element `IzLz` normalization is inspected and
the mapping corrected, with the data as arbiter ("test it, the data settles
it"). STOP-flag in code comment until the 21.8 MHz check passes.

### 3.4 A²Π rotational constant Be — RESOLVED (R²-form, no −2Λ²D)

Builder [hamiltonian_builders.py:145,165]: `Be·elements['N^2'] − D·(N²)²`,
where `elements['N^2']` is the **R²-form** operator (per commit `ef06b86` /
the authoritative RaF A0 inline comment: "B&C eq 9.138 N²-op minus Lam²·I").
Paper Eq. (3) is R²-form. RaF A0's `−2Λ²D` exists *only* because RaF's source
paper reported B in N²-form; **this paper's B is already R²-form**.

→ BaF A `'Be'` = **6347.847 MHz directly, NO −2Λ²D**.
Inline-comment must state: R²-form per arXiv:2511.06986 Eq. (3); deliberately
diverges from RaF A0 (RaF source N²-form → needs −2Λ²D; this source R²-form →
none). Applying −2Λ²D here would *introduce* a ∝D·J² high-J error.

### 3.5 High-J (deferred, not foreclosed)

`Be=6347.847` is correct at **all J**, not a low-J compromise. The paper's
rotational + centrifugal set is fixed from Effantin [18] (rotational levels
**J>100**) and reproduces Table I (SR-branch to J=15/2) at ≤1 MHz in the
R²-form framework the codebase engine implements. High-J *frequency* accuracy
is reachable now by transcribing the full Table III centrifugal set (D, H,
A_D, p_D, γ_D, γ_H) — cheap, no rework. Include them in the dict even though
N≤1 does not exercise them.

### 3.6 A²Π₁/₂–B²Σ⁺ perturbation (BR caveat)

Paper §III: measured QR/QQ branch intensities deviate from the single-state
pgopher simulation, attributed to A²Π₁/₂–B²Σ⁺ mixing (B²Σ⁺ at T₀=14040.163
cm⁻¹ [Effantin 18]); credited to a referee in Acknowledgements. The pipeline's
no-M BRs are idealized single-state values. **Expectation:** BR validation =
"sums to 1 + physically sensible structure", **not** quantitative agreement
with measured intensities. Frequencies are unaffected (Table I/II reproduced
to 1 MHz with the single-state Hamiltonian). Modeling A–B mixing is a separate
deferred task affecting BRs only.

## 4. Proposed `molecule_parameters.py` entries

`¹³⁸Ba¹⁹F` = `¹³⁸Ba` (I=0) + `¹⁹F` (I=1/2): structurally identical to RaF
boson `²²⁶Ra`(I=0)`¹⁹F`(I=1/2). Mirror the **RaF boson X0/A0** template
key-for-key. (`Molecule_Library` routes BaF through the YbOH backend — known
wart, CLAUDE.md; extend, do not refactor.)

Unit convention (per user): **MHz directly** from arXiv:2511.06986 Table III's
MHz column (cm⁻¹ in comment) — *not* the RaF `<cm⁻¹>*c` form. `muE` keeps the
codebase's intrinsic `<D>*0.503412` (Debye → MHz/(V·cm⁻¹)) form.

```python
molecules['BaF']['boson']['X0'] = {
    'Be': 6473.9588,        # Ryzlewicz1980 [22]; 138Ba19F; 0.21594802 cm^-1
    'Gamma_SR': 80.984,     # Ryzlewicz1980 [22]; 0.0027013 cm^-1
    'bF': b_2_bF(63.509, 8.224),  # =66.2503; Frosch-Foley b,c Ernst1986 [23]; bF=b+c/3 (B&C Eq.8.511; Steimle[24] tabulates bF(F)=66.25 directly)
    'c': 8.224,             # Ernst1986 [23]; 0.0002743 cm^-1
    'D': 0.0055250,         # Effantin1990 [18]; 1.84294e-7 cm^-1
    'muE': 3.170*0.503412,  # 3.170(3) D, Steimle2011 [24] Table V; Debye->MHz/(V/cm)
    # high-J only (§3.5; verify X2S builder consumes before adding): H=4.2e-10, Gamma_D=-0.0584, Gamma_H=0.000112 MHz
    }
molecules['BaF']['boson']['A0'] = {
    'Be': 6347.847,         # Steimle2011 [24] B(R^2); 0.2117414 cm^-1. arXiv:2511.06986 Eq.(3) R^2-form -> NO -2*Lam^2*D (unlike RaF A0 N^2-form source)
    'ASO': 18955512.5,      # This work [arXiv:2511.06986]; 632.287838 cm^-1
    'h1/2': 0,
    'a': 26.55 - 0.5*(-5.54),  # =29.32; Eq.(4)/B&C p.845 h1/2=a-(b+c)/2; Denis2022 [39]. STOP: gate on Table II 21.8 MHz A-hf splitting
    'bF': 0,
    'c': 0,
    'd': -3.58,             # Denis2022 [39]; pgopher->B+C needs minus (RaF A0 precedent)
    'p+2q': -7713.96 + 2*(-2.52),  # =-7719.00 (earlier draft -7718.99: hand-add slip, -7713.96-5.04=-7719.00; 0.01 MHz immaterial); Effantin1990 [18] combined (op-grounded builders:147; Steimle[24] (p+2q)=-7721 cross-check)
    'q': -2.52,             # Effantin1990 [18]; -8.40e-5 cm^-1
    'p2q_D': -0.00699,      # Effantin1990 [18] p_D; -2.332e-7 cm^-1
    'D': 0.006007,          # Effantin1990 [18]; 2.0036e-7 cm^-1
    'A_D': 0.93,            # Steimle2011 [24]; include iff A2Pi builder consumes it (Open #5)
    'muE': 1.50*0.503412,   # 1.50(2) D A2Pi1/2, Steimle2011 [24] Table V; Debye->MHz/(V/cm)
    'g_S': 2.0023,
    'Origin': <T0,0-based>, # see §5; absolute origin carries known convention offset
    }
```

`<T0,0-based>` is the one explicit unresolved input (§5), not a placeholder to
ship. `muE` values are now resolved (Steimle [24] Table V).

## 5. The `Origin` key

RaF A0 `'Origin'` = `Π₁/₂ electronic origin + 0·ASO + Be_X/c` ("B offset due
to code being R²"). The paper itself flags (§V) a **0.21 cm⁻¹** deviation of
T₀,₀ from Steimle [24] "attributed to different definitions in the effective
Hamiltonian" — i.e. the absolute origin is convention-sensitive at ~B. Build
`'Origin'` on the RaF A0 template using arXiv:2511.06986 T₀,₀; **document the
known absolute offset**; validate on **relative** splittings (§6), with a
generous absolute-origin tolerance. Exact construction resolved during
implementation against the RaF A0 formula + the Table II absolute frequency.

## 6. Verification plan

**Target (Table II, exact, MHz):** X²Σ⁺ v=0, N=0, J=1/2 → A²Π₁/₂ v=0,
J=1/2 (−):

| Lower F | Upper F | Frequency | Unc. |
|---|---|---|---|
| 1 | 0 | 348666402.6 | 0.3 |
| 1 | 1 | 348666424.4 | 0.3 |
| 0 | 1 | 348666490.0 | 0.3 |

**Empirical discriminators (decisive, ≫ uncertainty):**
- A²Π₁/₂(J=1/2) F=0↔1 splitting = 424.4 − 402.6 = **21.8 MHz** → tests the
  §3.3 A-hyperfine mapping (`'a'`, `'d'`).
- X²Σ⁺(N=0,J=1/2) F=0↔1 splitting = 490.0 − 424.4 = **65.6 MHz** → tests the
  §3.1 X Fermi-contact mapping (`'bF'`).

**Acceptance:**
1. No-M pipeline reproduces both splittings to **≤1 MHz** (paper's stated
   absolute accuracy; exp. uncertainties 0.3 MHz).
2. The three Table II (−)-line components reproduced **line-by-line, matched
   F_X→F_A** (F1→F0, F1→F1, F0→F1; F0→F0 is E1-forbidden), each ≤ ~1.5 MHz
   (paper abs. acc. ~1 MHz, fixed non-refit literature constants). *Validated
   2026-05-16: +0.6 / +0.8 / +1.4 MHz.* (Supersedes the old "centroid /
   documented-offset" check — a centroid that includes the forbidden F0→F0
   line is not apples-to-apples with Table II's 3 components.)
3. Cooling-line BR A²Π₁/₂(v=0,J=1/2,+) → X: **parity-closed** — Σ→X N=1 = 1.0,
   Σ→X N=0 ≤ 1e-3 (rotational closure, the defining property), normalized to
   1, spread over the X N=1 hyperfine manifold. Magnitudes are idealized
   single-state values — no quantitative intensity match expected (§3.6).
4. Sanity: computed A²Π₁/₂ J=1/2 Λ-splitting ≈ 7723 MHz (paper §IV).
5. Verification gate (CLAUDE.md): tutorial notebook runs end-to-end on a
   fresh kernel — `conda run -n Structure jupyter execute "Jupyter
   Notebooks/RaF_Calcs_Tutorial.ipynb"`.

If (1) fails for the A-state 21.8 MHz, the §3.3 residual is live: inspect the
`IzLz` ²Π₁/₂ matrix-element normalization and correct the `'a'`/`'d'` mapping;
do not tune blindly.

## 7. Open items / risks

1. **§3.3 A²Π₁/₂ `'a'` scalar** — RaF-pattern candidate `a−½(b+c)=29.32`;
   confirmed only by the 21.8 MHz Table II check. *Highest-risk item.*
2. **`muE` — RESOLVED.** Steimle [24] Table V (read citation-grade):
   μ(X²Σ⁺)=**3.170(3) D**, μ(A²Π₁/₂)=**1.50(2) D** (A²Π₃/₂=1.31(2) D, out of
   scope). Entered as `<D>*0.503412`. Not exercised by the freq/BR validation;
   required for field-on use.
3. **`'Origin'` construction** — convention offset documented (§5); resolve
   exact formula vs RaF A0 during implementation.
4. **YbOH-backend routing** — confirm BaF dispatch keys
   (`molecule_library_class.py`) follow the RaF boson path.
5. **`A_D`, `p2q_D` consumption** — confirm the boson A²Π builder consumes
   these keys (builder shows `p2q_D` and `D`; verify `A_D`).
6. **Unit-entry style — RESOLVED (user).** Use the arXiv:2511.06986 Table III
   **MHz column directly**, cm⁻¹ in comment; `muE` keeps `<D>*0.503412`.
   Diverges from RaF's `<cm⁻¹>*c` style by explicit user choice.

## 8. Verified-vs-asserted ledger

| Claim | Basis | Confidence |
|---|---|---|
| Table III values (§2) | PDF read this session | high |
| X `'bF'`=b+c/3=66.25 | B&C Eq.8.511 + B&C p.845 + Steimle[24] TblIII bF(F)=66.25 + builders.py:46/187 + 65.6 MHz Tbl II | high |
| μ(X)=3.170 D, μ(A²Π₁/₂)=1.50 D | Steimle[24] Table V (read this session) | high |
| A `'p+2q'`=combined | hamiltonian_builders.py:147 | high |
| A `'Be'` no −2Λ²D | Eq.(3) R²-form + builder L145/165 + ef06b86 | high |
| A `'a'` physics value 29.32 | B&C p.845 h1/2=a−(b+c)/2 + Eq.(4) | high |
| A `'a'` codebase-mapping exact | RaF pattern; Table II 21.8 MHz gates | moderate |
| `'d'`=−3.58 sign | RaF A0 comment + Denis pgopher-consistent | moderate |
| BR ≠ measured intensities | paper §III A–B anomaly | high |
| High-J improvable later | paper p.6/8 (J>100 fixed set) | high |
