# Nat. Phys. 2024 (s41567-023-02296-w) — RaF N² validation oracle

**Source:** O. O. Udrescu, S. M. Udrescu, R. F. Garcia Ruiz, *et al.*,
"Precision spectroscopy and laser cooling scheme of a radium-containing
molecule," **Nature Physics** (2024), DOI **10.1038/s41567-023-02296-w**.

Constants are the **PGopher / N² formulation** (paper p.6: "The fit was
performed using the publicly available software PGOPHER [29]"). Units **cm⁻¹**.
Isotopologue: **²²⁶Ra¹⁹F** (paper p.8, Conclusion: "...of ²²⁶Ra¹⁹F"; ²²⁶Ra has
I=0 ⇒ the **boson** entry). The paper explicitly defers ²²³,²²⁵RaF to future
work (p.8), so it does **not** constrain the `fermion` entry.

Verbatim transcription of **Table I** ("Fitted rotational parameters"),
PDF p.7 of the Zotero copy (`MSQLIKLJ`). Uncertainties: 1σ statistical in
(round) brackets, systematic in [square] brackets, in units of the last digit.
Footnotes: ᵇ Ref.[13], ᶜ Ref.[31], ᵈ Ref.[30] (theory); ᵃ theory reports
equilibrium params.

## X²Σ⁺   [cite: Table I, p.7]
| param | ν = 0 (cm⁻¹) | ν = 1 (cm⁻¹) | theory |
|---|---|---|---|
| B″_ν | **0.191985(5)[15]** | 0.19092(4)[6] | 0.192ᵇ / 0.1909ᶜ |
| 10⁷·D″_ν | 1.40(5)[10] | 1.2(3)[4] | — |
| γ_ν (spin-rot) | 0.00585(3)[7] | 0.00581(5)[15] | 0.006ᵈ |

Also (p.7): the two lowest X²Σ⁺ levels are separated by **263(4) MHz**.

## A²Π₁/₂   [cite: Table I, p.7]
| param | ν = 0 (cm⁻¹) | ν = 1 (cm⁻¹) | theory |
|---|---|---|---|
| B′_ν | **0.191015(5)[15]** | 0.18997(4)[6] | 0.192ᵇ / 0.1902ᶜ |
| 10⁷·D′_ν | 1.40(5)[10] | 1.5(3)[4] | — |
| p_ν (Λ-doubling) | **−0.41071(3)[7]** | −0.40978(10)[20] | — |
| 10⁷·p_D,ν | 1.9(2)[5] | 4.4(20)[15] | — |
| T_Π1/2,ν (band origin) | **13284.427(1)[20]** | 13278.316(1)[20] | 13300ᶜ |

`T_Π1/2,0 = 13284.427 cm⁻¹` is the **physical A²Π₁/₂–X²Σ⁺ band origin** (the
paper uses it verbatim as the x-axis shift in Figs. 2,3 and Methods, pp.6,9,10,13,14).

## Representative N=1 X–A line
The paper presents spectra as **figures** (Figs. 2, 4); no single N=1 (P(1)/Q(1)/
R(0)) transition frequency is **individually tabulated**. Per the plan's
acceptance criterion #3 fallback ("if only constants are tabulated, the
corrected `Origin` matches Tₑ within rounding"), validation is the
constants + band-origin comparison below.

## Codebase comparison (migrated `RaF/boson/A0` via `get_molecule_params`)
Converter (`formalism.py`, B&C Table 7.2, Λ=1): inverting R²→N² for the
paper comparison —
`Be(N²)=Be(R²)+2Λ²D`, `p+2q(N²)=p+2q(R²)−Λ²·p2q_D`,
`Origin_phys(N²)=Origin(R²)−Λ²·Be(N²)/c+Λ⁴·D/c`.

| quantity | paper (²²⁶RaF, N²) | migrated boson A0 → N² | result |
|---|---|---|---|
| Be (A²Π₁/₂, v=0) | 0.191015 cm⁻¹ | 0.191015 cm⁻¹ | ✅ exact |
| p+2q | −0.41071 cm⁻¹ | −0.41071 cm⁻¹ | ✅ exact |
| Origin / T_Π1/2,0 | 13284.427 cm⁻¹ | 13284.427 cm⁻¹ | ✅ exact |

## Provenance correction (Task 6 finding)
Pre-Task-6, `RaF/boson/A0['Be']` was **5743.96 MHz** (= 0.191598 cm⁻¹),
commented "B(N²) paper value" throughout the git history (bacb463 → ef06b86 →
b7c8548) and the design spec. It is **~17.5 MHz / 0.000583 cm⁻¹ high** vs this
paper's A²Π₁/₂ B′₀ = 0.191015 cm⁻¹ and did **not** originate here. This is
*not* an N²↔R² formalism artifact: B&C Table 7.2's B-row N²↔R² term is only
`−2Λ²D ≈ 0.008 MHz` (three orders below the gap; converter-handled, Task 5
byte-exact). It was a wrong-source value. Corrected to `0.191015*c` (user-
greenlit, 2026-05-16); derived `g_lp`/`g_l` and the Task 4/5 oracles updated in
lock-step. p+2q & Origin were entered correctly from this paper and were exact
already — they validate the converter's X-row and G-row (the Origin G-row fix
of the +11.6/+26.5 MHz hand-bug) against published spectroscopy.

`fermion/A0['Be']=5729.03` (²²⁵RaF) is **unchanged** — not measured by this
paper; its source is theory/isotope-scaling and is out of scope for this oracle.

Historical design/plan docs (`docs/superpowers/specs|plans/*`) still show the
superseded `5743.96` in worked examples and the historical "+11.6 MHz" bug
analysis (which was correct *relative to the then-current Be*); left as the
historical record by intent (not rewritten).
