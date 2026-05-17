# The N²↔R² formalism converter

> *"In my experience, there's no such thing as luck."* — Obi-Wan Kenobi

Every surprising number in this story — a +11.6 MHz line shift, an entire
fermion spectrum of `nan`, a rotational constant 17 MHz off the paper — turned
out to have a structural cause rather than a chance one. This note records what
the rotational-constant bookkeeping did before, what we changed, and what it
does now. We begin with the physics of the two conventions, narrate the
failure chain that motivated a converter, and close with the bidirectional
function the code now exposes.

## Why a convention even enters

The RaX program needs the A²Π₁/₂–X²Σ⁺ level structure of RaF correct in an
absolute sense — line positions feed the laser-cooling scheme and the
PT-violation sensitivity downstream, so a constant offset buried in the term
origin is not cosmetic. The subtlety is that a ²Π effective Hamiltonian can be
written with the rotational operator in two conventions, and the fitted
constants are *different numbers* in each. Modern spectroscopy, and the PGopher
default used by the Udrescu *et al.* RaF measurement,[^paper] writes the
rotational energy as `B·N²`, with **N** the rotational angular momentum that
still carries the on-axis electronic projection Λ. Our engine instead
diagonalizes the R²-form operator `B·(N²_op − Λ²)`, with the Λ² piece removed —
machine-verified equal to Brown & Carrington Eq. 9.138. The two operators
describe the same molecule; their constants do not coincide.

Before writing any formula, the intuition is worth stating: stripping Λ² is a
*constant* energy shift, not a change in how energy scales with rotation. So
the bulk of the N²↔R² difference must land on the term origin, of order
Λ²·B ≈ 0.19 cm⁻¹ ≈ 5.7 GHz for RaF A²Π — enormous — while the rotational
coefficient B itself moves only by the tiny quartic cross-term, 2Λ²D ≈ 0.008
MHz, utterly negligible. Brown & Carrington Table 7.2 makes this precise. We
quote it verbatim, quartic-truncated (the codebase carries no sextic *H*),
expressing the N² constants in terms of the R² ones:[^bc]

| row | B&C Table 7.2 (quartic), N² in terms of R² | the inverse used by the engine |
|---|---|---|
| B   | `B(N²) = B(R²) + 2Λ²D`              | `B(R²) = B(N²) − 2Λ²D` |
| D, X_D | identities at quartic order      | unchanged |
| X (p+2q, Γ_SR, q_lD) | `X(N²) = X(R²) − Λ²X_D` | `X(R²) = X(N²) + Λ²X_D` |
| G (Origin) | `G(N²) = G(R²) − Λ²B(R²) − Λ⁴D` | `G(R²) = G(N²) + Λ²B(N²) − Λ⁴D` |

The takeaway the table is meant to convey: only the term origin G carries a
large convention shift; B, D and the centrifugal partners barely move, and the
transformation is exactly invertible — it is, after all, just arithmetic.

## What it did before, and how it failed

The conversion used to be done by hand, inside the constant literals in
`molecule_parameters.py` — `'Be': 5743.96 - 2*1.4e-7*c`, `'Origin': 13284.427
+ 5755.56/c`, and so on. Hand arithmetic on a sign-sensitive transformation
fails the way hand arithmetic always does. First the B-row sign was wrong
(`+2Λ²D` where Brown & Carrington requires `−2Λ²D`); we caught and fixed that
in `ef06b86`. Then the term-origin shift was found to use the *X-state*
rotational constant, B_X = 5755.56 MHz, where the G-row demands the *A-state*
constant Λ²B_A — a clean systematic of +11.6 MHz on every ²²⁶RaF X–A line, and
+26.5 MHz on the ²²⁵RaF lines, sitting live in the code. The pattern is the
one this note opened with: not bad luck, a structural error in a step that
should never have been performed by hand.

The fix was to remove the hand step. We migrated the RaF A0 dicts to carry the
*raw paper N² values* tagged `formalism:'N2'`, `Lambda:1`, and built a
converter that applies Table 7.2 at load time. Verifying it surfaced two more
structural causes, both worth recording because neither was what it first
looked like. The fermion RaF spectrum came back entirely `nan`; the obvious
read was a broken hyperfine matrix element, and the phase
`(−1)^(2F₀−2M₀+J₁+i_H+F₁₀+J₀−P₀+S−Σ₀)` in the Frosch–Foley *d*-term was indeed
evaluating a negative base to a half-integer power. But instrumenting the
exponent showed the matrix element was correct in its intended domain — the
regression harness had simply omitted `I_nuclei`, defaulting the Ra spin to
zero and collapsing the two-nuclear-spin recoupling. The matrix element was
never the problem; the test configuration was. Separately, the external check
against the published constants showed the boson A²Π₁/₂ rotational constant the
code had carried for years — 5743.96 MHz, labelled "B(N²) paper value" — sitting
≈17.5 MHz above the Udrescu *et al.* value of 0.191015(5)[15] cm⁻¹. We traced
this far enough to be sure it is not a convention artifact (the N²↔R² B-row term
is the 0.008 MHz one, not 17 MHz) but a value from some other source; the prior
number's provenance is unclear and very plausibly just wrong, so we set it to
the published value and moved on without dramatizing the downstream
implications.

## What it does now

The converter is one pure function, `convert_formalism(params, c_cm) → dict`,
and it is bidirectional. Direction follows the data: a dict tagged `'N2'` is
converted to R² and returned tagged `'R2'`; a dict tagged `'R2'` is converted
the other way and returned tagged `'N2'`. The tag therefore always describes
the numbers actually in the dict, which means any caller can ask a parameter
set which convention it is in and the converter can always work out which way
to go. `Lambda` is read, never consumed, and left in place; a dict with no
`formalism` tag is returned untouched, which is exactly what keeps every legacy
(untagged, R²) entry byte-identical through the converter. It converts every
parameter Table 7.2 defines — the B-row on `Be`, the G-row on `Origin`, and the
generic-X row on every verified centrifugal pair present.[^partners] An
N²→R²→N² round trip returns the input to within 10⁻⁹, as it must for a
transformation that is pure arithmetic.

The engine wants R², so the two load-time call sites — `get_molecule_params`
and the user-dict path in `Energy_Levels` — convert only when the declared
formalism is `'N2'`, leaving R² and untagged sets alone. The `formalism` and
`Lambda` keys now ride along in the returned dict rather than being stripped;
the Hamiltonian builders read parameters by targeted key access, so the extra
keys are inert, and the tutorial notebook executing end-to-end on a fresh
kernel confirms this in practice rather than in principle.

Against the primary source the result is what we want: the migrated boson RaF
A0 reproduces the Udrescu *et al.* p+2q = −0.41071 cm⁻¹ and band origin
T_Π1/2,0 = 13284.427(1)[20] cm⁻¹ exactly, to Δ = 0!  The term-origin G-row, the
step that carried the +11.6 MHz hand-bug, is now validated against published
spectroscopy rather than merely against itself. The same converter, run
backward, is what lets us pull an engine-side R² set back to the paper's N²
convention for that comparison in the first place — which is the use the
bidirectional form was built for, and the natural lever for checking any future
state against a PGopher fit.

[^paper]: O. O. Udrescu, S. M. Udrescu, R. F. Garcia Ruiz, *et al.*, *Nature
Physics* (2024), DOI 10.1038/s41567-023-02296-w. ²²⁶Ra¹⁹F (I_Ra = 0 ⇒ the boson
entry); the paper defers ²²³,²²⁵RaF to future work, so it does not constrain the
fermion constants. Constants verbatim in
`docs/superpowers/plans/2026-05-16-n2-r2-formalism-converter-natphys-oracle.md`.

[^bc]: Brown & Carrington, *Rotational Spectroscopy of Diatomic Molecules*,
Table 7.2 (PDF p. 376; book p. 344). Origin is stored in cm⁻¹ while B and D are
in MHz, so the G-row carries a 1/c factor, c = 29979.2458 MHz/cm⁻¹.

[^partners]: The verified pairs in `CENTRIFUGAL_PARTNERS` are `p+2q`/`p2q_D`
(B&C eq. 7.190), `Gamma_SR`/`Gamma_D` (eq. 7.189), `q_lD`/`q_lD_D` (eq. 7.190),
and `ASO`/`A_D` (spin–orbit, eq. 7.187 — verified against B&C 2026-05-16). Each
enters the effective Hamiltonian as (X + X_D·N²) and is "any molecular parameter
other than G, B or D", so each obeys Table 7.2's generic-X row, X(N²) =
X(R²) − Λ²·X_D (B&C Table 7.2, PDF p. 376 / book p. 344). The 2026-05-15 design
spec describes the earlier one-directional, metadata-stripping contract; it now
carries a superseded banner pointing to the bidirectional module documented here.
