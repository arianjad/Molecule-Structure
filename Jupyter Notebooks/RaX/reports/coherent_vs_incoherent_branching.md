# Coherent vs. incoherent polarization sums in the TDM / branching code

**Status:** settled, 2026-05-15. Companion regression test: `validate_noM_zerofield.py`.
**TL;DR:** `branching_ratios` builds a *coherent* sum `Σ_p (−1)^p TDM_p` and squares it, whereas the
correct line strength (thesis Eq. 3.3) is the *incoherent* `Σ_p |TDM_p|²`. **For branching ratios the two
are identically equal** — not approximately — because of a one-line selection rule. They diverge only
for problems that coherently combine several p into the *same* final state (driving, CPT, interference,
non-axial fields). The code is correct for what it is used for.

---

## 1. The two quantities

Let `TDM_p(g,e) = ⟨g | T¹_p(d) | e⟩` be the lab-frame spherical component `p ∈ {−1,0,+1}` of the
transition-dipole operator between (eigenvector-contracted) states `g` (lower) and `e` (upper).

- **Code path.** `branching_ratios` / `XA_branching_ratios` use the operator
  `TransitionDipole_even_aBJ = Σ_p (−1)^p · TDM_p_even_aBJ(p, …)` and then square the eigenbasis
  contraction: `BR ∝ (G_evecs @ TDM @ E_evecs.T)²`. Call this **`D_coh = |Σ_p (−1)^p TDM_p|²`**.
- **Thesis definition.** The line strength / partial width is the *incoherent* sum
  (Springer thesis Eq. 3.1–3.3): `S = Σ_{M₀,M₁,p} |⟨g,i|T¹_p(d)|e,j⟩|²`, with
  `γ̃_{g,e} = S / (2F_e+1)`. Call the per-(g,e) piece **`S_inc = Σ_p |TDM_p|²`**.

The question: does using `D_coh` instead of `S_inc` corrupt branching ratios?

---

## 2. The single-p selection rule (the whole reason)

Wigner–Eckart in the lab frame (thesis Eq. A.35), for total angular momenta `F₀` (lower), `F₁` (upper):

```
⟨g; F₀ M₀ | T¹_p(d) | e; F₁ M₁⟩ = (−1)^{F₀−M₀} · ( F₀  1  F₁ ; −M₀  p  M₁ ) · ⟨g;F₀‖T¹(d)‖e;F₁⟩
```

The 3-j symbol is nonzero **only** when its bottom row sums to zero:

```
−M₀ + p + M₁ = 0      ⟹      p = M₀ − M₁ .
```

So between a *definite* pair of magnetic sublevels `(M₀, M₁)` there is **exactly one** allowed `p`
(namely `p* = M₀ − M₁`); for every other `p` the matrix element vanishes identically. In the `aBJ`
basis the reduced element and the 3-j/6-j factors are **real**, so `TDM_p` is a real signed number.

---

## 3. Theorem: for definite-M states, `D_coh ≡ S_inc`

Fix `(M₀, M₁)`. By §2, `TDM_p = 0` for all `p ≠ p*`, and `TDM_{p*} = A ∈ ℝ`.

- Coherent: `Σ_p (−1)^p TDM_p = (−1)^{p*} A`  ⟹  `D_coh = |(−1)^{p*} A|² = A²`.
- Incoherent: `S_inc = Σ_p |TDM_p|² = |A|² = A²`.

`D_coh = S_inc`. ∎

The `(−1)^p` phase is annihilated by the modulus square, and **every cross term `p ≠ p′` is zero
because at most one `TDM_p` is nonzero**. There is no approximation: the equality is exact whenever
both bra and ket are definite-M (pure `|…F M⟩`) — i.e. at zero field, or in any axial field, where
`M` stays a good quantum number (thesis p. 284: an axial field has only `p=0` components and does
**not** mix `M`; only off-axis `p≠0` components mix `M`).

### Corollary: why `_noM` is exact (the analytic M-collapse)

3-j orthogonality (B&C App. C; Zare):

```
Σ_{M₀, p} ( F₀  1  F₁ ; −M₀  p  M₁ )²  =  1 / (2F₁+1)      (sum over M₀ and p, fixed M₁).
```

Hence the orientation-summed strength factorizes analytically:

```
S = Σ_{M₀,p} |⟨g;F₀M₀|T¹_p|e;F₁M₁⟩|² = |⟨g;F₀‖T¹(d)‖e;F₁⟩|² · Σ_{M₀,p}(3j)²
  = |reduced|² / (2F₁+1)  =  |TransitionDipole_even_aBJ_noM|² ,
```

since `_noM` is built as `reduced / √(2F₁+1)`. No field, no eigenvectors, no labeling needed.
This is exactly thesis Eq. 3.3, and is what `validate_noM_zerofield.py` checks.

---

## 4. Physical reading

Branching = **spontaneous emission into all directions from a definite excited state**: you sum the
*rates* `|⟨g|T¹_p|e⟩|²` over emitted-photon polarization `p` and over final sublevels `M₀` — distinct
emitted photons carry no mutual phase, so the sum is incoherent by construction (thesis Eq. 3.1–3.3).
The single-p rule then makes the code's coherent square numerically identical to that incoherent sum.

---

## 5. Empirical evidence (rederived, RaF X²Σ⁺ – A²Π₁/₂, 226)

### Campaign 1 — `D_coh` (`branching_ratios`) vs explicit `S_inc = Σ_p |Calculate_TDMs(p)|²`

| Field | raw `max|D_coh − S_inc|` | normalized-BR max diff | median `S_inc/D_coh` |
|---|---|---|---|
| zero (B≈0)            | 7.9e-3 | 1.2e-2 | **1.0000** |
| axial B = 100 G       | 6e-16  | 4e-16  | **1.0000** |
| axial E = 500 V/cm    | 4.9e-3 | 1.4e-2 | **1.0000** |
| axial E + B           | 2e-9   | 2e-9   | **1.0000** |

**Reading.** The median ratio is *exactly* 1.0000 at every field — a genuine convention error would
bias the median, not merely scatter the max. Once the M-degeneracy is lifted (B = 100 G; Zeeman ∝ M),
agreement is **machine precision (6e-16)**. The ~1 % wobble survives *only* at exactly zero field (and
E-only, which does not split M): it is **degenerate-eigenvector basis ambiguity** — the solver returns
an arbitrary orthonormal basis within each exactly-degenerate M-subspace, and the two code paths land
on different rotations of it, so *per-eigenstate* elements differ while all *physical*
(degeneracy-summed) quantities agree. It is not the coherent/incoherent convention. A small quantizing
field, or summing degenerate partners, removes it entirely.

### Campaign 2 — independent zero-field validation (`validate_noM_zerofield.py`)

Independent from-scratch Wigner reimplementation of the thesis chain (Eqs. A.35 / A.39 / A.54;
branching 3.1–3.4), no fields / eigenvectors / labeling:

| Criterion | Result |
|---|---|
| **A** independent ME ≡ codebase `TDM_p_even_aBJ` | `1.1e-16` |
| **B** `S_ref / |_noM|² = 2F_e+1` exactly | `1.8e-15` (ratios 1, 3, 5 for F_e = 0,1,2) |
| **C** `r̃_noM ≡ r̃_ref`, closure `Σ_lower r̃ = 1` | `2e-16`, `1.0000000000` ∀ excited |

---

## 6. Where coherent vs. incoherent **does** matter

The equality is specific to incoherent decay. Cross terms are real and the distinction is essential
whenever multiple `p` pathways terminate on the **same** state and you square the *summed* amplitude:

- Polarization-resolved **driving** / Rabi coupling (`ε̂·d = Σ_p c_p(ε̂) TDM_p`, *then* square).
  x̂/ŷ light are coherent `(d_{+1} ∓ d_{−1})` superpositions; ŷ carries an `i` (complex amplitude).
- σ⁺/σ⁻ **interference** in transverse fields; **non-axial** fields where eigenvectors mix `M`
  (`p≠0`), so one eigenstate-pair element sums several `p` with `(−1)^p` phases — cross terms ≠ 0.
- **CPT / dark states**, two-photon EOM interference — the thesis Ch. 5 physics. Here coherence is
  the entire effect.

For those, keep the per-`p` complex amplitudes (`Calculate_TDMs` / `TDM_p_builders`) and combine
*before* squaring. Note `select_dipole` in the X–A notebooks does the coherent x̂ combination
correctly but carries a spurious `/√3` on the `'all'` path and lacks a `'y'` branch — irrelevant to
branching ratios, relevant if/when you use it for absolute strengths or coherent multi-`p` drives.

---

## 7. Practical guidance

- `branching_ratios` / `XA_branching_ratios` are **correct** for branching ratios at zero or any
  axial field. Add a small quantizing `Bz` (standard practice) to remove the zero-field
  degenerate-eigenvector ambiguity; or use the `_noM` path, which is exact at true zero field.
- `TransitionDipole_even_aBJ_noM` is a validated closed-form implementation of thesis Eq. 3.3.
- Re-run `validate_noM_zerofield.py` after any edit to the TDM / branching path — the project gate
  (`RaF_Calcs_Tutorial.ipynb`) exercises none of it.
- Do **not** reuse "coherent ≡ incoherent" outside incoherent decay; for driving/CPT/interference it
  is false and the cross terms are the physics.

**References:** Springer thesis Eqs. 3.1–3.6 (§3.2.4.1), A.35 (Wigner-Eckart), A.39 (spectator
theorem), A.47/A.54 (case-(a) reduced-D / TDM); Brown & Carrington App. C (3-j orthogonality).
Code: `Source Code/Energy_Levels.py` (`branching_ratios`, `Calculate_TDMs`),
`Source Code/matrix_elements.py` (`TDM_p_even_aBJ`, `TransitionDipole_even_aBJ_noM`).
