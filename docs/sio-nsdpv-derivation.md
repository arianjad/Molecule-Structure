# ²⁹SiO⁺ NSD-PV anapole operator C = (n̂×S)·I/I — case bβJ derivation note

**Task 4, step 1 (amended per spectroscopy-reviewer audit, 2026-06-11).** Derivation only;
no code in this note. The target is the matrix element of the dimensionless NSD-PV operator
**C = (n̂×S)·I/I** (literal `/I`: divide by the nuclear-spin quantum number; ×2 for I=1/2,
see §2.1) in the codebase's
case-bβJ basis `|K=0, N, J, F, M⟩` (S=1/2, I=1/2 in the last slot), in a form directly
implementable as `NSDPV_bBJ(K0,N0,J0,F0,M0,K1,N1,J1,F1,M1,S,I)` matching the codebase's
existing signatures, index convention, and Wigner-symbol style.

Every formula carries a `(source, eq, page)` tag or is derived in ≤3 visible steps from
tagged starting points. Numbers reported as **"verified"** were computed in the `Structure`
conda env against an independent first-principles decoupled-basis construction (see §4).

Conventions used throughout (per `docs/sio-conventions.md`, "B&C audit results"):
- **Index 1 = ket = primed; index 0 = bra = unprimed.**
- B&C = Brown & Carrington, *Rotational Spectroscopy of Diatomic Molecules*, Zotero CKZKCGXY.
  **PDF page = printed page + 32.** Pages cited below are **PDF pages**.
- Karthein = PRL 133, 033003 (2024), local copy arXiv:2310.11192v1, Zotero LYF85P8C.
- DeMille = PRL 100, 023003 (2008); local arXiv:0708.2925v1, `/tmp/demille2008_nsdpv.pdf`.

---

## 1. Operator definition and sign/phase convention

**Karthein** (p. 2, main text, "Effective Hamiltonian and Electroweak Properties"):

> An effective Hamiltonian acting only within the subspace of rotational and hyperfine
> levels … is given by **H_eff = κ′ W_A C**, where W_A … includes the expectation value of
> H_PV over the electronic wave function … **C = (n×S)·I / I** contains the angular momentum
> dependence of H_eff and its matrix elements can be calculated analytically using angular
> momentum algebra [21]. (Karthein's ref [21] = **Flambaum & Khriplovich, Phys. Lett. A 110,
> 121 (1985)** — the original NSD-PV effective-operator paper, which is also DeMille's ref [9];
> **DeMille 2008 is Karthein's ref [20]**, not [21]. Citation corrected per review.)

The PV matrix element is defined (Karthein p. 2):

> **iW(m′_N, m′_I, m_N, m_I) ≡ κ′ W_A ⟨Ψ⁻_↓(m′_N,m′_I)| C |Ψ⁺_↑(m_N,m_I)⟩**

and the two-level effective Hamiltonian (Karthein **Eq. (B1)**, p. 8):

> H± = [[ 0 , iW + Ω_R sin(ω_ext t) ], [ −iW + Ω_R sin(ω_ext t) , Δ ]],
> where **iW is the imaginary parity violating matrix element**.

**DeMille 2008** (p. 2) gives the identical operator and convention:

> the dimensionless operator **C ≡ (n × S) · I/I** encodes the angular momentum dependence
> of H_eff_P [9]. … iW(m′_N,m′_I,m_N,m_I) ≡ κ′ W_P ⟨ψ⁻_↓(m′_N,m′_I)|C|ψ⁺_↑(m_N,m_I)⟩.
> **Time-reversal invariance ensures that iW is pure imaginary.**

DeMille's selection-rule statement (p. 2), which is load-bearing for §3–§4:

> **C is a pseudoscalar, with non-zero m.e.'s C̃ between states with the same value of
> m_F ≡ m_N + m_S + m_I.**

### 1.1 Karthein ↔ DeMille mapping

The two papers are the same operator with renamed electronic constants:

| quantity | Karthein | DeMille | identical? |
|---|---|---|---|
| operator | `C = (n×S)·I/I` | `C ≡ (n×S)·I/I` | **yes, verbatim** |
| effective H | `H_eff = κ′ W_A C` | `H_eff_P = κ′ W_P C` | yes; `W_A ≡ W_P` (electronic m.e.) |
| PV m.e. | `iW = κ′ W_A ⟨−|C|+⟩` | `iW = κ′ W_P ⟨−|C|+⟩` | yes, verbatim |
| imaginary | Eq. B1: `iW` off-diagonal | "T-reversal ⇒ iW pure imaginary" | yes |
| `⟨C⟩` at crossing | `C = 0.5` **assumed** (Fig. 3 caption — single occurrence in the paper, alongside assumed κ′=0.05) | `C̃⁽ᵐ⁾` "at the level-crossing with the maximum value of m_F" (Table I caption, p. 4) | same operator; per-channel values in §4b |

**No factor differences.** Karthein's `W_A` and DeMille's `W_P` are the same electronic
reduced matrix element; only the subscript differs (`A` for anapole-inclusive, `P` for
parity). The angular operator `C` is identical character-for-character. The note therefore
derives the matrix element of the single operator `C` and the implementation multiplies by
`κ′ W_A` downstream (a scalar constant), exactly as `build_PTV_bBJ` multiplies the EDM
operator by its electronic coefficient.

**Sign/phase convention adopted (matches both papers + the codebase EDM path):** the
function `NSDPV_bBJ` returns the **real number** `Im⟨bra|C|ket⟩`; the *builder* carries the
explicit `i`:

```
⟨bra|C|ket⟩ = (i/I) · NSDPV_bBJ(bra, ket)            # pure imaginary; /I is the literal c-number (×2 for I=1/2)
H_PV_matrix[i,j] = κ′ W_A · (i/I) · NSDPV_bBJ(...)   # = iW, matching Karthein Eq. B1 / DeMille Eq. (1)
```

This mirrors `build_PTV_bBJ` (`hamiltonian_builders.py:302–313`), which builds the (real,
Hermitian) EDM operator as `H_PTV[i,j] = -EDM*Sz_bBJ(**q_args)`; the EDM operator `Σ·n̂` is
T,P-odd but **real** in this basis, whereas the NSD-PV `C` is T,P-odd and **imaginary** — so
the analogue here is `H_PV[i,j] = coeff * 1j * NSDPV_bBJ(**q_args) / I` (the `i` is the
operator's; the builder's overall sign is pinned to Karthein Eq. (B1), see §4d flag 2). See
§4a for why the matrix is imaginary-Hermitian and §4d for two implementation traps.

---

## 2. Spherical-tensor decomposition (each identity cited)

### 2.1 Cross product → compound rank-1 tensor

B&C couple two rank-1 spherical tensors with **Eq. (5.109)** (PDF p. 192):

> `T^K_P(R^{k1}, S^{k2}) = Σ_{p1,p2} ⟨k1, k2, p1, p2 | K, P⟩ T^{k1}_{p1}(R) T^{k2}_{p2}(S)`

and the spherical scalar product with **Eq. (5.111)** (PDF p. 193):

> `T^k(R) · T^k(S) = Σ_p (−1)^p T^k_p(R) T^k_{−p}(S)`.

B&C **p. 193** states the cross-product correspondence explicitly:

> "two vectors u and v can be coupled together to give a scalar product T⁰(u,v), a vector
> **T¹(u,v) which is proportional to the vector product u∧v**, and a second-rank spherical
> tensor T²(u,v)."

The proportionality constant (the one factor the whole derivation hinges on) follows from
B&C **Eqs. (5.113)–(5.118)** (PDF pp. 193–194), where `u_m`, `v_m` are the *spherical*
components (B&C Table 5.2 / Eq. 5.108): e.g.
`T¹_0(u,v) = (1/√2)(u₁v₋₁ − u₋₁v₁) = (1/√2)(u_X v_Y − u_Y v_X) = (i/√2)(u×v)_Z`
(transcription corrected per review: the CG evaluation of (5.114)'s middle expression carries
a **minus** between the two terms; the Cartesian endpoint is unchanged)
(the final `i` enters because the spherical components `u_{±1} = ∓(u_X ± i u_Y)/√2` carry it).
**Derived (≤3 steps) and numerically verified** (§4, 3 random vectors, machine precision):

> **(u × v) = −i√2 · T¹(T¹(u), T¹(v))**, i.e. **(u·v cross dotted) → (A×B)·C = −i√2 [T¹(T¹(A),T¹(B))]·T¹(C).**

Applied to C:

> **C = (n̂×S)·I / I = −i√2 · [ T¹( T¹(n̂), T¹(S) ) · T¹(I) ] / I** … (★)

where `T¹(n̂)` is the rank-1 direction-cosine tensor on the rotational space (`T¹_q(n̂) = C¹_q(θ,φ)`,
the normalized rank-1 spherical harmonic; B&C Eq. (5.120) class, here k=1), `T¹(S)` acts on
the electron spin, **`T¹(I)` is the rank-1 spherical tensor of the nuclear-spin operator Î**
(components `T¹_q(Î)`, reduced m.e. `√(I(I+1)(2I+1))`), and the trailing **`/I` is the
literal c-number division by the nuclear-spin quantum number** (I=1/2 ⇒ ×2).

**Normalization (RESOLVED — literal `/I` adopted; supersedes this note's earlier "no-/I"
reading).** Provenance, DeMille PRL 100, 023003 (verified at the local PDF):

- p. 2 defines `C ≡ (n×S)·I/I` with the explicit `/I`, and Karthein's operator is verbatim
  DeMille's (§1.1), so the same `/I` applies.
- **Table I (p. 4) is decisive.** The tabulated `C̃⁽ᵐ⁾` values across eight nuclei spanning
  **I = 9/2 → 1/2** are **flat in I**: ⁸⁷Sr −0.40, ⁹¹Zr −0.40, ¹³⁷Ba −0.44, ¹⁷¹Yb −0.52,
  ²⁷Al −0.42, ⁶⁹Ga −0.43, ⁸¹Br −0.42, ¹³⁹La −0.43 (caption: "Superscript (m) indicates the
  value at the level-crossing with the maximum value of m_F"). Only the literal-`/I` operator
  is I-independent at the max-m_F crossing: there the q=0 nuclear factor is `m_I = I`, which
  cancels the `/I` exactly, leaving the I-independent base `1/√6 ≈ 0.408` (verified, §4b;
  hyperfine/spin-rotation admixtures move it to the −0.40…−0.52 spread). **Without** the `/I`
  the column would scale ∝ I (e.g. ⁸⁷SrF, I=9/2: ≈ 9/2 × 0.41 ≈ **1.84**, not −0.40).

**Why the earlier "no-/I" resolution was wrong (post-mortem).** It anchored on Karthein's
"C = 0.5", treating it as a computed F=0-crossing matrix element. In fact "C = 0.5" occurs
exactly **once** in Karthein, **assumed** in the Fig. 3 caption alongside the assumed
κ′ = 0.05 — it is a representative round number, not a computed gate. And the zero-field
F=0 singlet used in that check does not exist at the 1.5 T operating point: the crossing
eigenstates there are **decoupled products** `|N,m_N⟩|S,m_S⟩|I,m_I⟩` (Karthein's own basis
statement, p. 2) — the S-Zeeman destroys the F-coupling. The no-/I/coupled-F=0 match to 0.5
was a coincidence of two wrong anchors. Consistency check with the literal `/I`:
κ′W_A·|⟨C⟩| = (0.05 × 16 Hz) × (0.408…0.577) = **0.33–0.46 Hz** over Karthein's decoupled
crossing channels, bracketing his quoted W/2π = 0.4 Hz within rounding.

The **−i√2** is the standard cross-product compound-tensor factor (e.g. Zare, *Angular
Momentum*, and matches DeMille's "iW pure imaginary"): it makes C intrinsically imaginary,
which is the spherical-tensor expression of T-reversal oddness.

### 2.2 Three-space decoupling for case bβJ

`(★)` is a **scalar product of two rank-1 tensors**: `T¹(n̂,S)` (acting on the rotational ⊗
electron-spin space, coupled as `J = N + S`) dotted with `T¹(I)` (nuclear space). This is
structurally **identical to B&C Eq. (8.220)** (the Fermi-contact `bF T¹(S)·T¹(I)` element,
PDF p. 472 — the codebase's audited `IS_bBJ`), except `T¹(S)` (pure spin, N-diagonal) is
replaced by the **compound** `T¹(n̂,S)` (built from n̂ on N and S on spin → ΔN=±1).

**Step A — F-level decoupling (I from J).** Exactly the B&C **Eq. (8.220)** F-reduction
(PDF p. 472): for a scalar product `T¹(A_{NS}) · T¹(I)` with A acting on the J-space,

> `⟨N S J I F M_F| T¹(A_{NS})·T¹(I) |N′ S J′ I F′ M′_F⟩`
> `= δ_{FF′} δ_{M_F M′_F} (−1)^{J′+I+F} {J I F; I J′ 1} ⟨N S J‖T¹(A_{NS})‖N′ S J′⟩ ⟨I‖T¹(I)‖I⟩`

with `⟨I‖T¹(I)‖I⟩ = √(I(I+1)(2I+1))` (B&C Eq. 8.220, p. 472). The 6j `{J I F; I J′ 1}` is
identical (verified, 0 mismatches over all J,J′,F) to the codebase `IS_bBJ` 6j
`{J1 I F0; I J0 1}` by the column-swap 6j symmetry the audit record already notes.

**Step B — J-level decoupling (S from N) of the compound `T¹(n̂,S)`.** The reduced matrix
element of a rank-1 compound tensor built from `T¹(n̂)` (on N) and `T¹(S)` (on S), between
J=N+S coupled states, is a **9-j**, now fully analytically anchored (review-confirmed):

> `⟨N S J‖T¹(C¹(n̂),S)‖N′ S J′⟩`
> `= √(3·(2J+1)(2J′+1)) · {N S J; N′ S J′; 1 1 1}(9j) · ⟨N‖C¹(n̂)‖N′⟩ · ⟨S‖T¹(S)‖S⟩`

per **B&C Eq. (5.139)** (PDF p. 198) — the compound-tensor reduced m.e., which carries
**no phase factor** (the `√3` is `√(2k+1)` for the compound rank k=1). The direction-cosine
reduced m.e. carries the phase of **B&C Eq. (5.148)** (PDF p. 199):

> `⟨N K‖C¹(n̂)‖N′ K′⟩ = (−1)^{N−K} √((2N+1)(2N′+1)) (N 1 N′; −K 0 K′)`   (K=K′=0 here)

and the scalar-product phase `(−1)^{J′+I+F}` of Step A is **B&C Eq. (5.140)** (PDF p. 198,
the equation underlying (8.220)). `⟨S‖T¹(S)‖S⟩ = √(S(S+1)(2S+1))`. With these three anchors
the previously-empirical `(−1)^{N′}·√2` is **derived**, not pinned: for K=0 and ΔN=±1,
`(−1)^{N−K} = −(−1)^{N′}`, and the global `−√2` is the real factor of the cross-product
`−i√2` (§2.1; the `i` goes to the builder). The empirically-validated form of the first
draft and this canonical form are algebraically identical — **verified exactly (difference
0.0) over every state pair in the N≤2 basis**. UNRESOLVED-1 is closed (§5).

The `(N 1 N′; 0 0 0)` 3-j vanishes unless **N + 1 + N′ is even**, i.e. **ΔN = ±1** (and
ΔN=0 forbidden) — this is the parity-odd, ΔN=±1 selection rule that makes H_PV connect
opposite-parity (`P = (−1)^N`) rotational levels.

---

## 3. Final formula (canonical, K-general, codebase notation)

Combining Step A and Step B with `(★)`: the function returns the real factor; the operator
`C` is `(i/I) ×` the function (the `i` and the `/I` are carried by the builder). Index 1 =
ket = primed. This is the **reviewer's canonical K-general form** — algebraically identical
to the first draft's empirically-pinned form (verified exactly, §4a), and letter-for-letter
in the codebase's direction-cosine convention (cf. `Sz_bBJ`, `StarkZ_bBJ`, `T2IS_bBJ`:
`(-1)**(N0-K0) * wigner_3j(N0,1,N1,-K0,0,K1)`).

```python
def NSDPV_bBJ(K0,N0,J0,F0,M0, K1,N1,J1,F1,M1, S=1/2, I=1/2):
    # Im< (n_hat x S) . I_operator >. The NSD-PV operator C = (n_hat x S).I/I has
    # matrix element  <bra|C|ket> = (i/I) * NSDPV_bBJ  (builder carries the i and the /I).
    # Rank-0 pseudoscalar: dK=0, dF=0, dM=0; dN=+-1 via the 3j (parity-odd channel).
    if not (kronecker(K0,K1)*kronecker(F0,F1)*kronecker(M0,M1)):
        return 0
    return -np.sqrt(2) \
        * (-1)**(J1 + I + F0) * wigner_6j(J0, I, F0, I, J1, 1) \
        * np.sqrt(I*(I+1)*(2*I+1)) \
        * np.sqrt(3*(2*J0+1)*(2*J1+1)) * wigner_9j(N0, S, J0, N1, S, J1, 1, 1, 1) \
        * (-1)**(N0-K0) * np.sqrt((2*N0+1)*(2*N1+1)) * wigner_3j(N0, 1, N1, -K0, 0, K1) \
        * np.sqrt(S*(S+1)*(2*S+1))
```

**Per-factor B&C provenance** (all PDF pages; review-confirmed):
- `(−1)^{J1+I+F0} · 6j · √(I(I+1)(2I+1))` — scalar-product F-reduction phase, **Eq. (5.140)**
  (PDF p. 198); realized for this coupling order as **Eq. (8.220)** (PDF p. 472, = `IS_bBJ`).
- `√(3(2J0+1)(2J1+1)) · 9j` — compound-tensor reduced m.e., **Eq. (5.139)** (PDF p. 198);
  **no phase factor** in this equation.
- `(−1)^{N0−K0} · √((2N0+1)(2N1+1)) · 3j(N0,1,N1,−K0,0,K1)` — direction-cosine reduced m.e.
  phase, **Eq. (5.148)** (PDF p. 199).
- global `−√2` — the real factor of the cross-product `−i√2` (§2.1, from Eqs. (5.113)–(5.118)).

Equivalently in closed form (bra unprimed, ket primed):

> **I · Im⟨N J F M| C |N′ J′ F′ M′⟩ = δ_{KK′} δ_{FF′} δ_{MM′} × (−√2)**
> **× (−1)^{J′+I+F} {J I F; I J′ 1} √(I(I+1)(2I+1))**
> **× √(3(2J+1)(2J′+1)) · {N S J; N′ S J′; 1 1 1} · (−1)^{N−K} √((2N+1)(2N′+1)) (N 1 N′; −K 0 K′) · √(S(S+1)(2S+1))**

and `⟨N J F M| C |N′ J′ F′ M′⟩ = (i/I) ·` NSDPV_bBJ.

**Antisymmetry:** `NSDPV_bBJ(a,b) = −NSDPV_bBJ(b,a)` (verified exact over the full N≤2
basis), so `(i/I)·NSDPV` is Hermitian, matching the `iW`/`−iW` pattern of Karthein Eq. (B1).
**Verification depth:** identical to the first-principles decoupled construction over the
full N≤2 basis (36 states, **46 nonzero matrix elements**, max|Δ| = 2.2×10⁻¹⁶ — re-run this
session; reviewer's independent run 1.7×10⁻¹⁶).

**Selection rules the formula enforces** (all verified, §4):
- **ΔK = 0, ΔF = 0, ΔM = 0** — `C` is a total rank-0 (pseudo)scalar (kronecker guards +
  the F-level 6j with rank 1 contracted to F-scalar). M-independent within an F-multiplet.
- **ΔN = ±1** — for K=0, from `(N 1 N′; 0 0 0)`; ΔN even (incl. 0) forbidden. This is the
  parity-odd channel (P = (−1)^N), the defining feature of H_PV.
- ΔJ = 0, ±1 — allowed by the 9j (J connects via rank-1); the F-conservation with ΔJ≠0 is
  carried by the `{J I F; I J′ 1}` 6j.

There is **no separate `√((2F+1)(2F′+1))` factor** (unlike the Zeeman/Stark elements): the
8.220-style scalar-product F-reduction puts the full F-dependence in the single 6j and the
δ_{FF′}, exactly as in the audited `IS_bBJ`. (`IS_bBJ` likewise has no `(2F+1)` factor.)

---

## 4. Expected properties (the gates the implementation must pass)

All gates below were checked by building `C = −i√2 [T¹(n̂,S)·T¹(I)]/I` in the fully decoupled
basis `|N,m_N⟩|S,m_S⟩|I,m_I⟩` from first principles (Gaunt/3j direction cosine + explicit
spin-1/2 spherical tensors + Clebsch–Gordan compound coupling), transforming to bβJ with the
codebase's own audited `decouple_b_even` (`matrix_elements.py`), and comparing against the
closed form of §3.

### (a) Anti-Hermitian / imaginary structure

- In the decoupled basis, **C is purely imaginary** (max|Re| = 0; max|Im| = O(1)) —
  **verified**. This is the spherical-tensor expression of the `−i√2` factor and of
  T-reversal oddness (DeMille: "iW pure imaginary").
- In the bβJ basis, the matrix `C` is **imaginary and Hermitian** (`C = C†`, max|C−C†| =
  2.8×10⁻¹⁷ — **verified**). An imaginary Hermitian matrix is real-antisymmetric × i:
  `C_{ba} = −C_{ab}` for the real coefficients, which is exactly the `iW` (upper) / `−iW`
  (lower) off-diagonal pattern of Karthein **Eq. (B1)**. So the implementation should
  build `H_PV[i,j] = κ′ W_A · i · NSDPV_bBJ(i,j)` with `NSDPV_bBJ` returning the real
  antisymmetric coefficient; the resulting matrix is anti-Hermitian-free (it IS Hermitian,
  because the `i` is part of the operator, not an ad-hoc factor).
- **How the codebase handles the analogue:** `build_PTV_bBJ` (`hamiltonian_builders.py:302`)
  builds `H_PTV[i,j] = -EDM*Sz_bBJ(**q_args)` — a **real** matrix, because the EDM operator
  `Σ·n̂` (registered via `PTV_builders` in `molecule_library_class.py:386–403`, the
  `build_PTV_bBJ`/`build_PTV_bBS` entries) is T,P-odd but **real** in this basis. The NSD-PV
  `C` differs: it is T,P-odd and **imaginary**, so its builder must carry an explicit `1j`
  that the EDM builder does not (see §4d flag 1 for the dtype trap this implies). For **real**
  eigenvectors `v` (the field-on Hamiltonian here is real-symmetric), `v·(iA)·v ≡ 0`
  identically for real antisymmetric `A` — so the diagonal PV expectation `evec@H_PV@evec`
  (Energy_Levels.py:486) is **identically zero, as P-oddness demands** (wording corrected per
  review; the earlier draft said "real"). The PV observable is the off-diagonal mixing `iW`
  between near-degenerate opposite-parity states, never a diagonal shift.
- **Full-matrix gate:** `(i/I)·NSDPV_bBJ` reproduces the first-principles `C` over the entire
  **N≤2 bβJ basis (36 states, 46 nonzero matrix elements) to max|ΔC| = 2.2×10⁻¹⁶** —
  **verified** (canonical form; also exactly equal to the first draft's pinned form, Δ = 0.0).

### (b) Gate values (literal `/I`; all verified this session, machine precision)

Karthein Fig. 3 caption (p. 4): W/2π = 0.4 Hz "when assuming **κ′ = 0.05 and C = 0.5**" —
both are **assumed representative values** ("C = 0.5" occurs exactly once in the paper, in
that caption; it is not a computed matrix element and is **not a gate**). Consistency: with
the literal `/I`, κ′W_A·|⟨C⟩| = (0.05 × 16 Hz) × (0.408…0.577) = **0.33–0.46 Hz** over the
decoupled crossing channels — Karthein's 0.4 Hz sits inside, within rounding.

**The operating-point basis.** At the 1.5167 T crossings the eigenstates are **decoupled
products** `|N,m_N⟩|S,m_S⟩|I,m_I⟩` (Karthein's basis statement, p. 2): **both S and I are
decoupled** there (the earlier draft's "S decoupled but I not" mechanism was wrong, as was
its no-/I normalization — see §2.1 post-mortem). The physically operative gates are the
decoupled-channel matrix elements:

**PRIMARY implementation gate — the canonical flip-flop crossing channel:**

> **⟨0, 0, +½, −½| C |1, 0, −½, +½⟩ = +i/√3**, i.e. **|⟨C⟩| = 1/√3 ≈ 0.57735** (**verified**)

(labels `|N, m_N, m_S, m_I⟩`; this is the task's named pair A/B — the ΔN=1, m_N-preserving,
S–I flip-flop channel, m_F = 0 on both sides.)

**Stretched / max-m_F channel — DeMille's Table I base value:**

> **⟨0, 0, +½, −½| C |1, +1, −½, −½⟩ = −i/√6**, i.e. **|⟨C⟩| = 1/√6 ≈ 0.40825** (**verified**)

This is the I-independent base of DeMille's `C̃⁽ᵐ⁾` column (−0.40…−0.52 across I = 9/2…1/2):
at max m_F the nuclear q=0 factor `m_I = I` cancels the `/I` exactly — the Table-I flatness
that pins the literal-`/I` normalization (§2.1).

**Zero-field coupled-basis cross-checks** (bβJ eigenstates of the B=0 Hamiltonian; useful as
formula gates even though they are not the 1.5 T operating basis; all **verified**):

| bra (N=0) | ket (N=1) | ⟨C⟩ | magnitude |
|---|---|---|---|
| J=½, F=0 | J=½, F=0 | **−1.0 i** | 1.0 |
| J=½, F=1 | J=½, F=1 | **+i/3** (all M) | 0.33333 |
| J=½, F=1 | J=1½, F=1 | **−i√2/3** | 0.47140 |

#### Addendum (2026-06-11, Task-4 implementation): field-dressed |⟨C⟩| at B_c — the "~0.5%" decoupling-correction estimate superseded

Measured with the implemented operator (`test_nsdpv_operator.py`, gate G3) on eigenvectors of
the full Hamiltonian at the canonical (0,0) flip-flop crossing, B_c = 15167.2 G:

> **dressed |⟨C⟩| = 0.59015** — **+2.22% above the decoupled 1/√3 = 0.57735**, and
> **offset-independent** (identical for off-crossing evaluation offsets 0.5–50 G, and equal
> to the at-crossing tracked-evec value). This is the true, stable field-dressed matrix
> element, not a numerical artifact.

**Mechanism.** Partner A is 99.991% pure `|0,0,+½,−½⟩`; partner B is 99.887%
`|1,0,−½,+½⟩` **plus a 0.105%-weight (≈3.2% amplitude) admixture of the stretched ket
`|1,+1,−½,−½⟩`**, mixed by the dipolar `c` coupling across the **hyperfine-scale** gap
inside the N=1, m_F=0, m_S=−½ manifold (the model has no rotational Zeeman term, so there is
no Zeeman-scale suppression of this intra-manifold mixing). The dressed element is the
coherent sum of the flip-flop and stretched channels, and the arithmetic closes:
1/√3 + 0.032×(1/√6) ≈ 0.577 + 0.013 = 0.590.

**Correction recorded:** this section's earlier "~0.5%" decoupling-correction estimate
implicitly assumed Zeeman-scale suppression of the channel mixing and is **wrong by ~4×**;
the actual correction at this crossing is +2.22%. Gate G3's magnitude tolerance is
recalibrated 2% → 3% accordingly (sign pin unchanged and strict; the pure-decoupled gate
values above remain exact at 1e-9). Consequence for the headline number:
**W/2π = κ′W_A·|⟨C⟩| = 0.05 × 16 Hz × 0.59015 = 0.4721 Hz** at the canonical crossing
(vs 0.462 Hz from the naive 1/√3), nudging the §2.1 consistency bracket's upper edge from
0.46 to 0.47 Hz. Verified independently by Arian (re-ran gate, same numbers; decision A,
2026-06-11).

**Full PV-active crossing map** (gate G5, field-dressed, ΔM_F = 0 only; the five
M_F(A) ≠ M_F(B) crossings of the Task-3 table at 15151.2 / 15155.9 / 15158.0 / 15160.4 /
15162.4 G all have ⟨C⟩ = 0 to machine precision — the ΔM_F kill of §(c)):

| B_c [G] | channel (partner B character) | dressed \|⟨C⟩\| | sign of ⟨+\|C\|−⟩/i | W/2π [Hz] (κ′=0.05, W_A/2π=16 Hz) |
|---|---|---|---|---|
| 15017.0 | (0,0)′ stretched-partner — B = `\|1,+1,−½,−½⟩` | 0.39282 | − | 0.3143 |
| 15167.2 | (0,0) canonical flip-flop — B = `\|1,0,−½,+½⟩` | 0.59015 | **+** | 0.4721 |
| 15296.9 | (+1,+1) stretched | 0.40411 | − | 0.3233 |

(The 15017.0 G crossing is the **second (0,0)** crossing — same partner A, the
hyperfine-shifted stretched partner of the same manifold; stretched-partner character
confirmed by decoupled decomposition in G5.)

### (c) Sign pattern / selection rules across the 7 crossings

**The decisive fact: C conserves m_F (ΔM_F = 0).** DeMille (p. 2): "C is a pseudoscalar,
with non-zero m.e.'s between states with the same value of m_F." **Verified**: every nonzero
decoupled and bβJ matrix element has m_F(bra) = m_F(ket); the element is M-independent within
an F-multiplet (F=1→F=1 gives +i/3 ≈ +0.33333 i for **all** M = −1, 0, +1 — verified).

The task lists six `(M_a, M_b)` crossing combinations:
`(+1,−1), (0,+1), (+1,0), (+1,−2), (0,−1), (+1,+1)`. Applying ΔM_F = 0:

| (M_a, M_b) | M_a = M_b ? | ⟨C⟩ |
|---|---|---|
| (+1, −1) | no | **0 (vanishes, ΔM_F = 2)** |
| (0, +1)  | no | **0 (vanishes, ΔM_F = 1)** |
| (+1, 0)  | no | **0 (vanishes, ΔM_F = 1)** |
| (+1, −2) | no | **0 (vanishes, ΔM_F = 3)** |
| (0, −1)  | no | **0 (vanishes, ΔM_F = 1)** |
| (+1, +1) | **yes** | can be ≠ 0 |

**Only the (+1,+1) combination survives**; the other five are exactly zero by the
M_F-conservation selection rule. **This kills the "crossing-sign-reversal handle" framed as
coming from (M_a, M_b) with M_a ≠ M_b** — those crossings have ⟨C⟩ = 0 and produce no PV
signal at all (consistent with the experiment selecting m_F-matched opposite-parity pairs).

The genuine **sign reversal across crossings** (DeMille p. 3: "the ratio C̃/d varies widely
but deterministically in magnitude and sign between nearby level crossings") comes instead
from the **intrinsic sign of ⟨C⟩ for the M_F-conserving crossings**, which alternates by
channel. The verified bβJ ⟨C⟩ values (the only nonzero ones, all with ΔM_F=0, ΔN=1):

| bra (N=0) | ket (N=1) | ⟨C⟩ (literal /I) |
|---|---|---|
| J=½, F=0 | J=½, F=0 | **−1.0 i** |
| J=½, F=1 | J=½, F=1 | **+i/3** (≈ +0.33333 i) |
| J=½, F=1 | J=1½, F=1 | **−i√2/3** (≈ −0.47140 i) |

(and the Hermitian conjugates N=1→N=0 with the sign flipped per `C_{ba} = −C_{ab}`.) The
sign **alternates** between the F=0 (−), F=1/J=½ (+), and F=1/J=1½ (−) channels — this is the
deterministic per-crossing sign reversal DeMille relies on for systematics rejection, and it
is **entirely within ΔM_F=0** crossings. The implementation must reproduce these three signed
magnitudes.

### (d) Implementation flags (added per review — traps for the Task-4 coding step)

1. **Complex dtype allocation.** `build_PTV_bBJ` allocates its matrix as `np.zeros((n,n))` —
   **real** dtype. Copying that pattern for NSD-PV silently drops the physics: assigning
   `1j * NSDPV` into a real array discards the imaginary part (numpy casts on assignment).
   The NSD-PV builder must allocate `np.zeros((n,n), dtype=complex)` **before** multiplying
   by `1j`.

2. **W-sign convention — an internal tension in both papers; pin to Eq. (B1).** Karthein's
   p. 2 definition `iW ≡ κ′W_A⟨Ψ⁻_↓|C|Ψ⁺_↑⟩` says `⟨−|H_PV|+⟩ = +iW`; but Eq. (B1) in the
   (+,−) basis ordering places `−iW` in that slot (i.e. B1 ⟹ `⟨+|H_PV|−⟩ = +iW`,
   `⟨−|H_PV|+⟩ = −iW`). The two statements differ by an overall sign. The **identical**
   tension exists in DeMille (Eq. (1) layout vs his iW definition — verified at the PDF,
   p. 2), so this is an inherited notational looseness, harmless for magnitudes and for the
   asymmetry observable. The repo's conventions record (`docs/sio-conventions.md` +
   `Source Code/test_sio_conventions.py`) anchors SiO⁺ conventions on Karthein's explicit
   equations, so the implementation **pins to Eq. (B1)**: a **positive real coefficient** for
   (bra = |+⟩ = even-N row, ket = |−⟩ = odd-N column) gives `⟨+|H_PV|−⟩ = +iW`.

---

## 5. Review items (status after the 2026-06-11 spectroscopy-reviewer audit)

1. **RESOLVED — the empirical `(−1)^{N′}` phase IS the B&C-derived one.** The reviewer
   anchored every factor analytically: **Eq. (5.139)** (PDF p. 198) gives the compound
   rank-1 9j reduced m.e. with **no phase**; **Eq. (5.140)** (PDF p. 198) gives the
   scalar-product phase `(−1)^{J′+I+F}`; **Eq. (5.148)** (PDF p. 199) gives the
   direction-cosine phase `(−1)^{N−K}`. For K=0, ΔN=±1, `(−1)^{N−K} = −(−1)^{N′}`, so the
   first draft's `+√2·(−1)^{N′}` equals the canonical `−√2·(−1)^{N0−K0}` exactly. §3 now
   states the formula in the canonical `(−1)^{N0−K0}·wigner_3j(N0,1,N1,−K0,0,K1)` form,
   letter-for-letter the codebase's direction-cosine convention (`Sz_bBJ`, `StarkZ_bBJ`,
   `T2IS_bBJ`). Equality of the two forms verified exactly (Δ = 0.0 over every N≤2 pair);
   first-principles agreement 2.2×10⁻¹⁶ (46 nonzero MEs).

2. **LAYOUT NOTE (stays open; non-blocking) — Karthein App. A/B layout vs the task
   description.** The task refers to "Karthein Eq. (B1)
   and Appendix B define … the C operator; Appendix A the basis." In the local copy
   (arXiv:2310.11192v1, Zotero LYF85P8C): **App. A** = "Effects of Time-Varying Electric
   Fields" (gives the decoupled basis `|N,mN⟩|S,mS⟩|I,mI⟩`, H_eff Eq. (A2), n̂ in spherical
   harmonics Eqs. (A5)–(A6)); **App. B** = "Asymmetry Analytical Formula" (gives **Eq. (B1)**,
   the 2×2 `H±` with `iW`). The **C operator definition and `H_eff = κ′ W_A C` are in the
   main text on p. 2**, not in an appendix, in this version. The published-PRL labels may
   differ; cited locations above are this local copy. (Non-blocking; the operator and Eq. (B1)
   are both quoted verbatim above.)

3. **RESOLVED — literal `/I` adopted (supersedes both of this note's earlier readings).**
   The first draft kept `/I` symbolic; an interim revision wrongly concluded "no `/I`" by
   anchoring on Karthein's `C = 0.5` — which turned out to be an **assumed** Fig.-3-caption
   value evaluated against a zero-field F=0 singlet that does not survive at the 1.5 T
   operating point (§2.1 post-mortem). The decisive evidence is **DeMille Table I flatness**
   (verified at the PDF): `C̃⁽ᵐ⁾` is I-independent across I = 9/2…1/2, which only the
   literal-`/I` operator produces (the max-m_F channel has nuclear factor `m_I = I`, exactly
   cancelling `/I`; base 1/√6 ≈ 0.408 — verified). Without `/I`, SrF (I=9/2) would read
   ≈1.84, not −0.40. So: `⟨C⟩ = (i/I)·NSDPV_bBJ`, i.e. **×2 for I=1/2**, and the
   generalization to I≠1/2 is now unambiguous — divide by the quantum number I. The
   `√(I(I+1)(2I+1))` reduced-m.e. placement inside `NSDPV_bBJ` matches B&C Eq. (8.220)'s
   `{I(I+1)(2I+1)}^{1/2}` and the audited `IS_bBJ`.

---

## Validation summary (executed in `Structure` conda env; amended values use literal `/I`)

- `(A×B)·C = −i√2 [T¹(T¹A,T¹B)]·T¹C` factor: verified to machine precision over 3 random
  vector triples (`(n̂×S)·I` real = `−i√2 [T¹(n̂,S)·T¹(I)]`).
- Decoupled `C` purely imaginary: max|Re| = 0.
- bβJ `C` Hermitian: max|C−C†| = 2.8×10⁻¹⁷.
- Basis change `decouple_b_even` unitary: max|UᵀU−I| = 4.4×10⁻¹⁶.
- **Canonical §3 form ≡ first draft's empirically-pinned form: Δ = 0.0 (exact) over every
  state pair of the N≤2 basis** (36 states, K=0).
- `(i/I)·NSDPV_bBJ` vs first-principles `C`: **max|Δ| = 2.2×10⁻¹⁶ over the full N≤2 basis,
  46 nonzero MEs** (reviewer's independent run: 1.7×10⁻¹⁶).
- Antisymmetry `NSDPV(a,b) = −NSDPV(b,a)`: max|sum| = 0.0 (exact), N≤2.
- F-coupling 6j ≡ audited `IS_bBJ` 6j: 0 mismatches over all (J,J′,F).
- **PRIMARY gate** ⟨0,0,+½,−½|C|1,0,−½,+½⟩ = +i/√3 (|⟨C⟩| = 0.57735): verified.
- Stretched gate ⟨0,0,+½,−½|C|1,+1,−½,−½⟩ = −i/√6 (0.40825, DeMille Table-I base): verified.
- Zero-field coupled checks: F=0↔F=0 → −1.0 i; F=1(J½)↔F=1(J½) → +i/3 (all M);
  F=1(J½)↔F=1(J1½) → −i√2/3: all verified.
- DeMille Table I read at the PDF (p. 4): C̃⁽ᵐ⁾ flat (−0.40…−0.52) across I = 9/2…1/2;
  caption "at the level-crossing with the maximum value of m_F": confirmed verbatim.
- 5 of 6 task (M_a,M_b) crossings vanish by ΔM_F=0; only (+1,+1) survives: verified.
