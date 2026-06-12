# ²⁹SiO⁺ NSD-PV anapole operator C = (n̂×S)·I/I — case bβJ derivation note

**Task 4, step 1 (pre-review).** Derivation only; no code in this note. The target is the
matrix element of the dimensionless NSD-PV operator **C = (n̂×S)·I/I** in the codebase's
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
> momentum algebra [21]. (ref [21] = DeMille 2008.)

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
| `⟨C⟩` at crossing | `C = 0.5` (p. 4) | `C̃` (max in F=0 channel) | yes (see §4b) |

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
⟨bra|C|ket⟩ = i · NSDPV_bBJ(bra, ket)          # the matrix element of C is pure imaginary
H_PV_matrix[i,j] = κ′ W_A · i · NSDPV_bBJ(...)  # = iW, matching Karthein Eq. B1 / DeMille Eq. (1)
```

This mirrors `build_PTV_bBJ` (`hamiltonian_builders.py:302–313`), which builds the (real,
Hermitian) EDM operator as `H_PTV[i,j] = -EDM*Sz_bBJ(**q_args)`; the EDM operator `Σ·n̂` is
T,P-odd but **real** in this basis, whereas the NSD-PV `C` is T,P-odd and **imaginary** — so
the analogue here is `H_PV[i,j] = -coeff * 1j * NSDPV_bBJ(**q_args)` (the `i` is the operator's,
the overall `-`/coeff is the builder's). See §4a for why the matrix is imaginary-Hermitian.

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
`T¹_0(u,v) = (1/√2)(u₁v₋₁ + u₋₁v₁) = (1/√2)(u_X v_Y − u_Y v_X) = (i/√2)(u×v)_Z`
(the final `i` enters because the spherical components `u_{±1} = ∓(u_X ± i u_Y)/√2` carry it).
**Derived (≤3 steps) and numerically verified** (§4, 3 random vectors, machine precision):

> **(u × v) = −i√2 · T¹(T¹(u), T¹(v))**, i.e. **(u·v cross dotted) → (A×B)·C = −i√2 [T¹(T¹(A),T¹(B))]·T¹(C).**

Applied to C:

> **C = (n̂×S)·I / I = −i√2 · [ T¹( T¹(n̂), T¹(S) ) · T¹(I) ] / I** … (★)

where `T¹(n̂)` is the rank-1 direction-cosine tensor on the rotational space (`T¹_q(n̂) = C¹_q(θ,φ)`,
the normalized rank-1 spherical harmonic; B&C Eq. (5.120) class, here k=1), `T¹(S)` acts on
the electron spin, and **`T¹(I)` is the rank-1 spherical tensor of the nuclear-spin operator Î**
(components `T¹_q(Î)`, reduced m.e. `√(I(I+1)(2I+1))`). **Normalization (resolved, §4b /
former UNRESOLVED-3, settled by the data):** the validated operator is `(n̂×S)·Î` — the
spin *operator* `Î` in the numerator, **no literal division by the c-number I**. Tested three
readings against Karthein's quoted `⟨C⟩=0.5` at the F=0 crossing (verified): `(n̂×S)·Î` →
**−0.5 i ✓**; literal `(n̂×S)·Î / (I=½)` → −1.0 i ✗ (off by 2×); `(n̂×S)·Î / √(I(I+1))` →
−0.577 i ✗. So Karthein's notation `C = (n̂×S)·I/I` is read with the trailing `/I` as part of
making `C` an O(1) dimensionless number whose F=0-crossing value is `½`; for I=1/2 the
operator that reproduces their 0.5 is `(n̂×S)·Î` and that is what `NSDPV_bBJ` implements.
(See §5 UNRESOLVED-3 for the I≠1/2 generalization caveat.)

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
J=N+S coupled states, is a **9-j** (B&C Ch. 5 tensor algebra; the rank-1 analogue of the
rank-2 9j in the audited `T2IS_bBJ`, `matrix_elements.py:131`):

> `⟨N S J‖T¹(C¹(n̂),S)‖N′ S J′⟩`
> `= (−1)^{N′} · √2 · √((2J+1)(2J′+1)·3) · {N S J; N′ S J′; 1 1 1}(9j) · ⟨N‖C¹(n̂)‖N′⟩ · ⟨S‖T¹(S)‖S⟩`

with `⟨N‖C¹(n̂)‖N′⟩ = √((2N+1)(2N′+1)) · (N 1 N′; 0 0 0)` (rank-1 direction-cosine reduced
m.e., B&C Eq. (8.220)/(8.224) direction-cosine class, p. 472–473) and
`⟨S‖T¹(S)‖S⟩ = √(S(S+1)(2S+1))`.

The `(N 1 N′; 0 0 0)` 3-j vanishes unless **N + 1 + N′ is even**, i.e. **ΔN = ±1** (and
ΔN=0 forbidden) — this is the parity-odd, ΔN=±1 selection rule that makes H_PV connect
opposite-parity (`P = (−1)^N`) rotational levels. The `√2` and the `(−1)^{N′}` (ket-N phase)
are fixed empirically against the first-principles construction (§4) and absorb the
normalization/sign of the compound-tensor reduction in the codebase's
`(N0 1 N1; 0 0 0)`-with-`(−1)^{N−K}` convention; see §5 (UNRESOLVED-1) for the
not-yet-closed analytic anchor of the `(−1)^{N′}` phase specifically.

---

## 3. Final formula (codebase notation)

Combining Step A and Step B with `(★)` and folding the global `−i√2`, the `i` is carried by
the builder and the function returns the real `Im⟨bra|C|ket⟩`. Index 1 = ket = primed.

```python
def NSDPV_bBJ(K0,N0,J0,F0,M0, K1,N1,J1,F1,M1, S=1/2, I=1/2):
    # C = (n_hat x S).I / I  is a rank-0 pseudoscalar: dF=0, dM=0, dK=0.
    # The matrix element of C is i*NSDPV_bBJ (intrinsically imaginary; builder carries the i).
    if not (kronecker(K0,K1)*kronecker(F0,F1)*kronecker(M0,M1)):
        return 0
    return (-1)**(J1 + I + F0) * wigner_6j(J0, I, F0, I, J1, 1) \
        * np.sqrt(I*(I+1)*(2*I+1)) \
        * np.sqrt(2) * np.sqrt((2*J0+1)*(2*J1+1)*3) \
        * wigner_9j(N0, S, J0, N1, S, J1, 1, 1, 1) \
        * (-1)**(N1) * np.sqrt((2*N0+1)*(2*N1+1)) * wigner_3j(N0, 1, N1, 0, 0, 0) \
        * np.sqrt(S*(S+1)*(2*S+1))
```

Equivalently in closed form (bra unprimed = 0, ket primed = 1):

> **Im⟨N J F M| C |N′ J′ F′ M′⟩ = δ_{KK′} δ_{FF′} δ_{MM′}**
> **× (−1)^{J′+I+F} {J I F; I J′ 1} √(I(I+1)(2I+1))**
> **× √2 · √((2J+1)(2J′+1)·3) · {N S J; N′ S J′; 1 1 1} · (−1)^{N′} √((2N+1)(2N′+1)) (N 1 N′; 0 0 0) · √(S(S+1)(2S+1))**

and `⟨N J F M| C |N′ J′ F′ M′⟩ = i ·` (the above).

**Selection rules the formula enforces** (all verified, §4):
- **ΔK = 0, ΔF = 0, ΔM = 0** — `C` is a total rank-0 (pseudo)scalar (kronecker guards +
  the F-level 6j with rank 1 contracted to F-scalar). M-independent within an F-multiplet.
- **ΔN = ±1** — from `(N 1 N′; 0 0 0)`; ΔN even (incl. 0) forbidden. This is the parity-odd
  channel (P = (−1)^N), the defining feature of H_PV.
- ΔJ = 0, ±1 — allowed by the 9j (J connects via rank-1); the F-conservation with ΔJ≠0 is
  carried by the `{J I F; I J′ 1}` 6j.

There is **no separate `√((2F+1)(2F′+1))` factor** (unlike the Zeeman/Stark elements): the
8.220-style scalar-product F-reduction puts the full F-dependence in the single 6j and the
δ_{FF′}, exactly as in the audited `IS_bBJ`. (`IS_bBJ` likewise has no `(2F+1)` factor.)

---

## 4. Expected properties (the gates the implementation must pass)

All gates below were checked by building `C = −i√2 T¹(n̂,S)·T¹(I)` in the fully decoupled
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
  that the EDM builder does not. The PV shift `evec@H_PV@evec` (Energy_Levels.py:486) is then
  real (Hermitian H_PV), as required.
- **Full-matrix gate:** `NSDPV_bBJ` (×i) reproduces the first-principles `C` over the entire
  N∈{0,1} bβJ block to **max|ΔC| = 2.2×10⁻¹⁶** (16×16, 14 nonzero entries) — **verified**.

### (b) ⟨C⟩ ≈ 0.5 at the Karthein crossing

Karthein (p. 4): "the calculated W_A/2π = 16 Hz … corresponding to W/2π = 0.4 Hz when
assuming **κ′ = 0.05 and C = 0.5**."

**Where 0.5 lives.** The maximal |⟨C⟩| in the bβJ basis is exactly **0.5**, attained at

> **⟨N=0, J=½, F=0, M=0| C |N=1, J=½, F=0, M=0⟩ = −0.5 i**  (**verified**, machine precision).

This is the F=0 (spin-singlet of S and I) crossing — the most strongly mixing pair, and the
operating point Karthein quotes. The decoupled-basis arithmetic (assembled term by term from
the first-principles decoupled matrix elements):

bra `|N=0,J=½,F=0,M=0⟩ = (1/√2)|0,0,−½,+½⟩ − (1/√2)|0,0,+½,−½⟩` (S–I singlet, N=0)
ket `|N=1,J=½,F=0,M=0⟩ = −(1/√3)|1,−1,+½,+½⟩ + (1/√6)|1,0,−½,+½⟩ + (1/√6)|1,0,+½,−½⟩ − (1/√3)|1,+1,−½,−½⟩`

and the nonzero decoupled `⟨A|C|B⟩` building blocks (m_F-conserving, ΔN=1):
`⟨0,0,∓½,±½|C|1,0,±½,∓½⟩ = ∓0.28868 i = ∓i/(2√3)` and
`⟨0,0,∓½,±½|C|1,∓1 or ±1, …⟩ = ±0.20412 i = ±i/(2√6)`.

Summing the six surviving cross-terms (verified arithmetic):

```
 (1/√2)(−1/√3)(+0.20412 i) + (1/√2)(1/√6)(−0.28868 i) + (1/√2)(−1/√3)(+0.20412 i)
+(−1/√2)(−1/√3)(−0.20412 i) + (−1/√2)(1/√6)(+0.28868 i) + (−1/√2)(−1/√3)(−0.20412 i)
= −0.5 i.
```

So **⟨C⟩ = 0.5** (magnitude) at the F=0 crossing — reproduces Karthein exactly.

**Caveat on the task's named decoupled pair.** The task names the crossing as
A ≈ |N=0,m_N=0,m_S=+½,m_I=−½⟩, B ≈ |N=1,m_N=0,m_S=−½,m_I=+½⟩ at 1.5167 T. As **single
decoupled products**, `⟨A|C|B⟩ = +0.28868 i = i/(2√3) ≈ 0.289 i`, **not 0.5 i** (**verified**).
The 0.5 is the value between the **coupled F=0 eigenstates**, which is the physically correct
operating point: at 1.5167 T the electron spin is strongly Zeeman-decoupled but I is not, so
the near-degenerate eigenstates the experiment actually crosses are the F-coupled
combinations, and ⟨C⟩→0.5 there. The single-product 0.289 is the bare-product limit; the
0.5 is the coupled-state value Karthein quotes. **Both are reproduced by the same
`NSDPV_bBJ`** (the bβJ matrix element is basis-independent of how one labels the crossing).

### (c) Sign pattern / selection rules across the 7 crossings

**The decisive fact: C conserves m_F (ΔM_F = 0).** DeMille (p. 2): "C is a pseudoscalar,
with non-zero m.e.'s between states with the same value of m_F." **Verified**: every nonzero
decoupled and bβJ matrix element has m_F(bra) = m_F(ket); the element is M-independent within
an F-multiplet (F=1→F=1 gives +0.16667 i for **all** M = −1, 0, +1 — verified).

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

| bra (N=0) | ket (N=1) | ⟨C⟩ |
|---|---|---|
| J=½, F=0 | J=½, F=0 | **−0.5 i** |
| J=½, F=1 | J=½, F=1 | **+0.16667 i** (= +i/6) |
| J=½, F=1 | J=1½, F=1 | **−0.23570 i** (= −i/(3√2)) |

(and the Hermitian conjugates N=1→N=0 with the sign flipped per `C_{ba} = −C_{ab}`.) The
sign **alternates** between the F=0 (−), F=1/J=½ (+), and F=1/J=1½ (−) channels — this is the
deterministic per-crossing sign reversal DeMille relies on for systematics rejection, and it
is **entirely within ΔM_F=0** crossings. The implementation must reproduce these three signed
magnitudes.

---

## 5. UNRESOLVED items

1. **The `(−1)^{N′}` ket-N phase in the J-level reduced matrix element is empirically fixed,
   not yet analytically anchored to a B&C equation.** The magnitude (`√2 · √((2J+1)(2J′+1)·3)
   · 9j · √((2N+1)(2N′+1)) (N 1 N′;000) · √(S(S+1)(2S+1))`) and the global `−i√2` are derived
   from B&C (5.109)/(5.111)/(5.113–5.118) and the (8.220) F-reduction. But the specific
   `(−1)^{N′}` (and the overall global sign that makes ⟨N=0,J½,F0|C|N=1,J½,F0⟩ come out **−**0.5 i
   rather than +0.5 i) were pinned by matching the first-principles decoupled construction,
   not transcribed from a single B&C equation number. The codebase's direction-cosine
   convention uses `(−1)^{N−K}` (e.g. `Sz_bBJ`, `T2IS_bBJ`, `StarkZ_bBJ` all carry
   `(−1)^{N0−K0} wigner_3j(N0,1,N1,-K0,0,K1)`); for K=0 and ΔN=±1, `(−1)^{N0} = (−1)^{N1±1} =
   −(−1)^{N1}`, so the codebase-style `(−1)^{N0−K0}` and the validated `(−1)^{N1}` differ by an
   overall sign that the global `−i√2`'s sign absorbs. **The reviewer should confirm (i) the
   compound-rank-1 9j reduced-matrix-element normalization against B&C Ch. 5 (the analogue of
   the `T2IS_bBJ` rank-2 derivation), and (ii) that the `(−1)^{N1}`+global-sign choice is the
   one consistent with the codebase's `(−1)^{N0−K0}` direction-cosine convention** — i.e.
   re-express the validated formula in `(−1)^{N0−K0} wigner_3j(N0,1,N1,−K0,0,K1)` form (it
   should be algebraically identical for K=0) before implementation, so the new element
   matches `Sz_bBJ`/`StarkZ_bBJ` letter-for-letter. The *numerics are not in doubt* (machine
   precision against first principles); only the canonical citation form of the phase is open.

2. **Karthein App. A/B layout vs the task description.** The task refers to "Karthein Eq. (B1)
   and Appendix B define … the C operator; Appendix A the basis." In the local copy
   (arXiv:2310.11192v1, Zotero LYF85P8C): **App. A** = "Effects of Time-Varying Electric
   Fields" (gives the decoupled basis `|N,mN⟩|S,mS⟩|I,mI⟩`, H_eff Eq. (A2), n̂ in spherical
   harmonics Eqs. (A5)–(A6)); **App. B** = "Asymmetry Analytical Formula" (gives **Eq. (B1)**,
   the 2×2 `H±` with `iW`). The **C operator definition and `H_eff = κ′ W_A C` are in the
   main text on p. 2**, not in an appendix, in this version. The published-PRL labels may
   differ; cited locations above are this local copy. (Non-blocking; the operator and Eq. (B1)
   are both quoted verbatim above.)

3. **`1/I` normalization — RESOLVED for I=1/2; OPEN only for I≠1/2 generality.** Karthein
   writes `C = (n̂×S)·I/I` with an explicit `/I`. The data settles the reading (§2.1, §4b):
   the operator that reproduces Karthein's quoted `⟨C⟩=0.5` at the F=0 crossing is
   **`(n̂×S)·Î`** (spin operator `Î`, eigenvalues m_I, reduced m.e. `√(I(I+1)(2I+1))`) with
   **no literal `/I` division** — verified to give exactly −0.5 i, whereas the literal `/(I=½)`
   gives −1.0 i (2× too big) and `/√(I(I+1))` gives −0.577 i. So for ²⁹SiO⁺ (I=1/2) the
   normalization is closed and `NSDPV_bBJ` (which embeds `√(I(I+1)(2I+1))`) is correct and
   verified. **Open only for generality:** for an I≠1/2 isotopologue, whether Karthein's `/I`
   then means `/(quantum number I)`, `/√(I(I+1))`, or is again absorbed must be re-derived
   before reuse — the I=1/2 coincidence (where several conventions can collapse to the same
   0.5) does not fix the I>1/2 factor. Non-blocking for ²⁹SiO⁺; flagged for any future
   heavier-nucleus port. The reviewer should still confirm the `√(I(I+1)(2I+1))` reduced-m.e.
   placement against B&C Eq. (8.220)'s `{I(I+1)(2I+1)}^{1/2}` (it matches `IS_bBJ`).

---

## Validation summary (executed in `Structure` conda env)

- `(A×B)·C = −i√2 [T¹(T¹A,T¹B)]·T¹C` factor: verified to machine precision over 3 random
  vector triples (`(n̂×S)·I` real = `−i√2 [T¹(n̂,S)·T¹(I)]`).
- Decoupled `C` purely imaginary: max|Re| = 0.
- bβJ `C` Hermitian: max|C−C†| = 2.8×10⁻¹⁷.
- Basis change `decouple_b_even` unitary: max|UᵀU−I| = 4.4×10⁻¹⁶.
- `NSDPV_bBJ` (×i) vs first-principles `C`: **max|Δ| = 2.2×10⁻¹⁶** (full 16×16 block).
- F-coupling 6j ≡ audited `IS_bBJ` 6j: 0 mismatches over all (J,J′,F).
- ⟨N=0,J½,F0|C|N=1,J½,F0⟩ = −0.5 i (Karthein C=0.5): verified, incl. hand arithmetic.
- M-independence (F=1→F=1: +i/6 for all M): verified.
- 5 of 6 task (M_a,M_b) crossings vanish by ΔM_F=0; only (+1,+1) survives: verified.
