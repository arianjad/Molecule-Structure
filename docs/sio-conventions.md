# 29SiO+ X2Sigma+ convention gates

Pins every sign/conversion the `molecules['SiO+']['boson']['X0']` model rests on,
with executable checks. Companion script: `Source Code/test_sio_conventions.py`
(exit 0 iff all executable checks pass). Run:

```
conda run -n Structure python "Source Code/test_sio_conventions.py"
```

> **STATUS 2026-06-11: AUDIT COMPLETE — SIGN-OFF.** All inline `PENDING-REVIEWER`
> slots below are RESOLVED in the "B&C audit results" section at the end of this
> file (citations filled; one code bug found and fixed; sole open item = Knight
> 1985 / Decision D5). The inline slot text is kept as the original record.

**Target model** (Karthein, PRL 133, 033003 (2024); local PDF = arXiv:2310.11192v1,
Zotero LYF85P8C — published-PRL equation labels unverified from this copy):
`H_eff = B*N^2 + gamma*N.S + b*I.S + c*(I.n)(S.n)` (Eq. A2, arXiv-v1 p7; the
`D*N^4` centrifugal term appears in the unnumbered main-text `H_0` on p2, NOT in
Eq. A2 — attribution corrected per B&C audit 2026-06-11);
`Hmag = -g_perp*mu_B*S.B - (g_par-g_perp)*mu_B*(S.n)(B.n) - g_I*mu_N*I.B` (Eq. A2).
Karthein Fig.1: positive-parity N=0 RISES with B; crossing with descending N=1 at
B0 ~ 1.5 T.

**Entry under test** (`Source Code/molecule_parameters.py:397-412`, all MHz; muE in
MHz/(V/cm)): `Be=21243.5, D=0.032, Gamma_SR=12., bF=-797., c=-192.,
muE=4.147*0.503412, g_N=-1.1106`.

**Code path** (X2Sigma+, even isotopologue, Hund's case b -> `bBJ` family):
- Builder: `H_even_X` in `Source Code/hamiltonian_builders.py:7-125`.
- Element-name -> matrix-element map: `bBJ_even_X_matrix_elements` in
  `Source Code/molecule_library_class.py:112-137`:
  `'N.S'->SR_bBJ` (115), `'I.S'->IS_bBJ` (122), `'T2_0(I,S)'->T2IS_bBJ` (123),
  `'ZeemanZ'->ZeemanZ_bBJ` (132), `'ZeemanLZ'->ZeemanLZ_bBJ` (134),
  `'ZeemanIZ'->ZeemanIZ_bBJ` (135).
- Matrix-element bodies in `Source Code/matrix_elements.py`:
  `SR_bBJ:115`, `IS_bBJ:122`, `T2IS_bBJ:131`, `ZeemanZ_bBJ:169`,
  `ZeemanIZ_bBJ:176`, `ZeemanLZ_bBJ:183`.
- Constants: `params_general` in `molecule_parameters.py:24-31`:
  `mu_B=1.399624494 MHz/G`, `mu_N=7.62259323e-4 MHz/G`, `g_S=2.0023`, `g_L=1`.

**Baseline reproduced** (Task 1, this session): `initialize_state('SiO+','X',0,
N_list=[0,1],fermion_or_boson='boson',I_nuclei=[0,1/2],P_values=[1/2])` ->
16 states; N=0 block `{-199.25 x3, +597.75}` MHz (splitting 797.0, F=1 triplet
lower); N=1 centered ~42.5 GHz. **High confidence** (ran the diagonalization).

---

## C1 - Hyperfine operator identity (sympy, exact) -- PASS

**Claim.** `bF*(I.S) + (c/3)*sqrt(6)*T2_0(I,S) == b*(I.S) + c*(I.n)(S.n)` with
`b = bF - c/3`, where `T2_0(I,S) = (3 Iz Sz - I.S)/sqrt(6)` in the molecule frame
(n = z). And the codebase actually applies this form.

**Code form.** `hamiltonian_builders.py:45-46`:

```python
H0[i][j] = params['Be']*elements['N^2'] + params['Gamma_SR']*elements['N.S'] + \
    params['bF']*elements['I.S'] + params['c']/3*np.sqrt(6)*elements['T2_0(I,S)']
```

The prefactor on the c term is exactly `c/3*sqrt(6)` multiplying `T2_0(I,S)` (the
rank-2 spherical tensor, mapped to `T2IS_bBJ`, `molecule_library_class.py:123`).

**Check (executed, sympy, I=S=1/2 4x4 matrices, basis |m_S,m_I>):**

```
LHS-RHS zero matrix: True
diff: [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
c-term identity (c/3*sqrt(6)*T20 == c*(IzSz - I.S/3)): True
```

The code's `c/3*sqrt(6)*T2_0(I,S)` term equals `c*(Iz Sz - I.S/3)`, i.e. exactly
the dipolar `c*(I.n)(S.n)` piece after absorbing `-c/3*(I.S)` into the Fermi-contact
coefficient. With `bF = b + c/3` this is the Karthein `b*I.S + c*(I.n)(S.n)`
parameterization. The code's prefactor IS c/3*sqrt(6)-equivalent -> not a FAIL.

**Verdict: PASS** (executable). The Knight-1985 b/c <-> bF/c numerics
(`c = 3t`, `bF = b + c/3` with `t = -64 MHz`, `c = -192 MHz`) are an algebraic
relation among the entry's own comments; the operator identity above is what the
code rests on and it holds exactly.

**PENDING-REVIEWER (citation slot):** the definition `T^2_0(I,S) = (3IzSz-I.S)/sqrt6`
and the case-b(beta-J) reduced matrix element in `T2IS_bBJ` (`matrix_elements.py:131-137`)
against Brown & Carrington -- dipolar hyperfine `T^2(I,S)` matrix element,
Ch. 8 / Ch. 9, eq. ____.

---

## C2 - g_S sign (numeric, executable) -- PASS

**Claim.** Rising N=0 branch has slope `+g_S*mu_B*(1/2) ~ +1.4012 MHz/G`
(electron-spin-dominated), matching Karthein Fig.1 (N=0 positive-parity rises);
descending N=1 branch ~ `-1.4 MHz/G`.

**Method.** `ZeemanMap` (`Energy_Levels.py:292`) over B = 0..200 G (linear regime,
state-tracked via `order=True`); linear-fit each tracked eigenstate's slope.
Builder term: `V_B += g_S*mu_B*ZeemanZ` (`hamiltonian_builders.py:53`).

**Check (executed):**

```
predicted g_S*mu_B*(1/2) = 1.40123 MHz/G  = 14.0123 GHz/T (per level)
rising N=0 max slope   = +1.40166 MHz/G  (14.017 GHz/T)   dev +0.030%
descending N=1 min slope = -1.40166 MHz/G (-14.017 GHz/T)
```

**Unit note (resolves a typo in the task spec).** The task says
"+28 GHz/T = +1.3996 MHz/G". 1 MHz/G = 10 GHz/T, so `1.3996 MHz/G = 14.0 GHz/T`
per level. The "28 GHz/T" is the m_S=+-1/2 *splitting* rate (`g_S*mu_B = 2*mu_B
= 28.03 GHz/T`); the per-level slope is half of that. The gate is checked against
the unambiguous MHz/G figure the task pins ("+1.3996 MHz/G within 5%"). Observed
+1.40166 MHz/G is within 0.03%.

**Verdict: PASS** (executable). Sign of the rising N=0 branch is +, consistent with
Karthein Fig.1 (positive-parity N=0 rises) and `g_S>0`.

---

## C3 - g_N sign + Larmor (numeric, executable) -- PASS (magnitude + sign)

**Claim (magnitude).** Within one electron-spin branch of N=0 at B=15000 G, the
nuclear splitting ~ `|g_N|*mu_N*B = 1.1106 * 7.62259e-4 * 15000 = 12.6985 MHz`.

**Claim (sign).** Code applies `V_B += -g_N*mu_N*ZeemanIZ`
(`hamiltonian_builders.py:57`). For a negative moment (`g_N<0`), the coefficient
`-g_N*mu_N = +|g_N|*mu_N > 0`, so `E(m_I) = -g_N*mu_N*m_I*B` puts `m_I=+1/2`
HIGHER. The entry's `g_N=-1.1106` must reproduce that ordering.

**Method.** Isolate via `get_molecule_params(...)` + an `overrides`-style scratch
dict (`params=...` passed to `initialize_state`): set `bF=c=1e-9` to kill hyperfine
(magnitude check), and additionally `g_S=1e-12` to kill electron Zeeman (sign
check, leaving pure nuclear Zeeman). Read dominant `M (=m_F)` of each eigenstate
from `q_numbers['M']`.

**Check -- magnitude (executed, hyperfine off, electron Zeeman on):**

```
predicted |g_N|*mu_N*B = 12.6985 MHz
N=0 hyperfine-OFF evals @15kG: [-21024.860, -21012.162, 21012.162, 21024.860]
within-branch nuclear splittings: [12.6985, 12.6985]   (0% deviation)
```

(The two electron-spin branches are split by ~42 GHz = g_S*mu_B*B; the 12.70 MHz
splitting sits inside each branch -- unambiguous, not contaminated by bF.)

**Check -- sign (executed, hyperfine + electron Zeeman off, pure nuclear Zeeman):**

```
predicted +-g_N(<0) shift = -g_N*mu_N*(+-1/2)*B = +-6.3492 MHz
highest level: E = +6.3492 MHz, dominant M = +1   (stretched m_I=+1/2)
lowest level : E = -6.3492 MHz, dominant M = -1
```

The `m_I=+1/2` stretched state lies HIGHER (E>0), exactly as predicted for `g_N<0`
under the code's `E = -g_N*mu_N*m_I*B` convention.

**Verdict: PASS** (executable, magnitude + sign). The code's nuclear-Zeeman
convention is `E = -g_N*mu_N*m_I*B` (operator `-g_N*mu_N*ZeemanIZ`); our entry's
`g_N=-1.1106` (the true signed g-factor of 29Si, mu=-0.55529 mu_N) produces the
physically correct ordering. **No sign error in the entry; do NOT flip.**

**PENDING-REVIEWER (citation slot):** the `ZeemanIZ_bBJ` reduced matrix element
(`matrix_elements.py:176-181`) and the overall sign of the nuclear-spin Zeeman
Hamiltonian against Brown & Carrington -- nuclear-spin Zeeman `g_N` / `-g_N mu_N I.B`,
Ch. 8 / Ch. 9, eq. ____. (Confirm B&C's `-g_N mu_N I.B` matches Karthein's
`-g_I mu_N I.B` and the code's `-g_N*mu_N*ZeemanIZ`.)

---

## C4 - gamma (spin-rotation) sign (numeric, executable) -- PASS (code side)

**Claim (code side).** With `gamma=+12 MHz`, the N=1 J=3/2 group sits ABOVE the
J=1/2 group. Case b: `<N.S> = [J(J+1)-N(N+1)-S(S+1)]/2` -> `+1/2` (J=3/2),
`-1` (J=1/2). Splitting `(3/2)*gamma = 18 MHz`. Builder term:
`Gamma_SR*elements['N.S']` (`hamiltonian_builders.py:45`), `'N.S'->SR_bBJ`.

**Method.** N=1-only build with I=0 (`I_nuclei=[0,0]`) to isolate gamma cleanly;
group eigenvalues by `q_numbers['J']`; subtract rotational center.

**Check (executed):**

```
q J list: [0.5, 0.5, 1.5, 1.5, 1.5, 1.5]
evals - center: J=1/2 -> -12.000 MHz (2-fold);  J=3/2 -> +6.000 MHz (4-fold)
split E(J=3/2)-E(J=1/2) = +18.000 MHz  (pred (3/2)*gamma = 18.0)
```

`<N.S>` per level: J=3/2 -> +gamma/2 = +6; J=1/2 -> -gamma = -12. J=3/2 ABOVE
J=1/2, exactly the gamma>0 ordering. The code enters spin-rotation as `+gamma*N.S`.

**Verdict: PASS** (code side, executable).

**PENDING-REVIEWER (literature side).** Zhu et al., JMS 384, 111582 (2022) fit
`gamma_eff = +12 MHz` (`gamma_eff = 0.012(25) GHz`, Tab. 2) with the standard
`+gamma*N.S` convention -- confirm Zhu's sign convention and that their single-state
gamma_eff is the quantity our entry uses (NOT Cameron 1995's deperturbed
gamma=344.16 MHz, a different model; see entry comment, Decision D1). Citation slot:
Zhu et al. 2022, sec/eq ____ (sign convention of gamma).

---

## C5 - g_l / anisotropic Zeeman (audit, partly executable) -- PASS (mechanism) + PENDING-REVIEWER

**What the code implements.** The SiO+ (case-b, X2Sigma+) Zeeman block lives in
`H_even_X`:

```
hamiltonian_builders.py:53   V_B += g_S*mu_B*elements['ZeemanZ']
hamiltonian_builders.py:54-55  if params.get('g_l') is not None:
                                   V_B += g_l*mu_B*elements['ZeemanLZ']
hamiltonian_builders.py:56-57  if params.get('g_N') is not None:
                                   V_B += -g_N*mu_N*elements['ZeemanIZ']
```

`'ZeemanLZ'->ZeemanLZ_bBJ` (`molecule_library_class.py:134`,
`matrix_elements.py:183-189`): an orbital/axis Zeeman matrix element carrying a
`wigner_3j(N0,1,N1,-K0,0,K1)` factor -- "Including coupling to internuclear axis"
(code comment, line 134). For X2Sigma+ (Lambda=0, K=0) this is the case-b axis
Zeeman channel.

(Note: there is a SECOND, distinct anisotropic-Zeeman operator in the code,
`ZeemanPerpZ_even_aBJ` (`matrix_elements.py:551`, comment block 543-550, "B&C
eq. 7.221, parameter g_l"), used by the case-a `H_even_A` 2Pi path. That aBJ g_l
is NOT on the SiO+ case-b path; SiO+ uses `H_even_X` only.)

**g_l is deliberately not set for SiO+.** Confirmed: `'g_l' not in
get_molecule_params('SiO+','X','0','boson')` returns True. The builder's
`params.get('g_l') is not None` guard (`hamiltonian_builders.py:54`) therefore
returns None -> the g_l term is never added. (CORRECTED per B&C audit 2026-06-11:
`H_even_X` never multiplies `ZeemanLZ` by the `params_general` `g_L=1` — with `g_l`
unset, the `ZeemanLZ` operator never enters the SiO+ model at all. The earlier
"part of the validated baseline" reading here was wrong.)

**Check (executed):**

```
'g_l' in SiO+ params: False
N=0,1 build finite spectrum (g_l term skipped): True
```

Build succeeds with a finite 16-level spectrum (Task-1 gates passing is the same
evidence). **Mechanism: PASS** (the `params.get(...) is not None` default-skip path).

**PENDING-REVIEWER (operator identity).** Whether `ZeemanLZ_bBJ`
(`matrix_elements.py:183-189`) is B&C's anisotropic spin-Zeeman form (eq. 9.71
class, or eq. 7.221 for the parallel/perp split) -- citation slot: B&C
anisotropic Zeeman, eq. ____. This is an audit-only question; nothing is gated
because g_l is omitted.

### Knight 1985 backfill -- (g_par - g_perp) slot

The Karthein `Hmag` anisotropic term `-(g_par - g_perp)*mu_B*(S.n)(B.n)` is
currently OMITTED for SiO+ (no g_l / no g_par-g_perp entry). This is consistent
with Karthein's own H_eff usage, where the anisotropic g-tensor enters only at the
~10^-4 level. Recorded as **pending (Decision D5, MIT library pull):** the Knight
et al. 1985 (JACS 107, 2857) neon-matrix ESR g-tensor (g_par, g_perp) for 29Si16O+
would supply `(g_par - g_perp)`. CORRECTED per B&C audit 2026-06-11: a future
backfill must NOT land on the `ZeemanLZ_bBJ` / g_l channel — `ZeemanLZ_bBJ` is
numerically identical to the direction cosine `<n_z>` (rank-1, molecule-frame q=0,
acting on the ROTATIONAL space with S and I spectators; the skeleton of B&C's
orbital Zeeman term (i) of Eq. (9.70), missing the Lambda weight that makes it
vanish for Sigma states). Karthein's `(S.n)(B.n)` acts on S; the correct operator
class is B&C Eq. (9.70) term (iv), the anisotropic spin-Zeeman
`g_l mu_B B_Z sum_{q=+-1} D^(1)_{0q}(omega)* T^1_q(S)` (the case-a analogue is the
codebase's `ZeemanPerpZ_even_aBJ`, "B&C eq. 7.221" — itself not audited). A case-b
element would need to be written/transformed when D5 lands. **Status: PENDING
(citation + decision D5).**

---

## Gate output (this session)

```
[PASS] C1 hyperfine operator identity: LHS-RHS == 0 (4x4): True; c-term == c*(IzSz-I.S/3): True [code prefactor c/3*sqrt(6) verified, hamiltonian_builders.py:46]
[PASS] C2 g_S sign (rising N=0 branch): rise=1.40166 MHz/G (14.017 GHz/T) vs pred 1.40123 (+0.030%); N=1 descend=-1.40166 MHz/G
[PASS] C3 g_N magnitude (Larmor @15kG): |g_N|muN B pred=12.6985 MHz; within-branch splittings=[12.6985, 12.6985]
[PASS] C3 g_N sign (m_I=+1/2 higher for g_N<0): highest level E=6.3492 MHz at dom M=+1.0, lowest at M=-1.0; pred +-6.3492 MHz [V_B += -g_N*mu_N*ZeemanIZ, hamiltonian_builders.py:57]
[PASS] C4 gamma sign (J=3/2 above J=1/2, gamma>0): E(J=3/2)-center=+6.0000, E(J=1/2)-center=-12.0000 MHz; split=18.0000 MHz (pred (3/2)*gamma=18.0)
[PASS] C5 g_l absence tolerated (params.get default): 'g_l' in params: False; build finite spectrum: True [guard: params.get('g_l') is not None, hamiltonian_builders.py:54]

OK: 6 passed, 0 failed
EXIT=0
```

## B&C audit results (spectroscopy-reviewer pass, 2026-06-11) — slots RESOLVED

Method: every element compared against an independent first-principles lab-frame
construction (Gaunt integrals + Clebsch-Gordan, decoupled basis) over all N<=2
pairs, PLUS targeted reads of the B&C page images (PDF page = book page + 32;
sign-bearing equations checked on the printed pages, not text extraction).
B&C copy: Zotero CKZKCGXY. Karthein copy: Zotero LYF85P8C (arXiv:2310.11192v1).

1. **C1 RESOLVED — CONFIRMED.**
   - `IS_bBJ` = B&C Eq. (8.220), p. 440 (restated Eq. (10.125), p. 803), bF-stripped,
     with code index 1=ket=primed; 6j's equal under column/row-pair swaps. Numeric
     max|dev| 2.8e-16.
   - Builder line `bF*I.S + (c/3)*sqrt(6)*T2_0(I,S)` is literally B&C Eq. (9.152),
     p. 662. `T2IS_bBJ` matches B&C Eqs. (8.229)-(8.233), pp. 440-441 /
     (10.129)+(10.131)-(10.133), pp. 803-804 (9j transposition + 3 odd permutations;
     phase residue (-1)^(N+N') = +1 by the 3j(N,2,N';0,0,0) parity rule). Numeric:
     `(sqrt6/3)*T2IS == (I.n)(S.n) - I.S/3` to 2.2e-16 incl. N=1, N=2 diagonals.
   - Which "c": B&C's c IS the Frosch-Foley c (Eqs. (10.133)/(8.237); B&C keeps g_S
     symbolic — "agree with Jette and Cahill if g_S = 2"; 0.1% / ~0.2 MHz on c=-192,
     sub-precision). `bF = b + c/3` verbatim at B&C p. 953.
2. **C3 RESOLVED — convention CONFIRMED; element had a BUG, now FIXED.**
   - Convention: code `-g_N*mu_N*ZeemanIZ` = B&C Eq. (9.70) term (v), p. 620
     (`-g_N mu_N B_Z T^1_{p=0}(I)`; also Eq. (1.39)) = Karthein Eq. (A2)
     `-g_I mu_N I.B`, g_N the SIGNED g-factor. Positive g_N puts m_I=+1/2 LOWER;
     gate measured the converse for g_N=-1.1106. Entry sign correct.
   - **BUG (STOP-the-line, found by the audit): `ZeemanIZ_bBJ` guard lacked
     `kronecker(N0,N1)`.** T^1(I) is rigorously N-diagonal; the formula is
     N-independent, so `<N=0,J=1/2,F=0,M|I_z|N=1,J=1/2,F=1,M>` returned -0.5
     (physical: 0) — spurious OPPOSITE-PARITY couplings ~|g_N| mu_N B ~ 13 MHz at
     1.5 T, a fake parity-violating Zeeman term ~6 orders above the PV signal at
     the Karthein N=0+/N=1- crossing. Fixed 2026-06-11 (delta_NN' added,
     `matrix_elements.py` ZeemanIZ_bBJ guard); repro: dN=1 element now exactly 0,
     dN=0 element unchanged (-0.5). Committed N=0 gates were unaffected (suppressed
     by 2B ~ 42.5 GHz). All gates re-run post-fix: conventions 6/6, Task-1
     instantiation (16 states / 797.0 MHz), RaF tutorial notebook EXIT=0.
3. **C4 RESOLVED — Zhu sign convention CONFIRMED at the source.** Zhu et al., JMS
   384, 111582 (2022) (arXiv:2111.03832, read 2026-06-11 via pdf-mcp): their
   Table 1 (p. 3) lists the case-(b)betaS matrix elements with the STRETCHED
   diagonal `G=1, F=N+1` carrying `+gamma*N/2` — i.e. `+gamma*N.S` with
   `<N.S> = +N/2` for the stretched state: SAME convention as the code's
   `+Gamma_SR*N.S`. Their Eq. (6) (p. 3) defines `t = t0` (the case-(b)betaJ
   dipolar constant), backing the entry's `c = 3t` conversion. Their N=0 note
   (Table 1) gives exactly two elements, `G=1: +bF/4`, `G=0: -3bF/4` — literally
   the Task-1 gate numbers (-199.25/+597.75 MHz). Note: Zhu describe 29SiO+ in
   case (b)betaS (hyperfine >> spin-rotation; G = I+S good at zero field) — a
   BASIS choice; eigenvalues are basis-independent, our bBJ build stands, but
   zero-field labels are better read as (G, F) than (J, F). Their p. 5 refit of
   Cameron's low-N lines with the single-state model gives gamma_X = 0.006(25) GHz,
   consistent with the entry's 12(25) MHz and confirming the Decision-D1 warning
   that Cameron's deperturbed 344.16 MHz is a different model.
4. **C5 RESOLVED — identity established (audit-only).** `ZeemanLZ_bBJ` ==
   direction cosine `<n_z>` exactly (numeric 3.3e-16): the skeleton of B&C
   Eq. (9.70) term (i) (orbital Zeeman) MISSING the Lambda weight — NOT the
   anisotropic spin term (iv) of Eq. (9.70)/(9.71), and not g_r (term iii).
   Disabled for SiO+; backfill pointer corrected above (Knight (g_par-g_perp)
   needs a term-(iv)-class case-b element, to be written when D5 lands).
5. **Knight 1985 backfill — still OPEN (Decision D5,** MIT library pull): the
   `(g_par - g_perp)` g-tensor from Knight et al. JACS 107, 2857 (1985); currently
   omitted, consistent with Karthein at the ~10^-4 level.

B&C citation cautions recorded by the reviewer: the printed closing line of B&C
Eq. (10.123) omits the 1/2 (book typo — cite the 6j line); the apparent
(8.232)/(10.129) phase mismatch resolves on the page images (equivalent forms).

**Audit verdict: SIGN-OFF** (after the ZeemanIZ_bBJ fix and the note corrections
above; sole open item = Knight 1985 / D5, non-blocking for Tasks 3-10).
