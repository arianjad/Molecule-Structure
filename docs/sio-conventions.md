# 29SiO+ X2Sigma+ convention gates

Pins every sign/conversion the `molecules['SiO+']['boson']['X0']` model rests on,
with executable checks. Companion script: `Source Code/test_sio_conventions.py`
(exit 0 iff all executable checks pass). Run:

```
conda run -n Structure python "Source Code/test_sio_conventions.py"
```

**Target model** (Karthein, PRL 133, 033003 (2024), Eq. A2):
`H_eff = B*N^2 + D*N^4 + gamma*N.S + b*I.S + c*(I.n)(S.n)`;
`Hmag = -g_perp*mu_B*S.B - (g_par-g_perp)*mu_B*(S.n)(B.n) - g_I*mu_N*I.B`.
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
returns None -> the g_l term is never added. (`params_general` does carry `g_L=1`,
the orbital g-factor for `ZeemanLZ`; for Lambda=0 this contributes through the
K=0 3-j and is part of the validated baseline, not the anisotropic g_l term.)

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
would supply `(g_par - g_perp)`. If added later, it would enter the case-b model
via the `ZeemanLZ_bBJ` / g_l channel (the `H_even_X` line 54-55 term), the case-b
analogue of Karthein's `(S.n)(B.n)` anisotropy. **Status: PENDING (citation +
decision).**

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

## Summary of PENDING-REVIEWER slots (Brown & Carrington pass)

1. **C1** -- dipolar hyperfine `T^2(I,S)` reduced matrix element / `T2_0(I,S) =
   (3IzSz-I.S)/sqrt6` definition (`matrix_elements.py:T2IS_bBJ:131`), B&C eq. ____.
2. **C3** -- nuclear-spin Zeeman matrix element + overall sign `-g_N mu_N I.B`
   (`matrix_elements.py:ZeemanIZ_bBJ:176`), B&C eq. ____; cross-check vs Karthein
   `-g_I mu_N I.B`.
3. **C4** (literature side) -- Zhu et al. JMS 384, 111582 (2022) sign convention of
   `gamma_eff = +12 MHz`, sec/eq ____.
4. **C5** -- `ZeemanLZ_bBJ` (`matrix_elements.py:183`) as B&C's anisotropic Zeeman
   form, eq. ____ (audit only; g_l omitted).
5. **Knight 1985 backfill** -- `(g_par - g_perp)` g-tensor from Knight et al. JACS
   107, 2857 (1985), pending MIT library pull (Decision D5); currently omitted,
   consistent with Karthein at the 10^-4 level.
