# Spec — Per-state N²/R² formalism converter

**Date:** 2026-05-15
**Branch:** `fix-for-distrib`
**Status:** SUPERSEDED (implemented, then redesigned 2026-05-16)
**Depends on:** `ef06b86` (RaF A0 Be sign fix — the bug this feature makes structurally impossible)

---

> **⚠️ Superseded — read this first.** This spec describes the *as-designed*
> converter: one-directional (N²→R² only), defaulting absent/`'R2'` to a
> pass-through, and **stripping** `formalism`/`Lambda` so the engine never saw
> them. The shipped module diverged on user direction (2026-05-16, commit
> `fc20867`) and is now **bidirectional** (`convert_formalism`: N²→R² *and*
> R²→N², exactly invertible), **tag-following** (`formalism` flipped to follow
> the data, never popped), with `Lambda` read-not-popped and absent-`formalism`
> → identity. Load call sites convert only when `formalism=='N2'` so the engine
> stays R²-fed. `CENTRIFUGAL_PARTNERS` gained the B&C-verified `ASO`/`A_D`
> pair. The §4 Table 7.2 math and §4.1 partner-map design below remain
> accurate; the §3/§5/§6 contract (strip, default-R2-passthrough,
> one-directional) does **not**. Authoritative now:
> [`docs/n2_r2_formalism_converter.md`](../../n2_r2_formalism_converter.md) and
> the `Source Code/formalism.py` module docstring. This file is kept as the
> historical design record, not rewritten.

---

## 1. Problem

`Source Code/` diagonalises an effective Hamiltonian whose **rotational operator is
R²-form** (`Rot_even_aBJ` == B&C eq. 9.138 N²-operator **minus Λ²·I**; machine-verified
to ~1e-8 MHz this session). Spectroscopic constants in modern papers and in PGopher's
**default** output are reported in the **N² formulation** (PGopher follows IUPAC / Hirota
et al. 1994; non-default `RSquaredH=True` switches PGopher to R²).

Today the N²→R² conversion is done **by hand** in `molecule_parameters.py`:

- `RaF/{boson,fermion}/A0` `'Be': 5743.96-2*1.4e-7*c` (B(N²) paper value − 2Λ²D)
- `RaF/{boson,fermion}/A0` `'p+2q': -0.41071*c+1.9e-7*c  # R^2 form` (X(N²) + Λ²·X_D)
- `RaF/{boson,fermion}/A0` `'Origin': 13284.427+0*1350 + 5755.56/c  # …B offset due to
  code being R^2`

Hand arithmetic produced the Be **sign error** fixed in `ef06b86` (`+2Λ²D` vs `−2Λ²D`;
~4 MHz at J~15) **and** a now-confirmed **second bug** in the `Origin` hand term (uses
B_X where the G-row requires Λ²B_A → ~+11.6 MHz on every RaF X–A line, live now; §6.1).
A user entering raw paper (N²) values has no safe path; this converter eliminates the
hand step and structurally fixes both bug classes.

**Goal:** let any state declare its source formalism; convert N²→R² (engine convention)
automatically at load time, implementing **the full B&C Table 7.2 interconversion** so the
rule is general, not a hand-maintained allowlist.

## 2. Architecture decision (locked)

**R² backend + N²→R² converter.** Engine operators **unchanged**. Dual-backend
(N²-native + cross-check) evaluated and **shelved** — it is a full rebuild of `Rot_*`
across bBJ/aBJ/bBS + re-verification of every Π parameter; the actual failure mode here is
parameter-convention confusion, which a converter fixes at far lower risk. Revisit only if
an automatic N²-vs-R² self-check becomes a priority.

## 3. Data model

Two optional keys on any state dict in `molecule_parameters.py`:

| Key | Values | Meaning |
|---|---|---|
| `'formalism'` | `'R2'` (default if absent) \| `'N2'` | Source convention of the constants in this dict |
| `'Lambda'` | int \|Λ\| (0/1/2 = Σ/Π/Δ) | **Declared case-(a) basis Λ label** (see §7) |

- **Default `'R2'` ⇒ every existing entry is byte-identical through the converter**
  (full back-compat). Σ entries (Λ=0) are no-ops regardless of formalism.
- `'Lambda'` is **required** when `formalism=='N2'` and the state is non-Σ; absent ⇒ hard
  error. It is the *declared case-(a) basis label*, never inferred (§7).

### 3.1 RaF A0 migration (reference fixture / regression oracle)

`RaF/boson/A0` and `RaF/fermion/A0` migrate from hand-converted R² values to raw N² +
declaration. Boson:

```python
molecules['RaF']['boson']['A0'] = {
    'formalism': 'N2', 'Lambda': 1,
    'Be': 5743.96,          # B(N²), raw paper value  (was 5743.96-2*1.4e-7*c)
    'D': 1.4e-7*c,
    'p+2q': -0.41071*c,     # X(N²), raw paper value  (was -0.41071*c+1.9e-7*c)
    'p2q_D': 1.9e-7*c,      # X_D
    # ASO, a, d, q, g_lp, g_l, muE, g_S, Origin unchanged (see §6 for Origin)
}
```

Converter **must** reproduce the exact post-`ef06b86` committed R² values:
`Be → 5743.96 − 2·1²·(1.4e-7*c)`; `p+2q → -0.41071*c + 1²·(1.9e-7*c)`
(fermion: `Be=5729.03`, `p+2q=-0.4109*c`). This is the regression oracle.

## 4. Conversion math — full B&C Table 7.2 (verbatim) inverted & quartic-truncated

B&C Table 7.2 (PDF p. 376; book p. 344), *verbatim*, expresses **N² params in terms of
R²** ("X represents any molecular parameter other than G, B or D"):

```
G(N²)   = G(R²) − Λ²B(R²) − Λ⁴D(R²) − Λ⁶H(R²) − ···
B(N²)   = B(R²) + 2Λ²D(R²) + 3Λ⁴H(R²) + ···
D(N²)   = D(R²) + 3Λ²H(R²) + ···
X(N²)   = X(R²) − Λ²X_D(R²) + Λ⁴X_H(R²) + ···
X_D(N²) = X_D(R²) − 2Λ²X_H(R²) + ···
X_H(N²) = X_H(R²) + ···
```

We have **N² inputs**, want **R² outputs** → invert. The codebase has **no sextic `H`**
term, so truncate at quartic: `H = X_H = 0`. Solving the system in dependency order:

```
D(R²)   = D(N²)                               # D row,  H=0
X_D(R²) = X_D(N²)                             # X_D row, X_H=0
B(R²)   = B(N²) − 2Λ²·D(N²)                   # B row,  using D(R²)=D(N²)
X(R²)   = X(N²) + Λ²·X_D(N²)                  # X row,  using X_D(R²)=X_D(N²)
G(R²)   = G(N²) + Λ²·B(N²) − Λ⁴·D(N²)         # G row,  substituting B(R²),D(R²)
```

**This is the complete rule for the entire table.** It is *general*, not an allowlist:

- `D` (and any X_D) — **identity** at quartic order.
- `B` (`Be`) — special rotational row.
- **Every** other parameter X obeys `X(R²)=X(N²)+Λ²·X_D(N²)`. A parameter with **no
  centrifugal partner X_D in the dict** has correction 0 ⇒ **pass-through automatically**.
  No per-parameter enumeration is required for correctness; only the (X ↔ X_D) *pairing*
  must be known.
- `G` (the `Origin` row) — **applied** (§6.1). `Origin` is stored in **cm⁻¹** while
  `Be`/`D` are **MHz**, so the G-row carries a `1/c` unit factor:
  `Origin(R²) = Origin(N²) + Λ²·Be(N²)/c − Λ⁴·D(N²)/c`. The user enters `Origin` as the
  **physical/paper band origin** (N² convention); the converter adds the R² compensation.
  This structurally fixes a confirmed live bug (§6.1).

Fixture check (RaF boson A0, Λ=1): `Be → 5743.96 − 2·(1.4e-7*c)` ✓ post-`ef06b86`;
`p+2q → −0.41071*c + 1·(1.9e-7*c)` (X_D = `p2q_D` = 1.9e-7*c) ✓ current hand value.

### 4.1 Centrifugal-partner map (verified against the actual schema)

Pairs **present in the schema** (`X` ↔ its `X_D`), confirmed by reading
`molecule_parameters.py` + `hamiltonian_builders.py`:

| X | X_D | Evidence |
|---|---|---|
| `p+2q` | `p2q_D` | RaF A0, YbOH A000 (the live fixture) |
| `Gamma_SR` | `Gamma_D` | YbOH X000/X010 carry both; B&C eq. 7.189 `(γ + γ_D N²)` |
| `q_lD` | `q_lD_D` | YbOH/RaOH carry `q_lD`; `q_lD_D` in commented schema example |
| `Be` | `D` | rotational special row (B(R²)=B(N²)−2Λ²D) |

**Resolved non-ambiguities (from data, not guessed):** the `H`/`Yb` suffixes in
`aH/bFH/cH/dH/h1/2H` and `bFYb/cYb/dYb` are **nucleus labels** (two-nucleus hyperfine,
reused from the YbOH backend wart), *not* centrifugal orders → generic-X, no X_D in
schema → pass-through. `q`, `p_lD`, `ASO`, `muE`, `g_S/g_L/g_N`, all hyperfine constants:
no X_D partner in schema → pass-through by the general rule (not special-cased).
`g_l`/`g_lp` are literal expressions on the raw-N² numbers (`-0.41071`, `5743.96`) — they
stay self-consistent under raw-N² entry and are pass-through (no X_D); the earlier
"author-managed exception" is now *entailed by* the general rule, not an exception.

The partner-map is a **declarative module-level dict**; adding a future pair is a one-line
data change, no code change.

## 5. Component & integration

**One pure function** `convert_params_to_engine_R2(params: dict) -> dict` in a new flat
module **`Source Code/formalism.py`** (matches the repo's flat-module convention; isolates
the pure function for unit testing). Logic:

1. `f = params.get('formalism','R2')`; `f=='R2'` → return shallow copy unchanged.
2. `f=='N2'`: `L2 = params.get('Lambda',0)**2`.
   - `L2==0` (Σ) → copy with bookkeeping keys stripped (no-op conversion).
   - else: `Be → Be − 2·L2·D` (if `Be`,`D` present); for each `(X,X_D)` in the
     partner-map with **both** keys present: `X → X + L2·X_D`; `Origin → Origin +
     L2·Be(N²)/c − L2²·D(N²)/c` (G-row, if `Origin` present — §6.1; uses pre-conversion
     `Be`/`D` and the `1/c` unit factor). `D` and every `X_D` unchanged.
3. Strip `'formalism'`,`'Lambda'`; **all other keys pass through untouched**.

**Two call sites, one function:** inside `get_molecule_params()`
([molecule_parameters.py] — `Energy_Levels.py:15,68`) immediately before return, **and**
on the user-supplied-dict path in `MoleculeLevels.__init__`
([Energy_Levels.py:67-70](Source Code/Energy_Levels.py#L67)). Engine,
`hamiltonian_builders.py`, `matrix_elements.py`: **zero changes**.

## 6. Scope & error handling

**In scope v1 (full Table 7.2 general rule):** rotational `Be`/`D`; the universal
`X→X+Λ²X_D` rule over the §4.1 partner-map; Σ identity; full back-compat for `'R2'`;
RaF A0 boson+fermion migrated to raw-N².

**Errors:**

- `formalism=='N2'`, non-Σ, `'Lambda'` absent → `ValueError` (never silently assume Λ=0
  — the exact silent-error class this feature kills).
- `formalism` ∉ `{'R2','N2'}` → `ValueError`.
- Any sextic `H`-type key on an `N2` dict → `ValueError` ("quartic-truncated converter").
- `X_D` present without its `X` partner on an `N2` dict → warn (likely data error).

### 6.1 `Origin` / G-row — CONFIRMED BUG; converter fixes it structurally

**Trace (verified this session).** `Origin` is read **nowhere** in `Source Code/`; it is
consumed in the notebook spectrum generators. `RaF_X_A_BR_diagnostic.ipynb`, verbatim:

```python
E_evals = A.eigensystem(Ez, Bz)[0] / c_cm
E_evals = E_evals + A.parameters['Origin']
nu = np.subtract.outer(E_evals, G_evals)        # ν = (E_A^{R²}/c + Origin) − E_X^{R²}/c
```

**Mechanism.** The R²-form operator is `N²_op − Λ²·I` (machine-verified), so
`Be·(N²_op − Λ²·I)` lowers **every level of a manifold** by `Be·Λ²` vs. the physical/N²
convention. X²Σ⁺ (Λ_X=0): no shift. A²Π (Λ_A=1): the entire A manifold sits `Be_A`
lower. For `ν` above to equal the physical value, `Origin` must add it back:
`Origin = Origin_phys + Be_A·Λ_A²/c`. This is exactly the Table 7.2 G-row
(`Origin(R²)=Origin(N²)+Λ²Be(N²)/c−Λ⁴D/c`) — the addend is the state's **own** Be.

**Bug.** The committed hand term is `+ 5755.56/c` = **B_X** (the *X-state* Be), not
B_A = 5743.96. Error = (5755.56 − 5743.96)/c ≈ **+11.6 MHz constant offset on every RaF
X–A transition frequency**, live in `RaF X-A.ipynb` and the diagnostic. Same class as the
Be bug (`ef06b86`); confirmed by two independent derivations (operator algebra + B&C
Table 7.2 G-row), both giving form *and* the correct state's-own-Be coefficient.

**Audit (all Π A-states, this session).** Active bug is **contained to RaF A0**:
boson `+5755.56/c` vs Be_A=5743.96 ⇒ **+11.6 MHz**; fermion same `+5755.56/c` vs
Be_A=**5729.03** ⇒ **+26.5 MHz**. YbOH/CaOH/DyO A-states carry **no `Origin` key**
(latent gap, not an active bug — covered once they declare `N2`+`Lambda`+`Origin`).
BaOH X000/X010 are Σ (`Origin:341.6`, no `/c`) — unaffected. "BaF X-A.ipynb" actually
runs `molecule_name='RaF', params=None` (no BaF params exist; BaF commented out). Dead
`# # …5726.48/c` YbOH-174 block: inactive, cleanup-only.

**Decision:** the G-row **is in scope** and the converter **fixes this bug structurally**.
Migrated entries store `Origin` as the **physical/paper band origin** (N²; e.g. the bare
`13284.427` + ASO bookkeeping, **no hand `+B/c`**); the converter applies
`Origin(R²) = Origin(N²) + Λ²·Be(N²)/c − Λ⁴·D(N²)/c` using the state's own declared Λ
and its `Be`/`D`. Caveats: (i) the `1/c` unit factor (Origin cm⁻¹, Be/D MHz) is part of
the row and must be unit-tested; (ii) the X²Σ ground state has Λ=0 ⇒ zero G-row shift, so
the A−X difference composes correctly with no ground-side term; a future non-Σ ground
state gets its own G-row and still composes via `E_A−E_X`; (iii) ASO/Π₁/₂ band-origin
bookkeeping (`+0*1350`, the `13284.427` base) is physical input the user owns — the
converter only adds the R² G-row compensation.

## 7. Λ is a *basis label*, not a good quantum number — case-(c) validity (verified)

B&C §7.5.3 / Table 7.2 are derived "for an electronic state |η, Λ⟩" with a definite Λ
(eq. 7.197: `R = Nₓi + Nᵧj + (N_z − L_z)k`; Table-7.2 scaling *is* L_z²=Λ²). The N²-vs-R²
distinction is therefore **only defined for a definite-Λ (case a/b) basis state**.

**Boundary (verified against B&C's actual case-(c) treatment, Ch 10 pp. 819–822 / eq.
10.144–10.146):** in a *native* Hund's case-(c)/Ω representation, where Λ is not a good
quantum number, B&C does **not** use N² or §7.5.3-R² at all — it writes
`Hrot ∝ R² = (J − Jₐ)²` in an `|v(Ω),Jₐ,J,Ω⟩` basis, absorbing `Bᵥ·Jₐ(Jₐ+1)` into the
term value (scaling with Jₐ(Jₐ+1), **not Λ²**). There is no Λ²-scaled N²↔R²
interconversion in that representation — structurally inapplicable, not merely omitted.
**The converter is undefined for a state modelled in a case-(c)/Ω basis** and must not be
applied there.

**Why this is nonetheless safe here:** the codebase engine builds the Hamiltonian in a
**Hund's case-(a) basis** (`matrix_elements.py`; every basis function carries a declared,
definite Λ). A physically case-(c) molecule is a *diagonalised superposition of
definite-Λ case-(a) basis functions*; the Table-7.2 reparametrisation acts on the
**per-manifold input constants scaled by that manifold's declared Λ²** (a fixed input
integer), never on the mixed eigenstate. So within this engine the converter is exact
regardless of physical coupling case. The precondition — a definite declared Λ — is
exactly what the §6 hard error enforces; the converter must never infer Λ from
eigenstates, and `'formalism':'N2'` is undefined (not "use basis-Λ") for any future
native case-(c)/Ω representation.

## 8. Testing / verification

1. **Regression oracle (`Be`, `p+2q`):** converter on RaF boson **and** fermion A0
   raw-N² input reproduces post-`ef06b86` R² values exactly (these hand values were
   *correct*, so byte-equality is the oracle). **`Origin` is the exception — its
   committed hand value is BUGGY** (§6.1), so the converter must **not** reproduce it;
   the oracle for `Origin` is the *physically corrected* value
   `Origin_phys + Be_A(N²)/c − D/c`, validated against an independent RaF X–A line /
   PGopher (see §10), not byte-equality.
2. **Back-compat sweep:** every existing entry (default/absent formalism) is
   dict-equal through the converter.
3. **Σ no-op:** `formalism:'N2', Lambda:0` ⇒ pass-through.
4. **Generality:** synthetic dict with each §4.1 pair converts by `X+Λ²X_D`; a param
   with no X_D partner is identity; `D`/`X_D` unchanged.
5. **`Origin` G-row:** N² dict with `Origin`+`Be`+`Lambda` ⇒
   `Origin += Λ²·Be/c − Λ⁴·D/c`; Σ (Λ=0) ⇒ `Origin` unchanged; unit factor `1/c`
   explicitly asserted. Cross-check: the corrected RaF A0 `Origin` differs from the
   buggy committed value by exactly `(Be_X − Be_A)/c ≈ 11.6 MHz`.
6. **Physics regression:** build RaF A²Π via the real `MoleculeLevels` path with the
   migrated N² entry; eigenvalue spacings match the pre-migration (post-bug-fix) R² entry
   to ~1e-8 MHz (reuse this session's verified harness).
7. **Error paths:** N²+non-Σ+no-Lambda, bad formalism, sextic key → each `ValueError`.
8. **Project gate (mandatory, no `--allow-errors`):**
   `conda run -n Structure jupyter execute "Jupyter Notebooks/RaF_Calcs_Tutorial.ipynb"`
   on a fresh kernel.

## 9. Out of scope (explicit)

Dual N²-native backend; sextic `H` terms; any change to engine
operators/builders/matrix elements; rewriting notebook spectrum generators (they keep
consuming `params['Origin']` — the converter changes only the stored value); molecules
other than RaF A0 stay at default `'R2'` until their source conventions are individually
confirmed. (Note: `Origin`/G-row conversion is now **in scope** — §6.1.)

## 10. Decisions & open items

**Resolved (user directive + verification):**
- **Q1 module** → new flat `Source Code/formalism.py`.
- **Q2 keys** → `'formalism'` / `'Lambda'` (no existing convention-key idiom to match).
- **Q3 scope** → **full Table 7.2 general rule** (user: "consider N²→R² for all variable
  types in the B&C table"), implemented as the universal `X+Λ²X_D` rule + special
  `Be`/`D` rows; `Origin`/G the one evidence-backed exception (§6.1).

**Flagged for the user / follow-up:**
1. **CONFIRMED second bug (resolved by this design):** the `Origin` hand term uses
   `+B_X/c` where the G-row requires `+Λ²B_A/c` — a ~+11.6 MHz constant offset on all
   RaF X–A lines, live now (§6.1). The converter fixes it. **Physics validation needed:**
   the corrected RaF A0 `Origin` and representative **N=1** X–A line(s) are validated
   against the published RaF X & A constants (PGopher / N² formalism) and line positions
   in **García Ruiz et al., *Nat. Phys.* (2024), s41567-023-02296-w**
   (`https://web.mit.edu/lns/news/s41567-023-02296-w.pdf`; constants table + N=1 lines in
   main figures / supplement). This is the regression oracle for the `Origin` G-row.
2. Confirm `Gamma_SR↔Gamma_D` and `q_lD↔q_lD_D` units/quartic-order partnering against
   B&C eq. 7.189-7.190 during implementation (fixture only covers `p+2q↔p2q_D`/`Be↔D`).
