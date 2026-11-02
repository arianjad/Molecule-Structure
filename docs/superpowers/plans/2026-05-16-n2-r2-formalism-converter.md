# N²/R² Formalism Converter Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers-extended-cc:subagent-driven-development (recommended) or superpowers-extended-cc:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a load-time converter that rewrites N²-formulation spectroscopic constants into the engine's R² convention per the full B&C Table 7.2, eliminating hand arithmetic and structurally fixing the confirmed RaF A0 `Origin` bug.

**Architecture:** One pure function `convert_params_to_engine_R2(params, c_cm)` in a new flat module `Source Code/formalism.py`, invoked at the two parameter-entry points (`get_molecule_params` return; `MoleculeLevels.__init__` user-dict path). Engine, builders, and matrix elements are untouched. Default `formalism='R2'` ⇒ exact pass-through (full back-compat).

**Tech Stack:** Python 3 (conda env `Structure`), no external test framework (plain-`assert` test scripts run via `python`), `numpy` only in the physics-regression task. Spec: `docs/superpowers/specs/2026-05-15-n2-r2-formalism-converter-design.md`.

---

## Critical implementation hazards (read before Task 0)

1. **The `'c'` key collision.** `params_general['c'] = 29979.2458` is the cm⁻¹↔MHz factor, but state dicts reuse `'c'` as the **hyperfine dipolar constant** (e.g. `RaF/boson/X0` `'c': 19`, `A0` `'c': 0`). `get_molecule_params` does `merged.update(state_params)`, so after merge `merged['c']` is the *hyperfine* value, not the unit factor. The G-row's `1/c` MUST use the module-global `c` (=`params_general['c']`), passed explicitly as `c_cm`. Never read the unit factor from the params dict.
2. **Snapshot N² values before overwriting.** The G-row needs `Be(N²)`/`D(N²)` (pre-conversion). `Be` is overwritten by the B-row earlier in the same call — capture `be_n2`/`d_n2` first.
3. **Exactly one conversion per path.** `get_molecule_params` converts the DB path. `MoleculeLevels.__init__` must convert ONLY the user-supplied-dict branch (add an `else:`), never re-convert the `params is None` branch (which already went through `get_molecule_params`).
4. **Spec §5-vs-§6 wording reconciliation.** §5 step 2 loosely says `params.get('Lambda',0)`; §6's error contract is authoritative: `formalism='N2'` with **no** `'Lambda'` ⇒ `ValueError` (never default to 0 — that is the silent-error class this feature kills). Σ states must set `'Lambda': 0` explicitly.
5. **`Origin` is NOT diagonalized.** It is a downstream A−X placement scalar (notebooks). Within-A-manifold eigenvalue *spacings* depend only on `Be`/`p+2q` (byte-exact post-conversion), so the physics regression (Task 5) validates Be/p+2q; `Origin` correctness is validated separately (Tasks 4, 6).

---

## File Structure

- **Create** `Source Code/formalism.py` — the converter + declarative partner-map. One responsibility: N²→R² parameter rewriting. No engine/numpy deps.
- **Create** `Source Code/test_formalism.py` — plain-`assert` test script (no pytest dependency), runnable via `conda run -n Structure python "Source Code/test_formalism.py"`. Grows task-by-task.
- **Modify** `Source Code/molecule_parameters.py` — import converter; replace `return merged` (line 391) with converted return; migrate `RaF/boson/A0` and `RaF/fermion/A0` to raw-N² + declaration.
- **Modify** `Source Code/Energy_Levels.py` — convert the user-supplied-dict branch (after line 70).
- **Create** `Source Code/test_formalism_physics.py` — MoleculeLevels eigenvalue-spacing regression (Structure env).
- **Create** `docs/superpowers/plans/2026-05-16-n2-r2-formalism-converter-natphys-oracle.md` — extracted Nat. Phys. validation numbers (Task 6 artifact).

---

### Task 0: Module scaffold + core dispatch (R2 pass-through, formalism errors, bookkeeping strip)

**Goal:** `Source Code/formalism.py` exists with `convert_params_to_engine_R2` handling the no-op/error paths; test harness established.

**Files:**
- Create: `Source Code/formalism.py`
- Test: `Source Code/test_formalism.py`

**Acceptance Criteria:**
- [ ] `formalism='R2'` (and absent) ⇒ shallow copy, value-equal, `formalism`/`Lambda` stripped.
- [ ] Unknown `formalism` ⇒ `ValueError`.
- [ ] `formalism='N2'` with no `'Lambda'` ⇒ `ValueError`.
- [ ] Caller's dict never mutated.

**Verify:** `conda run -n Structure python "Source Code/test_formalism.py"` → prints `OK: 4 passed` and exits 0.

**Steps:**

- [ ] **Step 1: Write the failing test**

Create `Source Code/test_formalism.py`:

```python
"""Plain-assert tests for formalism.convert_params_to_engine_R2 (no pytest dep)."""
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
from formalism import convert_params_to_engine_R2

_passed = 0
def check(name, fn):
    global _passed
    fn()
    _passed += 1
    print(f"  pass: {name}")

def expect_valueerror(fn):
    try:
        fn()
    except ValueError:
        return
    raise AssertionError("expected ValueError, none raised")

def test_r2_passthrough_and_strip():
    src = {'formalism': 'R2', 'Be': 100.0, 'p+2q': -3.0}
    out = convert_params_to_engine_R2(src)
    assert out == {'Be': 100.0, 'p+2q': -3.0}, out
    assert 'formalism' not in out and 'Lambda' not in out

def test_absent_formalism_is_r2():
    out = convert_params_to_engine_R2({'Be': 7.0})
    assert out == {'Be': 7.0}, out

def test_unknown_formalism_raises():
    expect_valueerror(lambda: convert_params_to_engine_R2({'formalism': 'X2'}))

def test_n2_without_lambda_raises():
    expect_valueerror(lambda: convert_params_to_engine_R2(
        {'formalism': 'N2', 'Be': 1.0}))

def test_caller_dict_not_mutated():
    src = {'formalism': 'R2', 'Be': 1.0}
    convert_params_to_engine_R2(src)
    assert src == {'formalism': 'R2', 'Be': 1.0}, src

if __name__ == '__main__':
    check('r2_passthrough_and_strip', test_r2_passthrough_and_strip)
    check('absent_formalism_is_r2', test_absent_formalism_is_r2)
    check('unknown_formalism_raises', test_unknown_formalism_raises)
    check('n2_without_lambda_raises', test_n2_without_lambda_raises)
    check('caller_dict_not_mutated', test_caller_dict_not_mutated)
    print(f"OK: {_passed} passed")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `conda run -n Structure python "Source Code/test_formalism.py"`
Expected: FAIL — `ModuleNotFoundError: No module named 'formalism'`.

- [ ] **Step 3: Write minimal implementation**

Create `Source Code/formalism.py`:

```python
"""N²→R² spectroscopic-parameter converter (B&C Table 7.2, quartic-truncated).

The engine's rotational operator is R²-form (N²_op − Λ²·I). Constants from
modern papers / default PGopher are N²-formulation. A state dict may declare:
  'formalism': 'N2'   (default 'R2' → exact pass-through)
  'Lambda':    |Λ|     (REQUIRED for any N2 state; 0 for Σ; never inferred)
This module rewrites N² constants into the engine's R² convention at load time.
Spec: docs/superpowers/specs/2026-05-15-n2-r2-formalism-converter-design.md
"""
import warnings

# cm⁻¹ ↔ MHz factor. MUST equal molecule_parameters.params_general['c'].
# Passed explicitly by callers to dodge the post-merge 'c'-key collision
# (state dicts reuse 'c' for the hyperfine dipolar constant). Hazard #1.
DEFAULT_C_CM = 29979.2458

# X ↔ its centrifugal-distortion partner X_D (B&C generic-X row). Declarative:
# add a future verified pair here, no code change (spec §4.1).
CENTRIFUGAL_PARTNERS = {
    'p+2q':     'p2q_D',
    'Gamma_SR': 'Gamma_D',
    'q_lD':     'q_lD_D',
}

# Sextic (H-order) keys: unsupported (quartic truncation; codebase has no H).
_SEXTIC_KEYS = {'H', 'p2q_H', 'Gamma_H', 'q_lD_H'}


def convert_params_to_engine_R2(params: dict, c_cm: float = DEFAULT_C_CM) -> dict:
    """Return a new dict with constants in the engine's R² convention.

    'formalism'/'Lambda' are consumed (stripped) so the engine never sees them.
    """
    out = dict(params)                       # shallow copy; never mutate caller
    formalism = out.pop('formalism', 'R2')
    lam = out.pop('Lambda', None)

    if formalism == 'R2':
        return out
    if formalism != 'N2':
        raise ValueError(
            f"Unknown 'formalism' {formalism!r}; expected 'N2' or 'R2'.")

    if lam is None:
        raise ValueError(
            "formalism='N2' requires 'Lambda' (|Λ|); set 'Lambda': 0 for Σ "
            "states. Λ is never inferred (spec §3, §6, §7).")

    return out  # full conversion added in Tasks 1–2
```

- [ ] **Step 4: Run test to verify it passes**

Run: `conda run -n Structure python "Source Code/test_formalism.py"`
Expected: `OK: 5 passed`, exit 0.

- [ ] **Step 5: Commit**

```bash
cd "Source Code" && /Library/Developer/CommandLineTools/usr/bin/git add formalism.py test_formalism.py && /Library/Developer/CommandLineTools/usr/bin/git commit -m "feat(formalism): scaffold N2->R2 converter (dispatch + error paths)"
```

---

### Task 1: Rotational + generic-X conversion (B/D rows, partner-map, Σ no-op, sextic guard)

**Goal:** Implement the `Be`/`D` row, the universal `X→X+Λ²X_D` rule over the partner-map, Σ identity, and the sextic-key guard.

**Files:**
- Modify: `Source Code/formalism.py`
- Test: `Source Code/test_formalism.py`

**Acceptance Criteria:**
- [ ] `Be(R²) = Be(N²) − 2Λ²·D(N²)`; `D` and every `X_D` unchanged.
- [ ] For each `(X,X_D)` in the partner-map present together: `X → X + Λ²·X_D`.
- [ ] Param with no `X_D` partner ⇒ identity (general rule, not allowlist).
- [ ] `formalism='N2', Lambda=0` (Σ) ⇒ identity (only bookkeeping stripped).
- [ ] Sextic key on an N2 dict ⇒ `ValueError`.
- [ ] `X_D` present without its `X` partner ⇒ `warnings.warn`, no conversion of that pair.

**Verify:** `conda run -n Structure python "Source Code/test_formalism.py"` → `OK: 11 passed`.

**Steps:**

- [ ] **Step 1: Write the failing tests** — append to `Source Code/test_formalism.py` before the `if __name__` block:

```python
import warnings as _w

def test_sigma_is_identity():
    src = {'formalism': 'N2', 'Lambda': 0, 'Be': 100.0, 'p+2q': -3.0,
           'p2q_D': 0.05}
    out = convert_params_to_engine_R2(src)
    assert out == {'Be': 100.0, 'p+2q': -3.0, 'p2q_D': 0.05}, out

def test_be_d_row():
    out = convert_params_to_engine_R2(
        {'formalism': 'N2', 'Lambda': 1, 'Be': 100.0, 'D': 0.5})
    assert out['Be'] == 100.0 - 2 * 1 * 0.5, out
    assert out['D'] == 0.5, out                     # D unchanged

def test_generic_x_row_lambda2():
    out = convert_params_to_engine_R2(
        {'formalism': 'N2', 'Lambda': 2, 'p+2q': -10.0, 'p2q_D': 0.25})
    assert out['p+2q'] == -10.0 + 4 * 0.25, out     # Λ²=4
    assert out['p2q_D'] == 0.25, out                # X_D unchanged

def test_all_partners_and_passthrough():
    src = {'formalism': 'N2', 'Lambda': 1,
           'Gamma_SR': 7.0, 'Gamma_D': 0.1,
           'q_lD': -4.0, 'q_lD_D': 0.2,
           'ASO': 5e6, 'muE': 1.23}             # ASO/muE: no X_D → identity
    out = convert_params_to_engine_R2(src)
    assert out['Gamma_SR'] == 7.0 + 1 * 0.1, out
    assert out['q_lD'] == -4.0 + 1 * 0.2, out
    assert out['ASO'] == 5e6 and out['muE'] == 1.23, out

def test_sextic_key_raises():
    expect_valueerror(lambda: convert_params_to_engine_R2(
        {'formalism': 'N2', 'Lambda': 1, 'Be': 1.0, 'p2q_H': 1e-9}))

def test_orphan_xd_warns():
    with _w.catch_warnings(record=True) as rec:
        _w.simplefilter('always')
        out = convert_params_to_engine_R2(
            {'formalism': 'N2', 'Lambda': 1, 'p2q_D': 0.05})
    assert any('p2q_D' in str(r.message) for r in rec), rec
    assert out['p2q_D'] == 0.05, out
```

Add these to the `__main__` runner block:

```python
    check('sigma_is_identity', test_sigma_is_identity)
    check('be_d_row', test_be_d_row)
    check('generic_x_row_lambda2', test_generic_x_row_lambda2)
    check('all_partners_and_passthrough', test_all_partners_and_passthrough)
    check('sextic_key_raises', test_sextic_key_raises)
    check('orphan_xd_warns', test_orphan_xd_warns)
```

- [ ] **Step 2: Run to verify failure**

Run: `conda run -n Structure python "Source Code/test_formalism.py"`
Expected: FAIL on `test_be_d_row` (Be not converted — current stub returns `out` unchanged).

- [ ] **Step 3: Implement** — replace the final `return out  # full conversion added in Tasks 1–2` in `formalism.py` with:

```python
    sextic = _SEXTIC_KEYS.intersection(out)
    if sextic:
        raise ValueError(
            f"formalism='N2' converter is quartic-truncated; sextic keys "
            f"{sorted(sextic)} not supported.")

    lam2 = int(lam) ** 2
    if lam2 == 0:
        return out                            # Σ: convention-independent

    be_n2 = out.get('Be')                     # snapshot N² values (hazard #2)
    d_n2 = out.get('D')

    if be_n2 is not None and d_n2 is not None:
        out['Be'] = be_n2 - 2 * lam2 * d_n2   # B row

    for x, xd in CENTRIFUGAL_PARTNERS.items():            # generic-X row
        if x in out and xd in out:
            out[x] = out[x] + lam2 * out[xd]
        elif xd in out and x not in out:
            warnings.warn(
                f"formalism='N2': '{xd}' present without partner '{x}'; "
                f"no conversion applied (likely data error).")

    return out                                # G row added in Task 2
```

- [ ] **Step 4: Run to verify pass**

Run: `conda run -n Structure python "Source Code/test_formalism.py"`
Expected: `OK: 11 passed`.

- [ ] **Step 5: Commit**

```bash
cd "Source Code" && /Library/Developer/CommandLineTools/usr/bin/git add formalism.py test_formalism.py && /Library/Developer/CommandLineTools/usr/bin/git commit -m "feat(formalism): B/D + generic-X Table 7.2 rows, Sigma no-op, sextic guard"
```

---

### Task 2: G-row / `Origin` conversion with the 1/c unit factor

**Goal:** Implement `Origin(R²) = Origin(N²) + Λ²·Be(N²)/c − Λ⁴·D(N²)/c`, using snapshotted N² `Be`/`D` and the explicit `c_cm`.

**Files:**
- Modify: `Source Code/formalism.py`
- Test: `Source Code/test_formalism.py`

**Acceptance Criteria:**
- [ ] N2 dict with `Origin`+`Be`+`Lambda=1` ⇒ `Origin += Be(N²)/c − D(N²)/c`.
- [ ] `Lambda=2` ⇒ `Origin += 4·Be/c − 16·D/c` (Λ²=4, Λ⁴=16).
- [ ] Σ (`Lambda=0`) ⇒ `Origin` unchanged.
- [ ] `Origin` present but `Be` absent ⇒ `Origin` unchanged (no spurious shift).
- [ ] Custom `c_cm` argument is honored.

**Verify:** `conda run -n Structure python "Source Code/test_formalism.py"` → `OK: 16 passed`.

**Steps:**

- [ ] **Step 1: Write the failing tests** — append before `if __name__`:

```python
def test_g_row_lambda1():
    c = 29979.2458
    out = convert_params_to_engine_R2(
        {'formalism': 'N2', 'Lambda': 1,
         'Be': 5743.96, 'D': 1.4e-7 * c, 'Origin': 13284.427})
    exp = 13284.427 + 1 * 5743.96 / c - 1 * (1.4e-7 * c) / c
    assert abs(out['Origin'] - exp) < 1e-12, (out['Origin'], exp)

def test_g_row_lambda2_powers():
    c = 1000.0
    out = convert_params_to_engine_R2(
        {'formalism': 'N2', 'Lambda': 2, 'Be': 10.0, 'D': 2.0,
         'Origin': 100.0}, c_cm=c)
    exp = 100.0 + 4 * 10.0 / c - 16 * 2.0 / c          # Λ²=4, Λ⁴=16
    assert abs(out['Origin'] - exp) < 1e-12, (out['Origin'], exp)

def test_g_row_sigma_unchanged():
    out = convert_params_to_engine_R2(
        {'formalism': 'N2', 'Lambda': 0, 'Be': 50.0, 'Origin': 9.0})
    assert out['Origin'] == 9.0, out

def test_g_row_no_be_no_shift():
    out = convert_params_to_engine_R2(
        {'formalism': 'N2', 'Lambda': 1, 'Origin': 9.0})
    assert out['Origin'] == 9.0, out

def test_g_row_custom_c():
    out = convert_params_to_engine_R2(
        {'formalism': 'N2', 'Lambda': 1, 'Be': 300.0, 'Origin': 0.0},
        c_cm=100.0)
    assert abs(out['Origin'] - 3.0) < 1e-12, out      # 300/100, no D
```

Add to runner:

```python
    check('g_row_lambda1', test_g_row_lambda1)
    check('g_row_lambda2_powers', test_g_row_lambda2_powers)
    check('g_row_sigma_unchanged', test_g_row_sigma_unchanged)
    check('g_row_no_be_no_shift', test_g_row_no_be_no_shift)
    check('g_row_custom_c', test_g_row_custom_c)
```

- [ ] **Step 2: Run to verify failure**

Run: `conda run -n Structure python "Source Code/test_formalism.py"`
Expected: FAIL on `test_g_row_lambda1` (`Origin` unchanged by current code).

- [ ] **Step 3: Implement** — replace `return out  # G row added in Task 2` with:

```python
    if 'Origin' in out and be_n2 is not None:             # G row (spec §6.1)
        d_over_c = (d_n2 / c_cm) if d_n2 is not None else 0.0
        out['Origin'] = (out['Origin']
                         + lam2 * be_n2 / c_cm
                         - (lam2 ** 2) * d_over_c)

    return out
```

- [ ] **Step 4: Run to verify pass**

Run: `conda run -n Structure python "Source Code/test_formalism.py"`
Expected: `OK: 16 passed`.

- [ ] **Step 5: Commit**

```bash
cd "Source Code" && /Library/Developer/CommandLineTools/usr/bin/git add formalism.py test_formalism.py && /Library/Developer/CommandLineTools/usr/bin/git commit -m "feat(formalism): G-row Origin conversion with 1/c unit factor"
```

---

### Task 3: Wire the two call sites + back-compat sweep

**Goal:** Invoke the converter at `get_molecule_params` return and the `MoleculeLevels` user-dict branch; prove every existing entry is unchanged (default R2).

**Files:**
- Modify: `Source Code/molecule_parameters.py:391`
- Modify: `Source Code/Energy_Levels.py` (after line 70)
- Test: `Source Code/test_formalism.py`

**Acceptance Criteria:**
- [ ] `get_molecule_params(...)` returns converted params, passing the module-global `c` as `c_cm` (hazard #1).
- [ ] `MoleculeLevels.__init__` converts ONLY the user-supplied-dict branch (hazard #3 — no double conversion).
- [ ] Back-compat: every `(molecule, spin, state)` in `molecules` is value-equal before vs. after conversion (all currently default R2).

**Verify:** `conda run -n Structure python "Source Code/test_formalism.py"` → `OK: 17 passed` (includes the sweep).

**Steps:**

- [ ] **Step 1: Write the failing test** — append before `if __name__`:

```python
def test_backcompat_all_entries_unchanged():
    import importlib
    mp = importlib.import_module('molecule_parameters')
    n = 0
    for mol, spins in mp.molecules.items():
        for spin, states in spins.items():
            for st, d in states.items():
                if not isinstance(d, dict):
                    continue
                before = dict(d)
                after = convert_params_to_engine_R2(dict(d), c_cm=mp.c)
                # all current entries are default-R2 ⇒ identical
                assert after == before, (mol, spin, st,
                    set(before) ^ set(after))
                n += 1
    assert n > 0
    print(f"    swept {n} entries")
```

Add to runner: `check('backcompat_all_entries_unchanged', test_backcompat_all_entries_unchanged)`

- [ ] **Step 2: Run to verify failure**

Run: `conda run -n Structure python "Source Code/test_formalism.py"`
Expected: FAIL — this passes trivially *only* once entries stay R2; it will currently PASS (no migration yet). To make it a true red: temporarily it is a guard for Task 4. Run it now and confirm `swept N entries` printed and `OK: 17 passed` (it is green by construction pre-migration; it becomes the regression guard once RaF A0 migrates in Task 4).

> Note for the implementer: this test has no red phase pre-migration (all entries are R2). Its purpose is the Task-4 regression guard. Proceed to Step 3 (wiring) and re-run.

- [ ] **Step 3: Wire call site 1** — in `Source Code/molecule_parameters.py`, add after the imports near the top (after line 3, the existing `from sympy...` import):

```python
from formalism import convert_params_to_engine_R2
```

Then replace line 391 `    return merged` with:

```python
    return convert_params_to_engine_R2(merged, c_cm=c)
```

(`c` here is the module global defined at line 32 `c=params_general['c']` — the unit factor, NOT `merged['c']`; hazard #1.)

- [ ] **Step 4: Wire call site 2** — in `Source Code/Energy_Levels.py`, add to the imports (next to line 15 `from molecule_parameters import get_molecule_params`):

```python
from formalism import convert_params_to_engine_R2
```

Replace the block at lines 67–70:

```python
        if params is None:
            params = get_molecule_params(molecule_name,elec_state,vib_state,fermion_or_boson)
        elif not isinstance(params,dict):
            raise ValueError('Params must be a dictionary of molecular parameters, see molecule_parameters.py for examples')
```

with:

```python
        if params is None:
            params = get_molecule_params(molecule_name,elec_state,vib_state,fermion_or_boson)
            # already converted inside get_molecule_params (hazard #3)
        elif not isinstance(params,dict):
            raise ValueError('Params must be a dictionary of molecular parameters, see molecule_parameters.py for examples')
        else:
            # user-supplied dict: convert here (the only place this path is converted)
            params = convert_params_to_engine_R2(params)
```

- [ ] **Step 5: Run to verify pass + smoke the import wiring**

Run: `conda run -n Structure python "Source Code/test_formalism.py"`
Expected: `OK: 17 passed`, with `swept N entries` printed (N ≈ 20+).

Run: `conda run -n Structure python -c "import sys; sys.path.insert(0,'Source Code'); from molecule_parameters import get_molecule_params; print(sorted(get_molecule_params('RaF','A','0','boson'))[:5])"`
Expected: a key list printed, no exception (RaF A0 is still R2 → unchanged).

- [ ] **Step 6: Commit**

```bash
cd "Source Code" && /Library/Developer/CommandLineTools/usr/bin/git add formalism.py test_formalism.py molecule_parameters.py Energy_Levels.py && /Library/Developer/CommandLineTools/usr/bin/git commit -m "feat(formalism): wire converter into get_molecule_params + MoleculeLevels user-dict path; back-compat sweep"
```

---

### Task 4: Migrate RaF A0 (boson + fermion) to raw-N²; regression oracle for Be/p+2q and corrected Origin

**Goal:** Replace hand arithmetic in both RaF A0 dicts with raw paper N² values + `formalism`/`Lambda`; prove `Be`/`p+2q` reproduce the (correct) committed values and `Origin` is the *corrected* (non-buggy) value.

**Files:**
- Modify: `Source Code/molecule_parameters.py` (`RaF/boson/A0` lines 50–69; `RaF/fermion/A0` lines 86–105)
- Test: `Source Code/test_formalism.py`

**Acceptance Criteria:**
- [ ] Boson: converted `Be == 5743.96 - 2*1.4e-7*c`; `p+2q == -0.41071*c + 1.9e-7*c` (byte-equal to pre-migration committed values).
- [ ] Fermion: converted `Be == 5729.03 - 2*1.4e-7*c`; `p+2q == -0.4109*c + 1.9e-7*c`.
- [ ] Boson `Origin` == `13284.427 + 5743.96/c - 1.4e-7` and differs from the OLD buggy committed `13284.427 + 5755.56/c` by exactly `(5743.96-5755.56)/c` (≈ −3.87e-4 cm⁻¹ ≡ −11.6 MHz).
- [ ] Fermion `Origin` differs from OLD buggy by `(5729.03-5755.56)/c` (≈ −26.5 MHz).
- [ ] `g_lp`/`g_l` literal expressions unchanged (pass-through; reference raw-N² numbers).

**Verify:** `conda run -n Structure python "Source Code/test_formalism.py"` → all pass incl. new oracle; back-compat sweep still green (RaF A0 now legitimately differs only because it declares N2 — see Step 1 note).

**Steps:**

- [ ] **Step 1: Write the failing oracle test** — append before `if __name__`:

```python
def test_raf_a0_regression_oracle():
    import importlib
    mp = importlib.import_module('molecule_parameters')
    importlib.reload(mp)
    c = mp.c
    # boson
    b = mp.get_molecule_params('RaF', 'A', '0', 'boson')
    assert abs(b['Be'] - (5743.96 - 2*1.4e-7*c)) < 1e-9, b['Be']
    assert abs(b['p+2q'] - (-0.41071*c + 1.9e-7*c)) < 1e-9, b['p+2q']
    old_buggy_b = 13284.427 + 5755.56/c
    corrected_b = 13284.427 + 5743.96/c - 1.4e-7
    assert abs(b['Origin'] - corrected_b) < 1e-9, b['Origin']
    assert abs((b['Origin'] - old_buggy_b) - (5743.96-5755.56)/c) < 1e-9
    assert 'formalism' not in b and 'Lambda' not in b
    # fermion
    f = mp.get_molecule_params('RaF', 'A', '0', 'fermion')
    assert abs(f['Be'] - (5729.03 - 2*1.4e-7*c)) < 1e-9, f['Be']
    assert abs(f['p+2q'] - (-0.4109*c + 1.9e-7*c)) < 1e-9, f['p+2q']
    old_buggy_f = 13284.427 + 5755.56/c
    corrected_f = 13284.427 + 5729.03/c - 1.4e-7
    assert abs(f['Origin'] - corrected_f) < 1e-9, f['Origin']
    assert abs((f['Origin'] - old_buggy_f) - (5729.03-5755.56)/c) < 1e-9
```

Add to runner: `check('raf_a0_regression_oracle', test_raf_a0_regression_oracle)`

Also update `test_backcompat_all_entries_unchanged`: RaF A0 now declares `formalism='N2'`, so it legitimately changes. Replace its assert body's loop guard with a skip for the two migrated entries:

```python
                if (mol, spin, st) in {('RaF','boson','A0'),
                                       ('RaF','fermion','A0')}:
                    continue   # migrated to N2 — covered by oracle test
```

- [ ] **Step 2: Run to verify failure**

Run: `conda run -n Structure python "Source Code/test_formalism.py"`
Expected: FAIL on `test_raf_a0_regression_oracle` — current `Origin` is the buggy `+5755.56/c`, `formalism` absent.

- [ ] **Step 3: Migrate `RaF/boson/A0`** — replace lines 50–69 of `Source Code/molecule_parameters.py`:

```python
molecules['RaF']['boson']['A0'] = {
    'formalism': 'N2',   # raw paper (N²) values; converter → engine R² (spec)
    'Lambda': 1,         # A²Π, declared case-(a) basis Λ
    'Be': 5743.96,       # B(N²) paper value
    'ASO': 1350*c,       # Fixed from 1350 cm^-1
    'h1/2': 0,           # Calc relayed from Silviu, h1/2 = a-(bf+2c/3)=A||/2
    'a':19/2,
    'bF': 0,
    'c': 0,
    'd': -9,             # Silviu; PGopher convention, B+C needs minus sign
    'p+2q': -0.41071*c,  # X(N²) paper value
    'q':0,
    'D': 1.4e-7*c,
    'p2q_D': 1.9e-7*c,   # X_D for p+2q
    'g_lp': -0.41071/(2*5743.96/c),#-0.865,
    'g_l': -6000/(2*5743.96),
    'muE': 1*0.503412,
    'g_S': 2.0023,
    'Origin': 13284.427, # PHYSICAL band origin (N²); converter adds R² G-row.
                         # NO hand +B/c term (that was the +11.6 MHz bug).
    }
```

- [ ] **Step 4: Migrate `RaF/fermion/A0`** — replace lines 86–105:

```python
molecules['RaF']['fermion']['A0'] = {
    'formalism': 'N2',
    'Lambda': 1,
    'Be': 5729.03,       # B(N²) paper value
    'ASO': 1350*c,
    'h1/2Yb': -2852/2,   # From Skripnikov theory
    'dYb': -1*-0.076*c,  # Wilkins expt; pgopher convention → add - sign
    'h1/2H': 19/2,
    'aH':19/2,
    'bFH': 0,
    'cH': 0,
    'dH': -9,
    'p+2q': -0.4109*c,   # X(N²) paper value
    'q':0,
    'e2Qq0': 0,
    'D': 1.4e-7*c,
    'p2q_D': 1.9e-7*c,
    'g_lp': -0.4109/(2*0.191100),#-0.865,
    'muE': 1.0*0.503412,
    'g_S': 2.0023,
    'Origin': 13284.427, # PHYSICAL band origin (N²); converter adds R² G-row
    }
```

- [ ] **Step 5: Run to verify pass**

Run: `conda run -n Structure python "Source Code/test_formalism.py"`
Expected: all checks pass incl. `raf_a0_regression_oracle`; sweep still green.

- [ ] **Step 6: Commit**

```bash
cd "Source Code" && /Library/Developer/CommandLineTools/usr/bin/git add molecule_parameters.py test_formalism.py && /Library/Developer/CommandLineTools/usr/bin/git commit -m "fix(molecule_parameters): migrate RaF A0 boson+fermion to raw-N2; converter fixes the Origin +11.6/+26.5 MHz hand bug"
```

---

### Task 5: Physics regression — eigenvalue spacings unchanged vs. pre-migration R²

**Goal:** Building RaF A²Π through the real `MoleculeLevels` path with the migrated entry yields eigenvalue spacings identical (≤1e-8 MHz) to the pre-migration post-bugfix R² constants (validates Be/p+2q conversion through the engine).

**Files:**
- Create: `Source Code/test_formalism_physics.py`

**Acceptance Criteria:**
- [ ] For RaF A0 boson and fermion, the sorted within-manifold eigenvalue *differences* from the migrated+converted params match those from an explicit pre-migration R² dict to ≤1e-8 MHz.

**Verify:** `conda run -n Structure python "Source Code/test_formalism_physics.py"` → prints `OK physics: boson max|Δ|=… fermion max|Δ|=…` with both < 1e-8.

**Steps:**

- [ ] **Step 1: Write the regression** — create `Source Code/test_formalism_physics.py`:

```python
"""RaF A²Π eigenvalue-spacing regression: migrated-N² (via converter) vs the
explicit pre-migration post-ef06b86 R² constants. Spacings are Origin-free."""
import os, sys, numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from Energy_Levels import MoleculeLevels
from molecule_parameters import c, params_general

def _old_r2(spin):
    be0 = 5743.96 if spin == 'boson' else 5729.03
    d = {'Be': be0 - 2*1.4e-7*c, 'ASO': 1350*c, 'p+2q': -0.41071*c + 1.9e-7*c
         if spin == 'boson' else -0.4109*c + 1.9e-7*c, 'q': 0,
         'D': 1.4e-7*c, 'p2q_D': 1.9e-7*c, 'muE': 1*0.503412, 'g_S': 2.0023}
    if spin == 'boson':
        d.update({'h1/2': 0, 'a': 19/2, 'bF': 0, 'c': 0, 'd': -9,
                  'g_lp': -0.41071/(2*5743.96/c), 'g_l': -6000/(2*5743.96)})
    else:
        d.update({'h1/2Yb': -2852/2, 'dYb': -1*-0.076*c, 'h1/2H': 19/2,
                  'aH': 19/2, 'bFH': 0, 'cH': 0, 'dH': -9, 'e2Qq0': 0,
                  'g_lp': -0.4109/(2*0.191100)})
    md = dict(params_general); md.update(d)
    return md

def _evals(params, spin):
    m = MoleculeLevels.initialize_state(
        molecule_name='RaF', elec_state='A', vib_state='0',
        fermion_or_boson=spin, N_range=[1, 2, 3], round=8,
        params=params, P_values=[1/2, 3/2])
    return np.sort(m.eigensystem(0.0, 1e-9)[0])

ok = True
for spin in ('boson', 'fermion'):
    old = _evals(_old_r2(spin), spin)                 # explicit R² (pre-migration)
    new = _evals(None, spin)                           # migrated N² via converter
    n = min(len(old), len(new))
    d_old = np.diff(old[:n]); d_new = np.diff(new[:n])
    m = float(np.max(np.abs(d_old - d_new))) if n > 1 else 0.0
    print(f"  {spin}: max|Δspacing| = {m:.3e} MHz  ({n} levels)")
    ok = ok and (m < 1e-8)
print("OK physics" if ok else "FAIL physics")
sys.exit(0 if ok else 1)
```

> Implementer note: if `initialize_state`'s signature differs, inspect `Source Code/Energy_Levels.py:29` (`def initialize_state`) and adjust kwargs; keep `params=` (explicit dict ⇒ user-dict path) vs `params=None` (DB path) — both must route through the converter.

- [ ] **Step 2: Run to verify it discriminates** — temporarily perturb `_old_r2` boson `be0` to `5743.96 + 1.0` and run:

Run: `conda run -n Structure python "Source Code/test_formalism_physics.py"`
Expected: `FAIL physics` (max|Δ| ≫ 1e-8). Then revert the perturbation.

- [ ] **Step 3: Run the true regression**

Run: `conda run -n Structure python "Source Code/test_formalism_physics.py"`
Expected: `OK physics`, both spins max|Δ| < 1e-8.

- [ ] **Step 4: Commit**

```bash
cd "Source Code" && /Library/Developer/CommandLineTools/usr/bin/git add test_formalism_physics.py && /Library/Developer/CommandLineTools/usr/bin/git commit -m "test(formalism): RaF A2Pi eigenvalue-spacing regression vs pre-migration R2"
```

---

### Task 6: External validation — corrected `Origin` vs. published Nat. Phys. RaF constants

**Goal:** Confirm the *corrected* RaF A0 `Origin` and a regenerated N=1 X–A line agree with the published PGopher-N² constants / line positions in Udrescu, Wilkins, … García Ruiz et al., *Nat. Phys.* 2024 (DOI 10.1038/s41567-023-02296-w) — i.e., the fix is physically right, not just internally consistent.

**Files:**
- Create: `docs/superpowers/plans/2026-05-16-n2-r2-formalism-converter-natphys-oracle.md`
- Create: `Source Code/test_formalism_natphys.py`

**Acceptance Criteria:**
- [ ] The paper + Supplementary Information are read; the RaF X²Σ⁺ and A²Π₁/₂ PGopher-N² constants (Tₑ/band origin, B, ASO, p+2q, D) and at least one explicit N=1 X–A transition frequency are transcribed into the oracle markdown with page/table citations.
- [ ] The migrated RaF A0 (via `get_molecule_params`) `Be`/`p+2q`/`Origin` match the paper's N² constants under the documented unit conversions (`Origin` compared as the physical band origin, i.e. *before* the R² G-row, since the paper is N²).
- [ ] A regenerated N=1 X–A line (built through `MoleculeLevels` + the notebook `ν = (E_A/c + Origin) − E_X/c` relation) agrees with the paper's measured line to within the paper's stated uncertainty (or, if only constants are tabulated, the corrected `Origin` matches Tₑ within rounding).

**Verify:** `conda run -n Structure python "Source Code/test_formalism_natphys.py"` → prints each comparison and `OK natphys` (exit 0).

**Steps:**

- [ ] **Step 1: Extract the oracle numbers.** The PDF is already cached at `/tmp/raf_natphys.pdf` (main article, 14 pp.); the **Supplementary Information** holds the PGopher constants table. Fetch it: the SI is linked from the article page `https://www.nature.com/articles/s41567-023-02296-w` (section "Supplementary information"). Download the SI PDF, then use `mcp__pdf-mcp__pdf_search`/`pdf_read_pages` to locate the "spectroscopic constants" / "PGopher" table and any tabulated N=1 (e.g. P(1)/Q(1)/R(0)) line. Record verbatim, with citations, into `docs/superpowers/plans/2026-05-16-n2-r2-formalism-converter-natphys-oracle.md` using this structure:

```markdown
# Nat. Phys. 2024 (s41567-023-02296-w) — RaF N² validation oracle

Source: Udrescu, Wilkins, … García Ruiz et al., Nature Physics (2024),
DOI 10.1038/s41567-023-02296-w. Constants are PGopher / N² formulation.

## X²Σ⁺ (v=0)   [cite: SI Table __, p.__]
- B(N²)   = … MHz (or cm⁻¹ → note conversion)
- D       = …
- γ (Gamma_SR) = …

## A²Π₁/₂ (v=0)  [cite: SI Table __, p.__]
- Tₑ / band origin = … cm⁻¹
- B(N²)   = …
- A_SO    = …
- p+2q    = …
- D       = …

## Representative N=1 X–A line  [cite: SI/Fig __, p.__]
- Assignment: … ;  measured ν = … cm⁻¹ (± …)
```

If the SI is not machine-readable for a number, record "not tabulated — validate via constants only" and proceed with the constants comparison.

- [ ] **Step 2: Write the validation test** — create `Source Code/test_formalism_natphys.py`. Fill the `PAPER_*` constants from the oracle markdown (exact values + the documented MHz↔cm⁻¹ handling). Skeleton with real assertions (no placeholders in logic — only the transcribed numbers come from Step 1):

```python
"""Validate migrated RaF A0 against the Nat. Phys. 2024 N² constants.
Numbers transcribed in docs/.../natphys-oracle.md (cite there)."""
import os, sys, numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from molecule_parameters import get_molecule_params, c

# --- transcribed from the oracle markdown (Step 1) ---
PAPER_A_Be_N2_MHz   = ...   # A²Π B(N²)
PAPER_A_p2q_N2_MHz  = ...   # A²Π p+2q (N²)
PAPER_BAND_ORIGIN_CM = ...  # A²Π₁/₂–X²Σ⁺ Tₑ / band origin, cm⁻¹
TOL_CONST_MHz = 5.0         # tighten to the paper's quoted precision

b = get_molecule_params('RaF', 'A', '0', 'boson')
# Be is stored R²; back out N² for the paper comparison: Be(N²)=Be(R²)+2Λ²D
be_n2 = b['Be'] + 2*1*b['D']
p2q_n2 = b['p+2q'] - 1*b['p2q_D']
# Origin is R² (physical + Λ²Be/c − Λ⁴D/c); back out the physical band origin:
origin_phys = b['Origin'] - 1*be_n2/c + 1*(b['D'])/c

print(f"  Be(N²)   converter={be_n2:.4f}  paper={PAPER_A_Be_N2_MHz:.4f}")
print(f"  p+2q(N²) converter={p2q_n2:.4f} paper={PAPER_A_p2q_N2_MHz:.4f}")
print(f"  origin_phys={origin_phys:.6f} cm⁻¹  paper={PAPER_BAND_ORIGIN_CM:.6f}")

ok = (abs(be_n2 - PAPER_A_Be_N2_MHz) < TOL_CONST_MHz
      and abs(p2q_n2 - PAPER_A_p2q_N2_MHz) < TOL_CONST_MHz
      and abs(origin_phys - PAPER_BAND_ORIGIN_CM) < (TOL_CONST_MHz / c))
print("OK natphys" if ok else "FAIL natphys")
sys.exit(0 if ok else 1)
```

If the SI tabulates an explicit N=1 line, extend with a `MoleculeLevels` build of X(N=1)/A(J=1/2) and the notebook `ν = (E_A/c + Origin) − E_X/c` relation, asserting agreement within the paper's uncertainty.

- [ ] **Step 3: Run the validation**

Run: `conda run -n Structure python "Source Code/test_formalism_natphys.py"`
Expected: `OK natphys`. If `FAIL`, do NOT loosen `TOL` to force green — investigate (most likely a unit or N²/R² back-out error) and report; a real disagreement is a finding, not a test to silence.

- [ ] **Step 4: Commit**

```bash
cd "Source Code" && /Library/Developer/CommandLineTools/usr/bin/git add test_formalism_natphys.py "../docs/superpowers/plans/2026-05-16-n2-r2-formalism-converter-natphys-oracle.md" && /Library/Developer/CommandLineTools/usr/bin/git commit -m "test(formalism): external validation vs Nat. Phys. 2024 RaF N2 constants"
```

---

### Task 7: Mandatory completion gate — tutorial notebook on a fresh kernel

**Goal:** The project's single completeness check passes with the converter live.

**Files:** none (verification only)

**Acceptance Criteria:**
- [ ] `Jupyter Notebooks/RaF_Calcs_Tutorial.ipynb` executes end-to-end on a fresh kernel with **no** `--allow-errors`, exit 0.
- [ ] Full unit suite green: `test_formalism.py`, `test_formalism_physics.py`, `test_formalism_natphys.py`.

**Verify:** the two commands below, both exit 0.

**Steps:**

- [ ] **Step 1: Run the full unit/physics/validation suite**

```bash
conda run -n Structure python "Source Code/test_formalism.py" && \
conda run -n Structure python "Source Code/test_formalism_physics.py" && \
conda run -n Structure python "Source Code/test_formalism_natphys.py"
```
Expected: `OK: N passed`, `OK physics`, `OK natphys`; combined exit 0.

- [ ] **Step 2: Run the mandatory project gate (fresh kernel, no --allow-errors)**

```bash
conda run -n Structure jupyter execute "Jupyter Notebooks/RaF_Calcs_Tutorial.ipynb"
```
Expected: exit 0, no cell errors. (Per memory `feedback_no_allow_errors` — never add `--allow-errors`; a failure here is a real failure to fix, not to suppress.)

- [ ] **Step 3: Final commit**

```bash
cd "Source Code" && /Library/Developer/CommandLineTools/usr/bin/git add -A && /Library/Developer/CommandLineTools/usr/bin/git commit -m "test(formalism): full suite + tutorial gate green; N2->R2 converter complete" --allow-empty
```

---

## Self-Review

**1. Spec coverage**

| Spec section | Task |
|---|---|
| §3 data model (`formalism`/`Lambda`, R2 default, Λ required) | 0, 1 |
| §3.1 RaF A0 migration + regression oracle | 4 |
| §4 full Table 7.2 (D/X_D identity, B row, generic-X, G row) | 1, 2 |
| §4.1 partner-map (declarative dict; resolved non-ambiguities) | 1 |
| §5 one pure fn + two call sites, engine untouched | 0–3 |
| §6 errors (bad formalism, missing Λ, sextic, orphan X_D) | 0, 1 |
| §6.1 Origin G-row bug fix (1/c factor, physical-origin entry) | 2, 4 |
| §7 case-(c) boundary (Λ = declared basis label; never inferred) | 0 (hard error on missing Λ enforces precondition) |
| §8.1 Be/p+2q byte-oracle; Origin corrected-not-buggy | 4 |
| §8.2 back-compat sweep | 3 |
| §8.3 Σ no-op | 1, 2 |
| §8.4 generality/passthrough | 1 |
| §8.5 Origin G-row + 11.6 MHz cross-check | 2, 4 |
| §8.6 physics regression ≤1e-8 MHz | 5 |
| §8.7 error paths | 0, 1 |
| §8.8 tutorial gate (no --allow-errors) | 7 |
| §10 flagged item 1 (Nat. Phys. external validation) | 6 |
| §10 flagged item 2 (confirm Gamma/q_lD partnering) | 6 oracle markdown notes; pairs already in §4.1 map, exercised in Task 1 |

No spec requirement is unmapped. Out-of-scope items (§9: dual backend, sextic H, engine edits, non-RaF molecules, notebook rewrites) are not implemented by design.

**2. Placeholder scan:** The only intentional fill-ins are the *transcribed paper numbers* in Task 6 (`PAPER_*`), which cannot be known until the SI is read — Step 1 produces them with citations. All code/logic is complete; no "TBD/handle errors/similar-to" placeholders.

**3. Type consistency:** `convert_params_to_engine_R2(params: dict, c_cm: float=DEFAULT_C_CM) -> dict` — identical signature in Tasks 0–6. `CENTRIFUGAL_PARTNERS`/`_SEXTIC_KEYS`/`DEFAULT_C_CM` names consistent. Call sites pass `c_cm=mp.c` (DB path) or default (user-dict path) consistently.

---

## Notes for the executor

- Run from repo root unless a step `cd`s. macOS git: use `/Library/Developer/CommandLineTools/usr/bin/git` (xcrun shim stalls — see user `rules/repl.md`).
- The Drive/network-FS git flush trap: keep `add && commit` chained in one shell invocation (already written that way).
- Module edits leave a running kernel stale — the Task 7 fresh-kernel re-execution is the verification, not autoreload.
- If `MoleculeLevels.initialize_state` kwargs differ from Task 5's call, read `Source Code/Energy_Levels.py:29` and adapt; the invariant is: `params=None` ⇒ DB path, explicit `params=dict` ⇒ user-dict path, both routed through the converter exactly once.
