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
