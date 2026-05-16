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
