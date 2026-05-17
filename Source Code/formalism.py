"""Bidirectional N²↔R² spectroscopic-parameter converter (B&C Table 7.2,
quartic-truncated).

The engine's rotational operator is R²-form (N²_op − Λ²·I). Modern papers and
the PGopher default are N²-formulation. A state dict declares which convention
ITS constants are in:
  'formalism': 'N2' | 'R2'
  'Lambda':    |Λ|            (REQUIRED whenever 'formalism' is set; 0 for Σ;
                               never inferred — spec §3, §6, §7)

`convert_formalism` rewrites the constants into the *other* convention — it is
just the (exactly invertible) Table 7.2 arithmetic, so N²→R² and R²→N² are the
same code with one sign — and sets 'formalism' to the new convention so the tag
always follows the data. 'Lambda' is read, never popped, left untouched. A dict
with no 'formalism' key is returned unchanged: no declared convention ⇒ nothing
to convert. That identity-on-absence is what keeps every legacy (untagged, R²)
entry byte-identical through the converter.

The engine wants R², so the load-time call sites convert only when the declared
formalism is 'N2' (molecule_parameters.get_molecule_params and the
Energy_Levels user-dict path). 'formalism'/'Lambda' are left in the returned
dict; the Hamiltonian builders read params by targeted key access, so the extra
keys are inert.

Spec: docs/superpowers/specs/2026-05-15-n2-r2-formalism-converter-design.md
(that spec describes the earlier one-directional, metadata-stripping contract;
this module is now bidirectional and tag-preserving — spec superseded here).
"""
import warnings

# cm⁻¹ ↔ MHz factor. MUST equal molecule_parameters.params_general['c'].
# Passed explicitly by callers to dodge the post-merge 'c'-key collision
# (state dicts reuse 'c' for the hyperfine dipolar constant). Hazard #1.
DEFAULT_C_CM = 29979.2458

# X ↔ its centrifugal-distortion partner X_D (B&C generic-X row). Declarative,
# verified pairs only (spec §4.1): add a pair here, no code change. The
# converter applies the X-row to every pair in this map present in the dict —
# this map is the sole limit on "convert everything convertible".
# Every entry is B&C-verified: the parameter enters the effective Hamiltonian
# as (X + X_D·N²) and is "any molecular parameter other than G, B or D", so it
# obeys Table 7.2's generic-X row X(N²)=X(R²)−Λ²·X_D (B&C p.376/Table 7.2;
# operator eqs. p.374). Eq. numbers below are the (X + X_D·N²) operator forms.
CENTRIFUGAL_PARTNERS = {
    'p+2q':     'p2q_D',     # Λ-doubling p,    B&C eq. 7.190
    'Gamma_SR': 'Gamma_D',   # spin-rotation γ, B&C eq. 7.189
    'q_lD':     'q_lD_D',    # Λ-doubling q,    B&C eq. 7.190
    'ASO':      'A_D',       # spin-orbit A,    B&C eq. 7.187 (verified 2026-05-16)
}

# Sextic (H-order) keys: unsupported (quartic truncation; codebase has no H).
_SEXTIC_KEYS = {'H', 'p2q_H', 'Gamma_H', 'q_lD_H'}


def convert_formalism(params: dict, c_cm: float = DEFAULT_C_CM) -> dict:
    """Convert a parameter dict between the N² and R² conventions.

    Direction follows the dict's own 'formalism' tag: 'N2' → R², 'R2' → N².
    The returned (new) dict carries the flipped 'formalism'; 'Lambda' is read,
    not popped, and left in place. A dict with no 'formalism' is returned
    unchanged. The caller's dict is never mutated.
    """
    out = dict(params)                          # copy; never mutate caller
    formalism = out.get('formalism')            # READ, never pop
    if formalism is None:
        return out                              # no declared convention ⇒ identity
    if formalism not in ('N2', 'R2'):
        raise ValueError(
            f"Unknown 'formalism' {formalism!r}; expected 'N2' or 'R2'.")

    lam = out.get('Lambda')                     # READ, never pop
    if lam is None:
        raise ValueError(
            "a set declaring 'formalism' requires 'Lambda' (|Λ|); set "
            "'Lambda': 0 for Σ. Λ is never inferred (spec §3, §6, §7).")

    sextic = _SEXTIC_KEYS.intersection(out)
    if sextic:
        raise ValueError(
            f"quartic-truncated converter; sextic keys {sorted(sextic)} "
            f"not supported.")

    # Quartic-truncated B&C Table 7.2. D(N²)=D(R²) and X_D(N²)=X_D(R²) are
    # identities at this order, so D / X_D are carried unchanged and serve as
    # the (direction-independent) shift coefficients. With sgn = −1 for
    # N²→R² and +1 for R²→N²:
    #   Be     : Be_out     = Be_in     + sgn·2Λ²·D
    #   X      : X_out       = X_in      − sgn·Λ²·X_D            (every partner)
    #   Origin : Origin_out  = Origin_in − sgn·Λ²·Be_in/c − Λ⁴·D/c
    # Origin is cm⁻¹ while Be/D are MHz ⇒ the 1/c. The Λ⁴D/c term is
    # direction-independent; the pair round-trips to the identity exactly.
    sgn = -1 if formalism == 'N2' else 1
    out['formalism'] = 'R2' if formalism == 'N2' else 'N2'

    lam2 = int(lam) ** 2
    if lam2 == 0:
        return out                              # Σ: convention-independent (tag flipped)

    be_in = out.get('Be')                       # snapshot input Be before B-row
    d = out.get('D')                            # D(N²)=D(R²) at quartic
    if be_in is not None and d is not None:
        out['Be'] = be_in + sgn * 2 * lam2 * d

    for x, xd in CENTRIFUGAL_PARTNERS.items():  # convert every verified pair present
        if x in out and xd in out:
            out[x] = out[x] - sgn * lam2 * out[xd]
        elif xd in out and x not in out:
            warnings.warn(
                f"'{xd}' present without partner '{x}'; no conversion applied "
                f"(likely data error).")

    if 'Origin' in out and be_in is not None:
        d_over_c = (d / c_cm) if d is not None else 0.0
        out['Origin'] = (out['Origin']
                         - sgn * lam2 * be_in / c_cm
                         - (lam2 ** 2) * d_over_c)

    return out
