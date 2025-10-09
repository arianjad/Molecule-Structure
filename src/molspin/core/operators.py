from __future__ import annotations
from typing import Protocol, Mapping, Any, Callable, Dict
import numpy as np
from .quantum_numbers import Basis

class Operator(Protocol):
    name: str
    def matrix(self, basis: Basis, params: Mapping[str, Any]) -> np.ndarray: ...

class RegisteredOperator:
    """Wrap any evaluator: eval_fn(basis, params) -> ndarray"""
    def __init__(self, name: str, eval_fn: Callable[[Basis, Mapping[str, Any]], np.ndarray]):
        self.name = name
        self._eval = eval_fn
    def matrix(self, basis: Basis, params: Mapping[str, Any]) -> np.ndarray:
        return self._eval(basis, params)

class SiteAwareOperator:
    """
    Dispatch to different evaluators based on nuclear spin metadata.
    Expects basis.metadata['nuclear_spins'] = {'M': I_M, 'H': I_H, ...}
    Uses keys: 'even' (I_M==0) vs 'odd' (I_M!=0).
    """
    def __init__(self, name: str, cases: Dict[str, Callable[[Basis, Mapping[str, Any]], np.ndarray]], default_case: str | None = None):
        self.name = name
        self._cases = dict(cases)
        self._default = default_case
    def _select_case(self, basis: Basis) -> str:
        info = getattr(basis, "metadata", {}).get("nuclear_spins", {})
        I_M = info.get("M", None)
        if I_M is not None:
            return "even" if abs(I_M) < 1e-12 else "odd"
        return self._default or ("even" if "even" in self._cases else next(iter(self._cases)))
    def matrix(self, basis: Basis, params: Mapping[str, Any]) -> np.ndarray:
        key = self._select_case(basis)
        return self._cases[key](basis, params)


from .operator_registry import register_operator
from .legacy_bridge import make_operator_from_element, make_siteaware_operator

# ---- Replace the placeholder imports with your real legacy functions ----
# from matrix_elements import (
#     spin_rotation_element_even,   # f(bra_qn, ket_qn, params, meta) -> float
#     spin_rotation_element_odd,
#     fermi_contact_M_even,
#     fermi_contact_M_odd,
#     electron_zeeman_element,      # site-agnostic (if yours is)
# )

# ---- Example placeholders (temporary zeros) so the file imports cleanly ----
def _zero_el(b, k, p, m): return 0.0
spin_rotation_element_even = _zero_el
spin_rotation_element_odd  = _zero_el
fermi_contact_M_even       = _zero_el
fermi_contact_M_odd        = _zero_el
electron_zeeman_element    = _zero_el

# Spin-rotation γ
register_operator(
    make_siteaware_operator(
        "SpinRotation",
        even_element_fn=spin_rotation_element_even,
        odd_element_fn=spin_rotation_element_odd,
    )
)

# Fermi contact b_F (metal-site example; you can add H-site as another operator if needed)
register_operator(
    make_siteaware_operator(
        "FermiContact_M",
        even_element_fn=fermi_contact_M_even,
        odd_element_fn=fermi_contact_M_odd,
    )
)

# Electron Zeeman (often site-agnostic)
register_operator(
    make_operator_from_element(
        "ElectronZeeman",
        electron_zeeman_element,
    )
)
