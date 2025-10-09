from __future__ import annotations
from typing import Mapping, Any, List, Tuple
import numpy as np
from .operator_registry import get_operator
from .quantum_numbers import Basis

class EffectiveHamiltonian:
    """Composable effective Hamiltonian: a list of (operator_name, params)."""
    def __init__(self, basis: Basis, terms: List[Tuple[str, Mapping[str, Any]]]):
        self.basis = basis
        self.terms = list(terms)
    def matrix(self) -> np.ndarray:
        H = np.zeros((len(self.basis), len(self.basis)))
        for name, pars in self.terms:
            H += get_operator(name).matrix(self.basis, pars)
        return H
    

from typing import Mapping, Any, List, Tuple

def X2Sigma(params: Mapping[str, Any]) -> List[Tuple[str, Mapping[str, Any]]]:
    """
    Minimal X²Σ model as a list of (operator_name, params) pairs.
    Expected params include (examples): gamma, bF_M, gS, Bz (etc.)
    """
    terms: List[Tuple[str, Mapping[str, Any]]] = []

    if "gamma" in params:
        terms.append(("SpinRotation", {"gamma": params["gamma"]}))

    if "bF_M" in params:
        terms.append(("FermiContact_M", {"bF": params["bF_M"]}))

    # Zeeman (electron)
    gS = params.get("gS", 2.0023)
    if "B" in params or "Bz" in params:
        # pass the whole params dict if your element function reads components
        terms.append(("ElectronZeeman", {"gS": gS, **{k: v for k, v in params.items() if k.startswith("B")}}))

    return terms


__all__ = ["EffectiveHamiltonian","X2Sigma"]
