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

__all__ = ["EffectiveHamiltonian"]
