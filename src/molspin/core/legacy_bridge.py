# src/molspin/core/legacy_bridge.py
from __future__ import annotations
from typing import Callable, Mapping, Any, Dict
import numpy as np
from .quantum_numbers import Basis
from .operators import RegisteredOperator, SiteAwareOperator

ElementFn = Callable[[Mapping[str, Any], Mapping[str, Any], Mapping[str, Any], Dict[str, Any]], float]

def make_operator_from_element(name: str, element_fn: ElementFn):
    """
    Wrap a legacy element function f(bra_qn, ket_qn, params, meta)->float
    into a RegisteredOperator that builds the full matrix by looping i,j.
    """
    def eval_(basis: Basis, params: Mapping[str, Any]) -> np.ndarray:
        n = len(basis)
        M = np.zeros((n, n), dtype=float)
        meta = basis.metadata
        # naive dense assembly (clarity first; optimize later)
        for i, bra in enumerate(basis.states):
            qi = bra.values
            for j, ket in enumerate(basis.states):
                qj = ket.values
                M[i, j] = element_fn(qi, qj, params, meta)
        return M
    return RegisteredOperator(name, eval_)

def make_siteaware_operator(
    name: str,
    even_element_fn: ElementFn,   # used when I_M == 0
    odd_element_fn: ElementFn,    # used when I_M != 0
):
    return SiteAwareOperator(
        name,
        cases={
            "even": lambda basis, p: make_operator_from_element(name+"[even]", even_element_fn).matrix(basis, p),
            "odd":  lambda basis, p: make_operator_from_element(name+"[odd]",  odd_element_fn ).matrix(basis, p),
        }
    )
