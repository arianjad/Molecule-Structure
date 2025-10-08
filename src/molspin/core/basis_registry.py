from __future__ import annotations
from typing import Callable, Dict, Tuple
import numpy as np
from .quantum_numbers import Basis

_BASIS: Dict[str, Basis] = {}
_CONVERTER: Dict[Tuple[str, str], Callable[[], np.ndarray]] = {}

def register_basis(basis: Basis) -> None:
    if basis.name in _BASIS:
        raise ValueError(f"Basis '{basis.name}' already registered")
    _BASIS[basis.name] = basis

def get_basis(name: str) -> Basis:
    try:
        return _BASIS[name]
    except KeyError as e:
        raise KeyError(f"Unknown basis '{name}'. Registered: {list(_BASIS)}") from e

def register_converter(src: str, dst: str, builder: Callable[[], np.ndarray]) -> None:
    key = (src, dst)
    if key in _CONVERTER:
        raise ValueError(f"Converter {src}->{dst} already registered")
    _CONVERTER[key] = builder

def has_converter(src: str, dst: str) -> bool:
    return (src, dst) in _CONVERTER or src == dst

def unitary(matrix: np.ndarray, tol: float = 1e-10) -> bool:
    if matrix.shape[0] != matrix.shape[1]: return False
    I = np.eye(matrix.shape[0])
    return np.allclose(matrix.conj().T @ matrix, I, atol=tol) and np.allclose(matrix @ matrix.conj().T, I, atol=tol)

def converter(src: str, dst: str) -> np.ndarray:
    if src == dst:
        n = len(get_basis(src))
        return np.eye(n)
    try:
        U = _CONVERTER[(src, dst)]()
    except KeyError as e:
        raise KeyError(f"No converter registered for {src}->{dst}") from e
    if not unitary(U):
        raise ValueError(f"Converter {src}->{dst} is not unitary within tolerance")
    return U

__all__ = ["register_basis","get_basis","register_converter","converter","has_converter","unitary"]
