from __future__ import annotations
from typing import Dict
from .operators import Operator

_OPS: Dict[str, Operator] = {}

def register_operator(op: Operator) -> None:
    if op.name in _OPS:
        raise ValueError(f"Operator '{op.name}' already registered")
    _OPS[op.name] = op

def get_operator(name: str) -> Operator:
    try:
        return _OPS[name]
    except KeyError as e:
        raise KeyError(f"Unknown operator '{name}'. Registered: {list(_OPS)}") from e

__all__ = ["register_operator","get_operator"]
