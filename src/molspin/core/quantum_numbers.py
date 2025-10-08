from __future__ import annotations
from dataclasses import dataclass
from typing import Mapping, Tuple, List, Dict, Any, Iterable, Sequence

QNOrder = Tuple[str, ...]

@dataclass(frozen=True)
class QuantumState:
    """Container for the quantum numbers of a single basis state."""
    values: Mapping[str, Any]
    def as_tuple(self, order: QNOrder) -> Tuple[Any, ...]:
        return tuple(self.values[k] for k in order)

class Basis:
    """Finite, ordered basis defined by a label order and a list of states."""
    def __init__(self, name: str, order: QNOrder, states: Iterable[QuantumState], *, metadata: Dict[str, Any] | None = None):
        self.name = name
        self.order = tuple(order)
        self.states: List[QuantumState] = list(states)
        self.metadata: Dict[str, Any] = dict(metadata or {})

        if not self.order:
            raise ValueError("Basis order must have at least one label")
        if not self.states:
            raise ValueError("Basis must contain at least one state")
        for s in self.states:
            missing = [k for k in self.order if k not in s.values]
            if missing:
                raise ValueError(f"State missing labels {missing}")

        self.index: Dict[Tuple[Any, ...], int] = {s.as_tuple(self.order): i for i, s in enumerate(self.states)}
        if len(self.index) != len(self.states):
            raise ValueError("Duplicate quantum number tuples detected in basis")

    def __len__(self) -> int: return len(self.states)
    def state(self, i: int) -> QuantumState: return self.states[i]
    def idx(self, qn_dict: Mapping[str, Any]) -> int:
        return self.index[tuple(qn_dict[k] for k in self.order)]

    def dict_of_lists(self) -> Dict[str, List[Any]]:
        out: Dict[str, List[Any]] = {k: [] for k in self.order}
        for s in self.states:
            for k in self.order:
                out[k].append(s.values[k])
        return out

def basis_from_dict(name: str, order: Sequence[str], qnumbers: Mapping[str, Sequence[Any]]) -> Basis:
    order = tuple(order)
    lengths = {k: len(v) for k, v in qnumbers.items()}
    if not lengths:
        raise ValueError("Empty qnumbers dictionary")
    n = lengths[next(iter(lengths))]
    if any(L != n for L in lengths.values()):
        raise ValueError("All qnumber lists must have equal length")
    for k in order:
        if k not in qnumbers:
            raise ValueError(f"Missing key '{k}' in qnumbers")
    states = [QuantumState({k: qnumbers[k][i] for k in order}) for i in range(n)]
    return Basis(name=name, order=order, states=states)

def dict_from_basis(basis: Basis) -> Dict[str, List[Any]]:
    return basis.dict_of_lists()

__all__ = ["QuantumState","Basis","basis_from_dict","dict_from_basis","QNOrder"]
