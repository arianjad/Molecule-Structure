from __future__ import annotations
from dataclasses import dataclass
from typing import Mapping, Tuple, List, Dict, Any, Iterable, Sequence

QNOrder = Tuple[str, ...]  # list/tuple of quantum-number label strings, e.g. ("N","S","J","M")

@dataclass(frozen=True)
class QuantumState:
    """Container for the quantum numbers of a single basis state."""
    values: Mapping[str, Any]
    def as_tuple(self, qstr_list: QNOrder) -> Tuple[Any, ...]:
        return tuple(self.values[k] for k in qstr_list)

class Basis:
    """
    Finite, ordered basis defined by a list of quantum-number labels and a list of states.

    Attributes
    ----------
    name : str
        Human-readable basis name.
    qstr_list : tuple[str, ...]
        The ordered quantum-number labels used for indexing (e.g. ("N","S","J","M")).
    states : list[QuantumState]
        The basis vectors (each is an immutable mapping of label -> value).
    index : dict[tuple, int]
        Map from the tuple of q.n. values (following qstr_list) to the integer index.
    qnumbers : dict[str, list[Any]]
        Legacy dict-of-lists view of the basis quantum numbers (cached at init).
    metadata : dict[str, Any]
        Extra info (e.g., {"nuclear_spins": {"M": I_M, "H": I_H}}).
    """
    def __init__(
        self,
        name: str,
        qstr_list: QNOrder,
        states: Iterable[QuantumState],
        *,
        metadata: Dict[str, Any] | None = None
    ):
        self.name = name
        self.qstr_list = tuple(qstr_list)
        self.states: List[QuantumState] = list(states)
        self.metadata: Dict[str, Any] = dict(metadata or {})

        if not self.qstr_list:
            raise ValueError("Basis must define at least one quantum-number label")
        if not self.states:
            raise ValueError("Basis must contain at least one state")

        # Validate each state has the required labels
        for s in self.states:
            missing = [k for k in self.qstr_list if k not in s.values]
            if missing:
                raise ValueError(f"State missing labels {missing}")

        # Build index: q-number tuple -> integer position
        self.index: Dict[Tuple[Any, ...], int] = {
            s.as_tuple(self.qstr_list): i for i, s in enumerate(self.states)
        }
        if len(self.index) != len(self.states):
            raise ValueError("Duplicate quantum number tuples detected in basis")

        # Build & cache the legacy dict-of-lists representation
        qn: Dict[str, List[Any]] = {k: [] for k in self.qstr_list}
        for s in self.states:
            for k in self.qstr_list:
                qn[k].append(s.values[k])
        self.qnumbers: Dict[str, List[Any]] = qn

    # --- Common conveniences ---
    def __len__(self) -> int:
        return len(self.states)

    def state(self, i: int) -> QuantumState:
        return self.states[i]

    def idx(self, qn_dict: Mapping[str, Any]) -> int:
        key = tuple(qn_dict[k] for k in self.qstr_list)
        return self.index[key]

# ---------- Adapters to/from legacy formats ----------

def basis_from_dict(name: str, qstr_list: Sequence[str], qnumbers: Mapping[str, Sequence[Any]]) -> Basis:
    """
    Create a Basis from the legacy dict-of-lists structure.
    Example qnumbers: {"N":[...], "S":[...], "J":[...], "M":[...]} (all lists of equal length).
    """
    qstr_list = tuple(qstr_list)
    lengths = {k: len(v) for k, v in qnumbers.items()}
    if not lengths:
        raise ValueError("Empty qnumbers dictionary")
    n = lengths[next(iter(lengths))]
    if any(L != n for L in lengths.values()):
        raise ValueError("All qnumber lists must have equal length")
    for k in qstr_list:
        if k not in qnumbers:
            raise ValueError(f"Missing key '{k}' in qnumbers")

    states = [QuantumState({k: qnumbers[k][i] for k in qstr_list}) for i in range(n)]
    return Basis(name=name, qstr_list=qstr_list, states=states)

def dict_from_basis(basis: Basis) -> Dict[str, List[Any]]:
    """Return the cached dict-of-lists view."""
    return basis.qnumbers

__all__ = ["QuantumState", "Basis", "basis_from_dict", "dict_from_basis", "QNOrder"]
