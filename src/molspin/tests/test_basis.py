import pytest
from molspin.core import Basis, QuantumState, basis_from_dict, dict_from_basis

def test_basis_roundtrip_from_legacy():
    qn = {"N":[0,1,1], "S":[0.5,0.5,0.5], "J":[0.5,1.5,0.5], "M":[-0.5,0.5,0.5]}
    order = ("N","S","J","M")
    b = basis_from_dict("caseB_doublet", order, qn)
    assert len(b) == 3
    assert dict_from_basis(b) == {k: qn[k] for k in order}

def test_basis_index_and_state_access():
    states = [
        QuantumState({"N":0,"S":0.5,"J":0.5,"M":-0.5}),
        QuantumState({"N":1,"S":0.5,"J":1.5,"M":0.5}),
    ]
    b = Basis("B", ("N","S","J","M"), states)
    assert b.idx({"N":0,"S":0.5,"J":0.5,"M":-0.5}) == 0
    assert b.idx({"N":1,"S":0.5,"J":1.5,"M":0.5}) == 1
    assert b.state(1).values["J"] == 1.5

def test_duplicate_states_raises():
    s = QuantumState({"N":0,"S":0.5,"J":0.5,"M":0.5})
    with pytest.raises(ValueError):
        Basis("dup", ("N","S","J","M"), [s, s])
