import numpy as np
from molspin.core import (
    QuantumState, Basis,
    RegisteredOperator, SiteAwareOperator,
    register_operator, get_operator,
    EffectiveHamiltonian,
)

def make_demo_basis(IM: float):
    states = [
        QuantumState({"N":0,"S":0.5,"J":0.5,"M":-0.5}),
        QuantumState({"N":0,"S":0.5,"J":0.5,"M": 0.5}),
    ]
    meta = {"nuclear_spins":{"M": IM, "H": 0.5}}
    return Basis("demo", ("N","S","J","M"), states, metadata=meta)

def test_registered_operator_and_model_sum():
    def diag_eval(basis, params):
        v = np.asarray(params.get("diag", [0.0]*len(basis)))
        return np.diag(v)
    op = RegisteredOperator("DiagOp", diag_eval)
    register_operator(op)

    b = make_demo_basis(0.0)
    terms = [("DiagOp", {"diag":[1.0, 2.0]})]
    H = EffectiveHamiltonian(b, terms).matrix()
    assert np.allclose(H, np.diag([1.0, 2.0]))

def test_siteaware_dispatch_even_vs_odd():
    calls = {"even":0, "odd":0}
    def eval_even(basis, params):
        calls["even"] += 1
        return np.eye(len(basis))
    def eval_odd(basis, params):
        calls["odd"] += 1
        return 2*np.eye(len(basis))

    op = SiteAwareOperator("SpinRotation_like", {"even": eval_even, "odd": eval_odd})
    register_operator(op)

    H_even = get_operator("SpinRotation_like").matrix(make_demo_basis(0.0), {})
    H_odd  = get_operator("SpinRotation_like").matrix(make_demo_basis(0.5), {})
    assert np.allclose(H_even, np.eye(2))
    assert np.allclose(H_odd, 2*np.eye(2))
    assert calls["even"] == 1 and calls["odd"] == 1
