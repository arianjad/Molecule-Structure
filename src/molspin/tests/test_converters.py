import numpy as np
from molspin.core import QuantumState, Basis, register_basis, register_converter, converter, unitary

def setup_demo_bases():
    b1 = Basis("B1", ("N","M"), [QuantumState({"N":0,"M":-0.5}), QuantumState({"N":0,"M":0.5})])
    b2 = Basis("B2", ("N","M"), [QuantumState({"N":0,"M":-0.5}), QuantumState({"N":0,"M":0.5})])
    register_basis(b1); register_basis(b2)

def test_register_and_use_converter():
    setup_demo_bases()
    U = (1/np.sqrt(2))*np.array([[1.0,1.0],[1.0,-1.0]])
    register_converter("B1","B2", lambda: U)
    U12 = converter("B1","B2")
    assert unitary(U12)
    psi_B1 = np.array([1.0, 0.0])
    psi_B2 = U12 @ psi_B1
    assert np.allclose(psi_B2, (1/np.sqrt(2))*np.array([1.0,1.0]))
