"""
molspin.core — Core data structures and operations for molecular Hamiltonian modeling.

This subpackage defines the low-level building blocks:
    - Quantum numbers and basis representations
    - Basis registries and converters
    - Operator abstractions and registries
    - Effective Hamiltonians and state models
"""

# Quantum numbers and bases
from .quantum_numbers import QuantumState, Basis, basis_from_dict, dict_from_basis, QNOrder
from .basis_registry import (
    register_basis,
    get_basis,
    register_converter,
    converter,
    has_converter,
    unitary,
)

# Operators and Hamiltonians
from .operators import Operator, RegisteredOperator, SiteAwareOperator
from .operator_registry import register_operator, get_operator
from .state_models import EffectiveHamiltonian

__all__ = [
    # Quantum numbers and bases
    "QuantumState",
    "Basis",
    "basis_from_dict",
    "dict_from_basis",
    "QNOrder",
    "register_basis",
    "get_basis",
    "register_converter",
    "converter",
    "has_converter",
    "unitary",

    # Operators and Hamiltonians
    "Operator",
    "RegisteredOperator",
    "SiteAwareOperator",
    "register_operator",
    "get_operator",
    "EffectiveHamiltonian",
]

# Optional version note
__version__ = "0.1.0"
