"""
molspin — Modular molecular structure and effective Hamiltonian toolkit.

This package provides a flexible, physics-friendly API for building and analyzing
effective Hamiltonians of diatomic and linear triatomic molecules.

Typical use:
    from molspin import Basis, QuantumState, EffectiveHamiltonian

Submodules:
    molspin.core      — main classes (Basis, Operators, Hamiltonians)
    molspin.molecules — molecule-specific parameter loaders (coming soon)
    molspin.io        — notebook/plotting utilities (coming soon)
"""

# Expose key API elements directly for convenience
from .core import (
    Basis,
    QuantumState,
    EffectiveHamiltonian,
    RegisteredOperator,
    SiteAwareOperator,
    register_operator,
    get_operator,
)

# Keep access to the whole core namespace for power users
from . import core

__all__ = [
    "Basis",
    "QuantumState",
    "EffectiveHamiltonian",
    "RegisteredOperator",
    "SiteAwareOperator",
    "register_operator",
    "get_operator",
    "core",
]

# Optional: simple version tag
__version__ = "0.1.0"
