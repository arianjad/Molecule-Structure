# Molecule-Structure

Python tools for computing energy levels, eigenvectors, and field shifts (Stark, Zeeman, AC-Stark, PT-violating) of paramagnetic open-shell molecules using effective Hamiltonians.

The library is built around a single class, `MoleculeLevels`, which encapsulates a vibronic state: its quantum-number basis, symbolic and numeric Hamiltonians, eigensystem, parities, and analysis methods.

## Setup

Use the `Structure` conda environment. Dependencies: `numpy`, `sympy`, `matplotlib`, `seaborn`, `torch`, `IPython`, `tabulate`.

The `Source Code/` directory is a flat module directory, not a Python package. Notebooks inject it onto `sys.path` via `config_path.py`. If imports fail, paste the following into the first cell:

```python
import os, sys

def add_to_sys_path(target_dir_name="Source Code"):
    current_dir = os.path.abspath('')
    while True:
        if target_dir_name in os.listdir(current_dir):
            path = os.path.join(current_dir, target_dir_name)
            if path not in sys.path:
                sys.path.append(path)
            break
        parent = os.path.dirname(current_dir)
        if current_dir == parent:
            break
        current_dir = parent

add_to_sys_path()
```

## Usage

The canonical reference is `Jupyter Notebooks/RaF_Calcs_Tutorial.ipynb`. Read it before writing new analysis code.

Minimal example:

```python
from Energy_Levels import MoleculeLevels

X0 = MoleculeLevels.initialize_state(
    molecule_name = 'RaF',
    elec_state    = 'X',
    vib_state     = 0,
    N_list        = [0, 1, 2, 3],
    fermion_or_boson = 'boson',
    I_nuclei      = [0, 1/2],
    isotope       = 226,
    P_values      = [1/2],
)

H = X0.H_function(0, 1e-3)            # numeric Hamiltonian, E=0 V/cm, B=1e-3 G
evals, evecs = X0.eigensystem(0, 0)
X0.StarkMap(Ez_array, plot=True)
```

## Conventions

### Units

All Hamiltonians and eigenvalues are in **MHz**. External fields: **V·cm⁻¹** for `E`, **Gauss** for `B`. Where a constant in `molecule_parameters.py` is sourced in cm⁻¹, it is multiplied inline by `c` (in cm/μs, defined in `params_general`) to convert to MHz.

### Notation

The code uses polyatomic conventions even for diatomic molecules:

- `K` is the projection of `N` on the molecular axis (the diatomic `Λ`).
- `P` is the projection of `J` on the molecular axis (the diatomic `Ω`).

Quantum-number arrays exposed via `MoleculeLevels.q_numbers` use this notation.

### Hund's cases

Three cases are implemented: `bBJ`, `bBS` (used when the metal nucleus has `I ≥ 1`, e.g. ¹⁷¹Yb, ²²⁵Ra), and `aBJ`. The case is selected automatically per `(molecule, isotope, electronic state)` via the `Molecule_Library.cases` dispatch table. Matrix elements live in `matrix_elements.py` keyed by case.

### Sign conventions

This code uses the `B+C` form for several effective-Hamiltonian terms. PGopher uses the opposite sign for `d`-type interactions, so values from PGopher-convention sources are sign-flipped at entry; see inline comments on `'d'`, `'dYb'`, etc. in `molecule_parameters.py`. Likewise some rotational constants are stored in the `R²` form rather than `R`; this is annotated where it matters.

Do not "fix" a minus sign in `molecule_parameters.py` without checking the inline source comment.

### Spin statistics

`fermion_or_boson` at `initialize_state` selects whether the total `F` is half-integer or integer. The flag controls which backend Hamiltonian builder runs (`H_even_*` vs `H_odd_*`) and which `q_numbers_*` builder produces the basis.

## Source Code map

| File | Purpose |
|---|---|
| `Energy_Levels.py` | `MoleculeLevels` class — the public API. |
| `molecule_library_class.py` | `Molecule_Library` — dispatch table mapping `(molecule, isotope, state)` to matrix elements, builders, q-numbers. |
| `molecule_parameters.py` | Effective-Hamiltonian constants, nested by molecule → spin-statistics → vibronic state. Inline-cited where possible. |
| `matrix_elements.py` | Per-Hund's-case matrix-element functions (numeric). |
| `matrix_elements_sym.py` | Symbolic variants. |
| `case_bBJ_MEs.py` | Additional `bBJ` matrix elements. |
| `hamiltonian_builders.py` | Assemble `H` from matrix elements + parameters; PT-violating and TDM builders. |
| `quantum_numbers.py` | Basis-vector dictionaries for each Hund's case. |
| `CPT_Sims.py` | CPT-related simulation helpers. |

## Known wart

`Molecule_Library` routes every molecule except `CaOH` through the YbOH backend (`molecule_library_class.py:40-43`). This is legacy from when the library targeted YbOH only. It is functional and extending it to new molecules is straightforward, but the dispatch is not symmetric.

## Notebook layout

Active analysis lives in `Jupyter Notebooks/RaX/`, with `RaF_Calcs_Tutorial.ipynb` as the reference. The `YbOH/`, `CaOH/`, `BaOH/`, `YbF/`, and `DyO/` subdirectories are older/archival projects; some use the pre-refactor API (`YbOH_Energy_Levels_symbolic`).

## Status

Personal fork (`arianjad/Molecule-Structure`). The codebase at `~/Documents/Molecular-Structure/` is a separate historical lab fork (HutzlerLab) and is not authoritative for this repository.
