# CLAUDE.md — Jupyter Notebooks zone

Local rules for analysis notebooks. The root CLAUDE.md still applies.

## Where work happens

Most active: `RaX/`, `RaF_Calcs_Tutorial.ipynb`, `Generic/`. New analyses default to `RaX/`.

Older / archival: `YbOH/`, `CaOH/`, `BaOH/`, `YbF/`, `DyO/`. Prior projects; occasionally edited but largely frozen. Two cautions when working there:

1. Patterns may lag the modern API. The current entry point is `from Energy_Levels import MoleculeLevels`. Anything importing `YbOH_Energy_Levels_symbolic` (e.g., `DyO/173YbOH.ipynb`) is pre-refactor — don't propagate that import back into active code.
2. Don't pattern-match those notebooks as canonical without checking the tutorial first.

## Reference

`RaF_Calcs_Tutorial.ipynb` is the canonical API reference. Read it before authoring new notebooks or extending the API surface.

## Boilerplate

First cell of every notebook:

```python
from config_path import add_to_sys_path
add_to_sys_path()
```

Each subdir has its own `config_path.py` (identical files). Keep the import; don't try to package-ify it.

## Helper drift

Functions like `simulate_spectra`, `select_dipole`, `lorentzian` are copy-pasted across active notebooks rather than imported from a shared module. If you fix one, check the others before assuming it's a one-off.

## Figure naming

Convention in `RaX/`: `{isotope}{molecule}_{observable}_{condition}_N{...}_E{...}_v{n}.{pdf,svg}` — e.g. `225RaF_BSense_ZeroG_N123_E5k_v2.pdf`. Emit both PDF and SVG.

## Verification

For library changes (`Source Code/`): root CLAUDE.md verification gate applies — fresh-kernel run of the tutorial via `jupyter execute`.

For notebook-only changes: "Restart & Run All" on the affected notebook (downstream cells depend on the ones you edited; partial re-runs hide stale state). When in doubt about *what* to verify physically, ask (per `~/.claude/rules/physics-verification.md`).
