# CLAUDE.md — Molecule-Structure (this fork)

Personal fork (`arianjad/Molecule-Structure`). The codebase at `~/Documents/Molecular-Structure/` is a historical lab fork (HutzlerLab) — separate codebase, not authoritative here. Don't import from it or cross-reference unprompted.

## Working copy & Drive snapshot

Git working copy (day-to-day): **`~/Code/Molecule-Structure`** — clone of `arianjad/Molecule-Structure`; GitHub is the source of truth. Do all edits and git here.

The Google Drive path (`…/EMA-data-server/RaX/Personal/ArianJadbabaie/Code/Molecule-Structure`) is a **git-less snapshot** of known-good code, decoupled 2026-05-18 — no `.git/`. Never run git there: it falls through to the `~/.git` home-config repo (which ignores the Drive tree, so it's inert, but it's still the wrong repo). Refresh the snapshot from the working copy with `scripts/refresh-drive-snapshot.sh` (dry-run by default; `--go` to apply). The refresh is **additive — it never deletes anything on Drive**: tracked code is overwritten with the working-copy version, new files added, and all gitignored Drive-only content (`Old Code/`, `Figures/`, RaX `figures/` & `reports/`, baselines, `*.zip`, `thinking/`) is left untouched. `Data/` and `From Others/` are gitignored but carried into the working copy too. Stale files removed/renamed in the repo linger on the snapshot — prune by hand if a clean mirror is ever wanted (no automated `--delete`: irreversible on Drive).

## Working frame

Define a molecule Hamiltonian → diagonalize → explore properties. Spectra, Stark/Zeeman maps, level diagrams, and PT-violating shifts are all downstream.

## Repo map

- `Source Code/` — the library. Flat module directory (not a package). Active edits land here.
  Public API: `Energy_Levels.MoleculeLevels`. Dispatch: `molecule_library_class.Molecule_Library`. Constants: `molecule_parameters.py`. Matrix elements per Hund's case: `matrix_elements.py`.
- `Jupyter Notebooks/` — analyses. Has its own CLAUDE.md governing that zone.
- `Old Code/`, `From Others/`, `CaOH_Supplemental/`, `*.zip` — reference / archival material. Useful to read; not authoritative — check `Source Code/` before treating any pattern as canonical.
- `Data/` — large reference data. Don't read into context unless the task needs it.
- `Figures/` — figure outputs (gitignored).

## Branches

Default `main`. `fix-for-distrib` is the current working branch. `BaOH`, `YbOH`, `dev`, `molspin` are archival.

## Environment

- `Structure` — for running this codebase (notebooks, source-code imports). Use `conda run -n Structure ...`.
- `claude-code` — your sandbox for unrelated utility commands. Not for running this codebase.

## Verification gate

The single completeness check for any `Source Code/` change: the tutorial notebook runs end-to-end on a fresh kernel.

```
conda run -n Structure jupyter execute "Jupyter Notebooks/RaF_Calcs_Tutorial.ipynb"
```

Module edits leave the user's running kernel stale — fresh-kernel re-execution is the verification, not autoreload.

## Source Code rules

**`molecule_parameters.py` integrity.** When adding or changing a constant: inline-comment the source (paper / NIST / theory note), name the isotope, confirm units. Most entries are MHz; constants given in cm⁻¹ are multiplied by `c` inline. Don't propagate a value across isotopes/states without flagging it.

**YbOH-backend wart.** `Molecule_Library` routes everything except `CaOH` through the YbOH backend (`molecule_library_class.py:40-43`). Known structural wart; extending it to new molecules is fine, but don't refactor unprompted.

**Physics conventions.** See `README.md` (units MHz/V·cm⁻¹/Gauss; polyatomic K/P notation; Hund's cases; PGopher-flipped sign conventions; spin-statistics flag). Don't fabricate or guess on coupling cases, matrix elements, or selection rules — stop and ask.

**Library API currency.** Before nontrivial SymPy, QuTiP, or SciPy-linalg work, consult the context7 MCP (`mcp__context7__*`) for the installed version's API — these are version-volatile. Not for stable NumPy/stdlib.

## References

Brown & Carrington, *Rotational Spectroscopy of Diatomic Molecules* is the canonical reference for diatomic structure (matrix elements, sign conventions, basis transformations, coupling cases). Two copies, same edition:

- Local (prefer): `find ~/Zotero/storage -maxdepth 2 -iname 'Brown and Carrington*.pdf'` — hash subdir varies per machine.
- Cross-machine: `/Users/arianjadbabaie/Library/CloudStorage/GoogleDrive-arianjad@mit.edu/Shared drives/EMA-data-server/Books/AMO/Brown and Carrington, Rotational Spectroscopy of Diatomic Molecules.pdf` — may be a Google Drive cloud-only placeholder; `pdf_info` returns `Failed to open file` until hydrated.

Consult via `mcp__pdf-mcp__pdf_search` / `pdf_read_pages` before answering diatomic-physics questions. The `brown-carrington` skill provides the chapter map, resolution order, and workflow.
