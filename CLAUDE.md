# CLAUDE.md — Molecule-Structure (this fork)

Personal fork (`arianjad/Molecule-Structure`).

## Working copy & Drive snapshot

Git working copy (day-to-day): **`~/Code/Molecule-Structure`** — clone of `arianjad/Molecule-Structure`; GitHub is the source of truth. Do all edits and git here.

The Google Drive path (`…/EMA-data-server/RaX/Personal/ArianJadbabaie/Code/Molecule-Structure`) is a **git-less snapshot** — no `.git/`. Never run git there: it falls through to the `~/.git` home-config repo (inert, but the wrong repo). Refresh it with `scripts/refresh-drive-snapshot.sh` (dry-run by default, `--go` to apply; additive only, never deletes on Drive — the script header has the details). `Data/` and `From Others/` are gitignored but carried into the working copy too.

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

- `structure` — for running this codebase (notebooks, source-code imports). Use `conda run -n structure ...`.
- `claude-code` — your sandbox for unrelated utility commands. Not for running this codebase.

## Verification gate

The single completeness check for any `Source Code/` change: the tutorial notebook runs end-to-end on a fresh kernel.

```
conda run -n structure jupyter execute "Jupyter Notebooks/RaF_Calcs_Tutorial.ipynb"
```

Module edits leave the user's running kernel stale — fresh-kernel re-execution is the verification, not autoreload.

**Windows encoding dep:** `jupyter execute` (nbclient) re-saves the notebook on completion; without `PYTHONUTF8=1` user-wide it reads under cp1252 and writes UTF-8, silently baking mojibake into source. Setup and repair recipe: `repl-platform` skill (PYTHONUTF8 bullet).

## Source Code rules

**`molecule_parameters.py` integrity.** When adding or changing a constant: inline-comment the source (paper / NIST / theory note), name the isotope, confirm units. Most entries are MHz; constants given in cm⁻¹ are multiplied by `c` inline. Don't propagate a value across isotopes/states without flagging it.

**YbOH-backend wart.** `Molecule_Library` routes everything except `CaOH` through the YbOH backend (`molecule_library_class.py:40-43`). Known structural wart; extending it to new molecules is fine, but don't refactor unprompted.

**Physics conventions.** See `README.md` (units MHz/V·cm⁻¹/Gauss; polyatomic K/P notation; Hund's cases; PGopher-flipped sign conventions; spin-statistics flag). Don't fabricate or guess on coupling cases, matrix elements, or selection rules — stop and ask.

**Library API currency.** Before nontrivial SymPy, QuTiP, or SciPy-linalg work, consult the context7 MCP (`mcp__context7__*`) for the installed version's API — these are version-volatile. Not for stable NumPy/stdlib.

## RaF-specific notes (per Arian, 2026-05-20)

- `molecule_parameters.py:58` `RaF.boson.A0 'ASO': 1350` is **stale**; canonical = 2067.6 cm⁻¹ (matches `~/Downloads/225_RaF.pgo` Excited `A`). The codebase value works as a held-fixed input but isn't the right spin-orbit constant.
- `molecule_parameters.py:72` `'Origin': 13284.427` cm⁻¹ is **T_Π1/2** (Udrescu Tbl I fit). PGopher convention is **T_e** (electronic origin); they differ by ½·A_SO. `225_RaF.pgo` Excited `Origin=14318.33` is T_e — same physical ²Π₁/₂ origin, different convention. Reference: `~/obsidian-vault/04-resources/references/raf-effective-hamiltonian/source.md` Appendix D.
- `'bF':96.3`, `'c':19` (X0); `'a':19/2` (A0); `'muE':3.91 D` (X0), `1.0 D` (A0) are **Skripnikov ab initio**, personal/unpublished communication — there is no paper to cite.
- `~/Downloads/225_RaF.pgo` is the **²²⁵RaF (fermion) file**; its single `<LinearNucleus>` block is **Ra hyperfine, not ¹⁹F** (magnitudes ~16 GHz vs codebase F-hyperfine ~100 MHz; signs opposite). Apples-to-oranges trap if mistaken for an ¹⁹F cross-check.

## References

Brown & Carrington, *Rotational Spectroscopy of Diatomic Molecules* is the canonical reference for diatomic structure (matrix elements, sign conventions, basis transformations, coupling cases). Consult it via `mcp__pdf-mcp__pdf_search` / `pdf_read_pages` before answering diatomic-physics questions; the `brown-carrington` skill holds the PDF locations, chapter map, resolution order, and workflow.
