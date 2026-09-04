# zeeman_slower_2 — Phase 4 follow-ups design

**Date:** 2026-05-14
**Branch:** fix-for-distrib
**Target:** `Jupyter Notebooks/RaX/zeeman_slower_2.ipynb`
**Predecessor:** [`docs/superpowers/specs/2026-05-13-zeeman-slower-2-comprehensive-hardening-design.md`](./2026-05-13-zeeman-slower-2-comprehensive-hardening-design.md) (Phases 1–3, complete)

## Motivation

Three follow-ups surfaced at Phase 3 close:

1. **Cell 46 (`s_remix` sweep) was not refactored to the new helper trio.** Its sweep axis (`s_remix`), bespoke physics-rich print table (`s_remix | B_⊥ | Ω_L | F_max | ratio | v_capture`), and overlay-style v(z)/v(t) plot don't match the grid-style helpers' shape. Phase 3b skipped it; the deviation needs a principled resolution and minimal de-duplication.
2. **Regression baseline (cell 54 + 55) has unclear semantics.** The file `zs_baselines_pre_refactor.json` is named as if it were an anchored snapshot but actually behaves as a continuously-overwritten last-run record. There is no audit trail when an intentional change updates the baseline, and the bless workflow (delete file or comment cell 54) is unergonomic. Continuous-detection semantics also allow a buggy state to silently become "the new truth".
3. **Cell 46 retains the `F_max_default = None` lazy-init from before Phase 1 task 1e.** Phase 1 added a comment but didn't refactor; the implicit dependency on `remix_values` being sorted ascending is a footgun.

These three are not independent. Fixes 1 and 3 both touch cell 46 and combine into one coherent edit. Fix 2 is independent.

## Goals

1. **Honest abstraction boundary at cell 46.** Make the helper trio's intentional non-coverage of the `s_remix` axis explicit, not accidental. Eliminate the small bit of duplication that *can* be cleanly factored.
2. **Anchored, audit-trailed regression baseline** matching the established golden-file pattern from pytest-regressions / Jest / Insta. Acceptance must live outside the notebook source code (no flag-in-cell footgun).
3. **Cell 46 reads top-down.** Remove the lazy-init; make the F_max-default anchor explicit before the print loop.

## Non-goals

- Generalizing the helper trio to be axis-agnostic (Brainstorming approach A — rejected: cell 46's duplication is too small to justify the abstraction cost).
- Building parallel `run_remix_sweep` / `format_remix_table` / `plot_remix_overlay` helpers (approach B — rejected: one caller doesn't justify a parallel API).
- Statistical-tolerance or property-based regression checks. Exact-value-with-rtol is the right primitive for this notebook (capture sets are integer-valued; F_max comparisons can absorb FP noise via `rtol=1e-6`).
- Hashing source cells to auto-bless on source change. Brittle (cosmetic edits trigger re-blessing) and pattern analysis shows nobody does this.
- Adding metadata to the baseline JSON. YAGNI; revisit if a `git blame` is ever unclear.

## Three-fix structure

Each fix has its own commit. Fixes are presented below in conceptual order (1, 2, 3); the **Sequencing** section reorders to (1, 3, 2) for execution because Fix 3 consumes Fix 1's helper and Fix 2 is independent.

### Fix 1 — `_captured_at` helper extraction

Add a fourth function to the existing helper-trio cell (cell 38):

```python
def _captured_at(results, key_fn, v0_iter):
    """Sorted list of v0 captured (exit_reason=='v_stop') for a given key_fn.
    Used by format_capture_table and by cell 46's bespoke print table."""
    return sorted([float(v0) for v0 in v0_iter
                   if results[key_fn(v0)]["exit_reason"] == "v_stop"])
```

Refactor `format_capture_table` to call `_captured_at` instead of the inline list comprehension. (Sub-1-line change; eliminates the only piece of cell 46's loop that can be cleanly factored.)

**Cell 46 stays inline for its main loop** — its `s_remix` axis, bespoke physics-rich table, and overlay-style v(z)/v(t) plot don't fit the grid-style helpers' shape. The "Out of scope" section below records this as a deliberate deviation, not an oversight.

### Fix 2 — Env-var bless + git-tracked baseline (G')

Three coordinated changes:

1. **Rename file:** `zs_baselines_pre_refactor.json` → `zs_baselines.json`. Update both cell 54 and cell 55 paths.
2. **Untrack from `.gitignore`:** remove the line `Jupyter Notebooks/RaX/zs_baselines_pre_refactor.json`. Add the renamed file to git.
3. **Cell 54 logic:**
   ```python
   import os as _os
   _bless = bool(_os.environ.get('ZS_BASELINE_BLESS'))
   if _bless:
       print("⚡ ZS_BASELINE_BLESS set — skipping check; cell 55 will write the new baseline.")
   elif not _os.path.exists(_baseline_path):
       print(f"⚠ {_baseline_path} not found — skipping (cell 55 will create it).")
   else:
       # ... existing diff-and-raise logic ...
   ```
   Cell 54's docstring also documents the bless workflow.
4. **Cell 55 unchanged** — always writes (so a blessed run produces the new committable baseline).

**Bless workflow** (mirrors `jest --updateSnapshot`):

```bash
# 1. Make code change
# 2. Run normally → cell 54 FAILS, deltas printed
# 3. Confirm intentional, then bless:
ZS_BASELINE_BLESS=1 conda run -n Structure jupyter nbconvert --to notebook \
    --execute --inplace "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"
# 4. Audit:
git diff "Jupyter Notebooks/RaX/zs_baselines.json"
# 5. Commit:
git add "Jupyter Notebooks/RaX/zs_baselines.json"
git commit -m "physics: <reason>"
```

**Why this beats a cell-source `BASELINE_BLESS = False` flag:** the flag would be source state — forget to flip back and every subsequent run silently writes whatever it produces. Pattern survey of pytest-regressions, Jest, and Insta confirms all three use external (CLI flag / env var) acceptance, not source state. Self-cleaning UX, no footgun.

### Fix 3 — Cell 46 two-phase split

Restructure cell 46's main loop into two phases:

```python
# Phase 1 — precompute all F-grids and F_max values
F_grids_remix = {}
for rmx in remix_values:
    F_grids_remix[rmx] = precompute_F_grid_remix(TM, Bd_remix, v_design,
                                                    sb_template_remix, rmx,
                                                    je_ref, ig_ref, v_grid, B_grid)
F_max_remix   = {rmx: float(np.abs(Fg).max()) for rmx, Fg in F_grids_remix.items()}
F_max_default = F_max_remix[remix_values[0]]   # explicit "ascending-sorted" anchor

# Phase 2 — trajectory sweep + print table
print(...)   # header
results_remix = {}
for rmx in remix_values:
    F_interp_r = RegularGridInterpolator((v_grid, B_grid), F_grids_remix[rmx],
                                            method="linear", bounds_error=False, fill_value=0.0)
    for v0 in v0_remix_scan:
        sol = simulate_trajectory(v0, B_of_z_remix, F_interp_r, mass_RaF, L_remix,
                                    B_grid_lo=float(B_grid[0]),
                                    B_grid_hi=float(B_grid[-1]))
        results_remix[(rmx, v0)] = sol
    captured = _captured_at(results_remix, lambda v0, rmx=rmx: (rmx, v0), v0_remix_scan)
    v_cap = max(captured) if captured else None
    # ... rest of print row (B_⊥, Ω_L, F_max, ratio, v_capture) unchanged ...
```

Same total work, no lazy-init, no implicit ordering dependency. Composes with Fix 1's `_captured_at`.

The original-cell comment block above `F_max_default = None` is also removed — the new code is self-documenting.

## Sequencing

| Step | Commit | Notebook re-execute |
|---|---|---|
| Fix 1 | `_captured_at` extraction in cell 38; refactor `format_capture_table` | clean (no behavior change) |
| Fix 3 | Cell 46 two-phase split using `_captured_at` | clean; cell 46 output identical |
| Fix 2 | Rename file, untrack, cell 54 env-var logic; commit `zs_baselines.json` | clean; regression PASSes (no env var, baseline matches) |

After Fix 2's commit, run a one-off **bless smoke test** (no committed change):

```bash
ZS_BASELINE_BLESS=1 conda run -n Structure jupyter nbconvert --to notebook \
    --execute --inplace "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"
git diff "Jupyter Notebooks/RaX/zs_baselines.json"   # expect empty diff
git checkout "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"   # discard re-exec outputs
```

Empty diff confirms cell 55 is reproducible and the bless path doesn't accidentally drift.

## Verification

Per project CLAUDE.md, notebook-only changes verify via `jupyter nbconvert --to notebook --execute --inplace` (NEVER `--allow-errors`; cf. memory `feedback_no_allow_errors.md`). Each fix's commit must produce `output_type=='error'` count == 0, validated by reading the .ipynb JSON, not by trusting nbconvert's exit code.

After all three fixes:

1. Notebook re-executes clean. Cell 54 PASSes (no env var set).
2. `_captured_at` is defined exactly once and called from both `format_capture_table` and cell 46.
3. Cell 46's main loop reads top-down; no `F_max_default = None` line remains.
4. `zs_baselines.json` exists, is tracked in git, contains current sweep results.
5. `.gitignore` no longer mentions `zs_baselines_pre_refactor.json` or any equivalent.
6. Bless smoke test produces empty `git diff` on the JSON.

## Out of scope (explicit deferrals)

- **`run_remix_sweep` parallel helper family** — rejected. Cell 46's duplication footprint after Fixes 1+3 is just the v0 inner loop with `simulate_trajectory` (~5 lines). A second helper family for one caller is over-abstraction.
- **Axis-agnostic helper redesign** — rejected. Would re-touch the four already-refactored cells (39, 43, 48, 53) without enough payoff. The grid-shape vs single-axis distinction is genuine and worth honoring.
- **Per-sublevel `g_F` in `b_perp_to_s_remix`** — Open Question 1 from the Phase 1–3 spec; still standing.
- **`max_step` integrator-tolerance validation** — Open Question 4 from the Phase 1–3 spec; baseline ran clean, so still deferred.
- **Cross-notebook `simulate_spectra` / `select_dipole` cleanup** — out of scope; separate conversation.
- **`Source Code/` changes** — none required. Notebook-only.

## Pattern provenance

Approach G' (env-var bless + git-tracked baseline) was selected after a structured pattern survey:

| Tool | Storage | Acceptance | Semantic |
|---|---|---|---|
| pytest-regressions | git-tracked golden files | `--force-regen` CLI flag | anchored |
| Jest | `__snapshots__/*.snap` in git | `--updateSnapshot` CLI flag | anchored |
| Insta (Rust) | `*.snap` in git | `cargo insta accept` command | anchored |

All three: file in git + acceptance mechanism external to source. Env var for this notebook (no CLI test runner) is the natural equivalent.
