# zeeman_slower_2 — comprehensive hardening design

**Date:** 2026-05-13
**Branch:** fix-for-distrib
**Target:** `Jupyter Notebooks/RaX/zeeman_slower_2.ipynb`
**Motivation:** Consolidate the cleanup actions surfaced by two code reviews (`649cd36` and `9d675df..7dca1da`) plus the physics investigation of `remix`, ν_eg(B) monotonicity, and the PB σ⁻ F_max-at-low-B artifact, into a single staged change that hardens the analysis without changing physics conclusions.

## Goals

1. **Physics interpretability.** `remix` is currently named as if it were a fraction but enters the rate equation as a saturation parameter (s) of a fictitious always-on-resonance isotropic pseudo-laser. Rename + recalibrate so a future reader can map remix sweeps to physical off-axis B-field magnitude.
2. **Robustness.** Add guards on the silent failure modes the reviews surfaced: non-monotonic ν_eg(B) inversion, multi-N basis breaking the energy-sort m_S identification, hyperfine-group clustering with arbitrary tolerance.
3. **De-duplication.** Five trajectory sweeps share ~30 lines of identical boilerplate. Extract into helpers so a future fix is a one-place change.
4. **No physics regression.** Capture-velocity tables and per-config F_max must be bit-identical before and after the refactor (or differences must be explainable in terms of explicit physics changes, not coincidence).

## Non-goals

- Adding pytest infrastructure (the project's verification gate is full-notebook re-execution).
- Backward compatibility with sibling notebooks — the renamed `s_remix` parameter and helper function signatures are internal.
- Implementing Commit B (N=3 leakage / repumper) or the full Liang-Yan 12-sideband scheme. Those are downstream physics work, gated on this cleanup.
- Adding stochastic recoil heating to the trajectory integrator.

## Three-phase structure

The work is broken into three sequential commits. Each commit's verification gate is
`conda run -n Structure jupyter execute --inplace "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"`
exiting 0 with zero `output_type=='error'` cells (validated by reading the .ipynb JSON, NOT by `--allow-errors`).

### Phase 1 — Documentation + cheap fixes (warmup)

Establishes the regression baseline before bigger changes. No behavior change anywhere.

| # | Action | File / cell |
|---|---|---|
| 1a | Upgrade `steady_state_rate_eq` docstring: spell out that `remix` is a saturation parameter (not a fraction); cross-reference the Larmor mapping introduced in Phase 2b. | cell 27 |
| 1b | Add monotonicity check at the top of `make_B_of_z_nonlinear`: if `np.diff(nu_eg_curve)` changes sign, print a warning naming the (je, ig) pair and the B range, and proceed with the current `argsort + interp1d` behavior. Does not halt. | cell 33 |
| 1c | Add `Total captures: X/Y` summary line to the capture-velocity tables in cells 36 and 40, matching the existing summary in cells 44 and 48. | cells 36, 40 |
| 1d | After `find_X_N1_hyperfine_groups(g)` in cell 38, print `⚠ ...` warning if `len(groups) != 4` (RaF X²Σ⁺ N=1 should give four (J, F) groups: J=½ F=0, J=½ F=1, J=3/2 F=1, J=3/2 F=2). Continue execution either way. | cell 38 |
| 1e | Replace the `F_max_default = None ... if F_max_default is None:` lazy-init pattern in cell 42 with `F_max_default = F_max_remix[remix_values[0]]` computed after the loop, OR with explicit `F_max_default = F_max_remix[min(remix_values)]` if iteration order isn't guaranteed. | cell 42 |
| 1f | New cell at the end of the notebook: serialize all sweep results (capture-velocity range tables and per-config F_max values) to `Jupyter Notebooks/RaX/zs_baselines_pre_refactor.json`. Cell prints a one-line summary of what was saved. Used by Phase 3e for regression check. | new last cell |

**Commit message:** `zeeman_slower_2: Phase 1 — documentation upgrades + monotonicity warning + regression baseline`

### Phase 2 — Physics diagnostic cells

Additive only — no existing code is modified except the `remix` → `s_remix` rename (mechanical) and cell 42's secondary axis.

| # | Action | New cell location |
|---|---|---|
| 2a | New cell after cell 33 (the `je_ref, ig_ref` pick): plot `ν_eg(B)` over `Bz = [1e-3, 1000]` G for the chosen `(je_ref, ig_ref)`. Identify local extrema by sign changes of `np.diff(nu_eg_curve)`. Print "MONOTONIC" or "NON-MONOTONIC at B ≈ [...] G". | after cell 33 |
| 2b | Rename `remix` → `s_remix` in `steady_state_rate_eq` signature, `precompute_F_grid_remix` signature, and all call sites (cells 27, 42, and any others surfaced by grep). Add helper `b_perp_to_s_remix(B_perp_G, Gamma_FWHM_MHz, g_F=0.5)` returning `s_remix` such that the pseudo-laser pumping rate equals the Larmor mixing rate `Ω_L = g_F μ_B B_⊥ / ℏ`. Derivation:  `s_remix · (Γ_angular/2) = Ω_L`  ⟹  `s_remix = 2 g_F μ_B B_⊥ / (ℏ Γ_angular) = 2 g_F μ_B B_⊥ / (h Γ_FWHM)` (cyclic units). New markdown cell before cell 42 with the derivation and a numeric mapping table (e.g. B_⊥ = 0.1, 1, 3.3, 10, 30 G ↔ s_remix). | helper cell after cell 27; markdown before cell 42 |
| 2c | New cell after cell 48 (the PB F-grid): "comb selectivity" diagnostic. For B ∈ [50, 1000] G at v_design = 80 m/s, decompose the total scattering rate from `steady_state_rate_eq` (via the existing `R_per_pol` return) into "intended" (transitions terminating in the m_S=+½ manifold per `g_char_full['mS']`) vs "off-target". Plot intended-fraction vs B. Quantifies the F_max-at-low-B artifact we observed at s_sb=100. | after cell 48 |
| 2d | Modify cell 42 (existing remix sweep): add a second top x-axis showing B_⊥ in Gauss derived from `b_perp_to_s_remix`. Print the physical Larmor frequency Ω_L (in MHz cyclic) alongside each `s_remix` value in the capture table. | cell 42 |
| 2e | Markdown updates: a short paragraph in each Step (6, 7, 7b, 7c, PB) referencing what the new diagnostics show. Specifically: Step 7b notes the s_remix↔B_⊥ mapping and the captured-trajectory threshold corresponds to B_⊥ ≈ 3 G; PB section notes the off-target fraction explains the s=100 anomaly. | existing markdown cells |

**Commit message:** `zeeman_slower_2: Phase 2 — physics diagnostics, s_remix rename + B_perp calibration, PB comb selectivity`

### Phase 3 — Refactor

Highest-risk phase. Regression-gated against the Phase 1f baseline.

| # | Action | File / cell |
|---|---|---|
| 3a | New "sweep helpers" cell after cell 35 (`simulate_trajectory`) defining three functions: `run_trajectory_sweep(F_interps_dict, configs, je_ref_lookup, ig_ref_lookup, B_grid, v_i_sweep, v_f_sweep, v_design, v0_scan, a_max_est, mass)` → flat dict keyed by `(B_design, s_per_sb, L, v0)`; `format_capture_table(results, configs)` → str; `plot_v_of_z_grid(results, L_show, configs, v0_cmap, v0_norm, suptitle)` → fig. Configs argument is a list of `(B_design, s_per_sb, L_values)` tuples (or a single tuple) — handles both Step-6/7/7c (3×B × 2×s × 4×L) and Step-7b/PB (1×B × 1×s × 4×L) shapes uniformly. | new cell after cell 35 |
| 3b | Replace duplicated trajectory-sweep + table + plot boilerplate in cells 36, 40, 42, 44, 48 with thin wrappers around the helpers. Each cell shrinks to ~10 lines: prepare `F_interps_dict` (already done in earlier cells) + `je_ref / ig_ref` mapping → call `run_trajectory_sweep` → call `format_capture_table` → call `plot_v_of_z_grid`. Keep all step-specific markdown headers and narrative comments. | cells 36, 40, 42, 44, 48 |
| 3c | Switch `build_paschen_back_sigma_minus` (cell 46) manifold selection from energy-sort to `g_char_full['mS']`: `ig_list = sorted(int(i) for i in np.where(g_char_full['mS'][iB_design] > 0)[0])` for upper manifold, `< 0` for lower. Assert `len(ig_list) == n_sidebands` with an informative message if it fails (e.g. multi-N basis). | cell 46 |
| 3d | In cell 47 (PB F-grid), build `sb_pb_template`, `je_pb_ref`, `ig_pb_ref` ONCE (already done at cell-46 print step). For each `s_per_sb` value, construct the per-call sideband list by varying only the `s` field of each entry in the template. Drop the `ref_pb[s_per_sb]` dict and pass the single `(je_pb_ref, ig_pb_ref)` pair to both `precompute_F_grid_with_sidebands` calls and to `run_trajectory_sweep`. | cell 47 |
| 3e | New verification cell at end of notebook: load `zs_baselines_pre_refactor.json` from Phase 1f, compute current capture-velocity tables and F_max from the refactored cells, diff. Print zero-diff confirmation OR itemized deltas with attempted physics attribution. Fail-fast if any delta > 1% in F_max or any change in capture count. | new cell near end (before Notes & next steps markdown) |
| 3f | Update this spec doc with the post-refactor capture table and any flagged deltas. Commit. | docs/superpowers/specs/2026-05-13-... |

**Commit message:** `zeeman_slower_2: Phase 3 — extract sweep helpers, PB m_S-from-decoupled-basis, regression check`

## Open questions

(None blocking — defaults captured above. Recorded for review during execution.)

1. **`b_perp_to_s_remix` g_F default = 0.5.** This is a rough approximation for X²Σ⁺ N=1 in the case-(b) limit; the actual g_F varies by hyperfine sublevel from 0 to ~2/3. A better helper would accept the specific sublevel. For Phase 2b we use 0.5 as the order-of-magnitude default with the caveat in the docstring.
2. **`format_capture_table` config-detection.** Auto-handling both nested (3×2×4) and degenerate (1×1×4) shapes requires the helper to inspect the configs argument. Alternative: always pass `B_designs`, `s_values`, `L_values` as separate args, with degenerate ones as single-element lists. Cleaner; will use this in Phase 3a unless complications surface.
3. **`v_stop` capture classification.** Review 2 Important #3 flagged this. Investigation showed `solve_ivp`'s event root-finding handles the v=5 crossing exactly, and the binary "captured if v_stop fired" criterion is the right physical proxy for downstream MOT capture. No fix planned. Documented as a known-acceptable simplification in cell 35's docstring (added in Phase 1a alongside the `remix` docstring upgrade).
4. **`max_step=1e-4` in `simulate_trajectory`.** Review 2 Important #4 asked for a validation against `max_step=1e-5`. Will add as a one-shot validation cell in Phase 2 only if Phase 1f's baseline shows any sensitivity to integrator tolerance during the eventual refactor.

## Out of scope (explicit deferrals)

- Commit B (N=3 leakage / repumper): deferred to a separate spec. Will pick up `build_paschen_back_sigma_minus`'s `g_char_full['mS']`-based selection from Phase 3c as a prerequisite that's now satisfied.
- Full Liang-Yan 12-sideband (6 σ⁻ + 6 σ⁺) Type-II scheme: separate downstream work.
- Chirped-detuning slower variant: separate downstream work; the `run_trajectory_sweep` helper from Phase 3a will host it as a one-liner addition.
- Stochastic recoil heating in trajectory integrator: separate.
- Refactoring `Source Code/` — this entire spec is notebook-only. The verification gate from project CLAUDE.md (`jupyter execute RaF_Calcs_Tutorial.ipynb`) does NOT apply.
