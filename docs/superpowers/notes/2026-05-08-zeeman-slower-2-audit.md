# `zeeman_slower_2.ipynb` audit — 2026-05-08

Read of [`Jupyter Notebooks/RaX/zeeman_slower_2.ipynb`](../../../Jupyter Notebooks/RaX/zeeman_slower_2.ipynb) (untracked, 31 cells) for physics issues and bugs. Companion implementation plan: [`../plans/2026-05-08-zeeman-slower-2-audit-fixes.md`](../plans/2026-05-08-zeeman-slower-2-audit-fixes.md).

## Bugs / errors

**B1. `a_recoil` is dimensionally wrong** (cell 2)
```python
v_recoil = hbar_J * k_XA / mass_RaF        # m/s per scattered photon
a_recoil = hbar_J * k_XA / mass_RaF        # m/s^2 / scattering rate factor
```
Same RHS, but `a_recoil` is labeled m/s². Acceleration is `ℏk·R_scatter/m`, not `ℏk/m`. Variable is unused later (`a_max_est = ℏk·Γ/(2m)` in cell 19 is computed independently and correctly). Delete or fix.

**B2. `Calculate_TDM_evecs` print-spam loop** (cell 7)
303 calls to `Calculate_TDM_evecs` (3 polarizations × 101 B-points), each internally re-converting case-(b) ground evecs to aBJ and printing two "Successfully converted eigenvectors" lines. ~600 prints + ~10× slower than necessary. Same bug as the helper notebook; fix is identical: hoist the case-conversion outside the (B, p) loop and wrap with `contextlib.redirect_stdout(io.StringIO())`. Inner loop becomes a single batched `E_a @ TDM_matrix @ G_a.transpose(0,2,1)` per p.

**B3. `find_doppler_trackable_bands` clusters by slope only, not frequency** (cell 17)
Two transitions at very different frequencies (e.g., 600 MHz apart) but matching `dν/dB` get grouped. For "single laser tracks B" the cluster is meaningless; for "each transition has its own laser-sideband but shares a B(z) profile" the cluster is fine — the markdown's wording suggests the former, the math implements the latter. Either re-word the description or add a frequency-coincidence check.

**B4. Closed-subspace BR silently inflates the force** (cell 27)
```python
weights_e2g = TDM2[iB].sum(axis=0)        # subset only
Z = weights_e2g.sum(axis=1, keepdims=True)
BR_e2g = weights_e2g / Z                  # forced to 1 within subspace
```
Divides by the subset-only normalization, so spontaneous decay is renormalized back into the cycle. If `leakage` from cell 22 is e.g. 5%, the rate-eq force is over-estimated by ~Γ_leak/Γ. Markdown notes this is the assumption, but `steady_state_rate_eq` should at minimum print/assert the leakage value at the queried B and warn above some threshold (≥1%).

**B5. `B_SR_PB`/`B_HF_PB` shading is misleading** (cell 15, defaults `(50,200)` and `(500,2500)`)
Original audit text understated the constants. Actual values from [`Source Code/molecule_parameters.py`](../../../Source Code/molecule_parameters.py) for ²²⁶RaF X0 are `Gamma_SR = 175.38 MHz`, `bF = 96.3 MHz`. PB onsets (μ_B·B ≈ coupling) are ~125 G (SR) and ~69 G (HF). Original `B_SR_PB=(50, 200)` is roughly correct; `B_HF_PB=(500, 2500)` is way too high — fully decoupled by ~500 G. Replace with values derived from the constants.

## Methodological caveats (not bugs but worth noting)

**C1. Laser auto-retunes per B-query** (cell 27)
`omega_L_MHz = dE[iB, je0, ig0] + detuning_MHz` is recomputed at each B. So `F(v, B)` plots represent "perfectly tuned laser per B", not a fixed-frequency laser sweep. Fine for diagnosing F_max envelope vs B, wrong for modeling an actual slower experiment with fixed lasers. Add caveat to cell-26 markdown.

**C2. Reference normalization label** (cell 8)
`TDM2_ref = TDM2[0].max()` is the max over the (subset, polarization) at the lowest B in the grid (1 mG). Calling it "case-(b) stretched ME" is approximate — A²Π_{1/2} is case-(a)/aBJ, so what gets picked up is whichever (b)→(a) mixed pair maximizes |TDM|² at that B. Numerically usually close, but the *label* is unverified.

**C3. `Bz[0] = 1e-3 G`** (cell 5)
Avoids the diagonalization singularity at B=0 — fine — but `np.gradient` at the endpoint uses a one-sided difference over the (1 mG, 10 G) interval, which spans the entire SR coupling regime. Slope at index 0 is unreliable; consider starting at e.g. 0.1 G and discarding `slope[0]` from clustering, or extending Bz with a couple of low-B points.

**C4. Unused import** `branching_ratios` (cell 2) — drop or use for a cross-check.

## Physics correctness — verified ✓

| Item | Check | Status |
|---|---|---|
| Polarization mapping `[σ⁻, π, σ⁺] ↔ [-1, 0, +1]` | Standard for absorption Δm = M_e − M_g = p | ✓ |
| q_electronic = (−1, +1) for X²Σ → A²Π_⊥ | Coherent q-sum inside `TDM_p_builders` (confirmed in [`Source Code/Energy_Levels.py:1287`](../../../Source Code/Energy_Levels.py#L1287)) | ✓ |
| Doppler sign for counter-propagating beam | `Δ = ω_L − ω_eg − kv` ≡ counter-prop. beam at −ẑ vs molecule at +ẑ | ✓ |
| Force sign `F = −ℏk·R` | Beam in −ẑ → kick in −ẑ → decelerates +v_z | ✓ |
| Rate-matrix structure | Σ_u r_ul net rate out of l, +BR_ul·Γ feedback, total decay −Γ from u | ✓ |
| Saturation lineshape | `R = (Γ/2)·sb / (1 + 4Δ²/Γ² + 2sb)` matches Tarbutt 2018 | ✓ |
| `pop_diff = N_l − N_u` for stimulated force | Net absorption − stim-emission per cycle | ✓ |
| `remix` term applies to lineshape only, not to force | `s_vec_force` excludes remix → no spurious thrust from fictitious dressing | ✓ |
| `v_recoil ≈ 2.16 mm/s` | ℏk/m for ²²⁶Ra¹⁹F at 753 nm | ✓ |
| `mu_B = 1.399624494 MHz/G` | NIST | ✓ |

## Top fixes, ordered (mapped to plan tasks)

| # | Audit ID | Plan Task |
|---|---|---|
| 1 | B2 (print spam) | Task 1 |
| 2 | B4 (leakage diagnostic) | Task 2 |
| 3 | B3 (cluster criterion) | Task 3 |
| 4 | B5 (PB shading) | Task 4 |
| 5 | B1 (a_recoil) + C1 (caveat) + C2 (label) + C4 (import) | Task 5 |
| 6 | end-to-end execute | Task 6 |

C3 (endpoint slope) is explicitly NOT addressed — minor cosmetic.

No physics-correctness errors that change conclusions. Slope/clustering/force results are quantitatively trustworthy modulo the closed-subspace BR over-count (B4).
