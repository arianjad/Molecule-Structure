# `zeeman_slower_2.ipynb` — deferred extensions spec

Captured 2026-05-10. Follow-ups to the BR/sidebands/leakage refactor that landed today.
Commit log so far on `fix-for-distrib`:

| SHA | Subject |
|---|---|
| `57e409e` | Doppler sign + saturation form fixes |
| `7204830` | re-execute |
| `f50fb73` | Multi-sideband laser API (Commit A) |
| `640cf4c` | Rotational leakage plot — A(J=½,+) → X(N=3) vs B |

This spec captures Commits B/C/D from the prior plan plus a couple of related extensions surfaced in conversation.

## Commit B — Expanded Hamiltonian (X N=1,3; A J=½ and J=3/2 from start)

**Goal**: have all the relevant cycling + leakage + repumper states in `g`/`e` from the start, so all downstream diagnostics (Plot A through F_max(B)) display the broader picture without per-cell extended-basis hackery (which is what the new leakage cell currently does).

**Changes**:
- Cell 4: `N_g = [1, 3]` and keep `N_e = [1, 2]` (already covers J=½ and J=3/2 of A).
- Cell 12 `subset_TM` calls (and downstream): the "default" cycling subset stays X(N=1, parity=−) → A(J=½, parity=+); add a second subset for the repumper X(N=3, parity=−) → A(J=3/2, parity=+).
- The leakage cell (current cell 31) collapses to a slice of the canonical TM_full — drop `g_ext`/`TM_ext` boilerplate.
- Plot A/B/C/D/E auto-scale to the larger transition map. Manually verify they remain readable; may need to filter by polarization or use facetting (per-J or per-N) if too dense.
- Cell 8 closure check — should still pass at 1.000 (TDM2_BR sums correctly per excited state regardless of subset).

**Cost**: ~2-3× slower on cell 7 (`build_transition_map` does an `e.size × g.size` matmul; X N=1,3 with I=½ HF gives g.size=40, vs current 12). Notebook runtime probably grows from ~15s → ~30-40s. Acceptable.

**Deferred questions**:
- Should we include both Λ-doublet parities of A (J=½, J=3/2 ±parity)? Useful for completeness; doubles the excited basis.
- Plot D (BR map) becomes hard to read with 40 ground states × 8 excited states; consider per-N facetting.

## Commit C — Full multi-sideband Zeeman slower config + repumper

**Goal**: realize the actual MOT/slower laser configuration with all four hyperfine sideband frequencies addressing X(N=1) → A(J=½, +parity) plus an optional repumper sideband from X(N=3) → A(J=3/2, +parity).

**Changes**:
- New cell defining a `slower_sidebands_RaF` list: 4 sidebands at the X(N=1) hyperfine offsets (J=½ F=0, J=½ F=1, J=3/2 F=1, J=3/2 F=2 — frequencies relative to the brightest channel at low B). Polarizations chosen per Tarbutt's "two-frequency MOT" trick: alternating σ⁺/σ⁻ with one detuning per pair. Initial intensities: each at `s = 1` per Tarbutt 2015 numbers, total ~4·I_s ≈ 5.6 mW/cm² per beam.
- A `repumper_sidebands` list addressing X(N=3) hyperfine groups. Polarization either σ⁺ or remixed; intensity comparable to main slower light.
- Re-run F(v) plot with combined sidebands; should show much-broader velocity acceptance window than single-sideband case.
- Re-run F_max(B) — peak should remain high to higher B because the repumper unblocks the leakage channel.
- Add a "scattering rate per scatter incl. repumper" diagnostic — at high B, what fraction of cycle photons go through the main vs repumper transitions? Sets the repumper power requirement.

**Deferred questions**:
- Decreasing-B (Bz from high → low) or increasing-B (low → high) slower? The sign of the slope determines polarization choices (σ⁺ vs σ⁻). Set by user's experimental geometry.
- How precisely does the user know A(J=3/2, +parity) hyperfine structure? RaF measurements probably exist (Garcia Ruiz group) but should double-check before locking in repumper sideband detunings. Currently in `molecule_parameters.py`: `'h1/2': 0` for A0 — implies hyperfine collapsed in the model. Either accept this approximation or add `h3/2` from theory/experiment.

## Commit D — Trajectory simulation

**Goal**: numerically integrate $m \dot v = F(v(t), B(z(t)))$ along $z(t) = \int v\,dt$ for a designed B(z) profile and a chosen sideband configuration. Outputs phase-space trajectory `(z, v)` and velocity-vs-time `v(t)`. This is the load-bearing physical test of whether the slower works.

**Changes**:
- New utility function `simulate_trajectory(initial_v, B_of_z_func, sidebands, dt, t_max, mass)`: RK4 (or `scipy.integrate.solve_ivp` with adaptive step) integrator. Returns arrays `t, z, v`.
- B(z) profile: either (a) the linearized profile from `required_Bz` (cell 19), (b) an analytic shape (Bessel-like), (c) a tabulated experimental field. Build a `make_B_of_z(...)` that returns a callable.
- Plot 1: `v(z)` overlaid against the kinematic `v_target(z) = sqrt(v_i² - 2η·a_max·z)` from `slowing_profile`.
- Plot 2: `v(t)` for several initial velocities; identify the capture velocity where the molecule fails to track B(z).
- Optional: add stochastic recoil heating (random momentum kicks per emitted photon) to estimate residual velocity spread.

**Deferred questions**:
- Tarbutt's rate equations track populations adiabatically; does that hold for the timescales of slowing? The molecule traverses 1 m at 100 m/s in 10 ms; with Γ ~ 4.6 MHz, populations equilibrate in ~1 µs. So adiabatic is fine for ~10000× slower dynamics. ✓
- Should we model laser intensity as constant along the slower or use a profile (Gaussian beam waist along z)? Constant is fine for first pass; profile matters for capture velocity tail.

## Tarbutt-vs-thesis convention discrepancies (informational, not action items)

Documented during the saturation research:

- Tarbutt 2014 (NJP 17, 015007), Eq 2: rate denominator is `1 + 4Δ²/Γ²` (linear, no saturation in denom). Population balance handles saturation.
- Steck Eq 5.250 / your thesis Eq 3.22: closed 2-level OBE steady state has `1 + s + 4Δ²/Γ²` — but this is the saturated population, NOT to be used inside a rate-equation framework.
- Notebook (post-fix `57e409e`): follows Tarbutt convention. ✓

Other minor convention notes:
- $I_s$ scale convention factor of 3: thesis Eq 3.21 uses `π h c Γ/(3λ³)` (isotropic-orientation averaged); Steck Eq 5.253 uses `π h c Γ/λ³` (linearly polarized random orientation). Differ by factor of 3. Notebook uses thesis convention. For RaF: $I_s ≈ 1.41$ mW/cm².
- Tarbutt's `f_{l,u,p}`, your thesis `r_{ij}`, my notebook's `b_{lu,p}` (= `TDM2_BR`) are all the same thing: per-channel BR with denominator = sum over all (p, g_full) per excited state.

## Open physics questions deferred

- **Closure of A(J=3/2, +parity) cycle.** When we add the repumper, A(J=3/2, +parity) excited states decay to X(N=1, −) AND X(N=3, −) AND X(N=5, −). Is the BR back to N=1, 3 high enough to keep the repumper closed without further repumpers (N=5)? Need to compute, similar to the J=½+ → N=3 leakage but for J=3/2+ → N=5.
- **Multi-color slower configurations.** With both A(J=½,+) and A(J=3/2,+) addressed simultaneously by separate sidebands, can we find a "chirped multi-color" slower that handles a wider velocity class than the single-color version? Use `find_doppler_trackable_bands` on the union of transitions.
- **Photon recoil spread.** Once trajectory sim is in place, add Poissonian scattering noise per integration step to estimate the final velocity spread. Compare to recoil temperature limit.
- **Saturation form for inhomogeneous broadening.** Steck Eq 5.134 has a different form for inhomogeneous broadening (collisions): `s = (Ω²/γ⊥Γ)/(1 + Δ²/γ⊥²)`. Probably not relevant for cold beams (no collisions, narrow Doppler) but worth flagging.

## Outstanding deferred items from prior conversation

- The thesis Eq 3.21 has a vibrational-branching version of `r̃` scaled by Franck-Condon factors. Currently the notebook ignores vibrational branching (`vib_state=0` only). For RaF $q_{v=0\to v=0} \approx 0.99$ this is < 1% error; for other molecules (worse FC) it would matter.
- The `BR_e2g` in `steady_state_rate_eq` (cell 27) renormalizes spontaneous emission within the SUBSET to sum to 1. This is the closed-cycle assumption. With Commit B (full Hamiltonian + repumper), the subset is closed by construction so this approximation becomes exact.
- The dark-state count in Plot E uses a 0.01 BR threshold. With Commit B the BR scale changes (more ground states means smaller per-state BR); may need to re-tune the threshold.
