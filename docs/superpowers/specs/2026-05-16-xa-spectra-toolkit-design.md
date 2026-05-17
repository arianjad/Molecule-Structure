# X–A spectra + averaged-BR toolkit — design

Status: designed 2026-05-16, awaiting user spec review.
Supersedes nothing. Builds on the validated BaF X–A pipeline
(`docs/superpowers/specs/2026-05-16-baf-xa-pipeline-design.md`).

## 1. Goal

Generalize the copy-pasted, drifting X–A spectrum helpers into one reusable,
physics-correct module, and add an excited-state-averaged branching ratio.
Drive it off the converter (not the hand G-row) for BaF A0.

Four deliverables:

1. **Item 3** — `Source Code/molecule_parameters.py` BaF A0: converter-driven
   `Origin`.
2. **`Jupyter Notebooks/RaX/xa_spectra.py`** — new molecule-agnostic module.
3. **`Jupyter Notebooks/RaX/BaF_spectrum_plot.py`** — regenerate its simulation
   from molecular structure via the module instead of pasted arrays.
4. Wire the cooling-notebook generator to the module; deliver a one-cell
   snippet for the user to wire `RaF X-A.ipynb` / `BaF X-A.ipynb` themselves.
   Those two notebooks carry uncommitted user WIP and are never touched here.

## 2. Item 3 — converter-driven Origin

Replace the hand G-row in BaF A0 with the converter:

```python
'Origin': 11946.109609,   # physical N2 band origin T0,0 (arXiv:2511.06986 Tbl III)
'formalism': 'N2',
'Lambda': 1,
```

Delete the `+ 6347.847/c` hand term. The load-time call site
(`molecule_parameters.py:426` / `Energy_Levels.py:76`) runs
`convert_formalism` because `formalism == 'N2'`, producing the R² backend.

Tag BaF X0 `'formalism': 'N2'`, `'Lambda': 0` for a consistent convention
declaration. This is cosmetic: the converter returns early for Λ²=0 (Σ states
are convention-independent), so X0 physics is unchanged.

**Verified deltas** (converter output minus the validated hand path, measured
this session on the live A0 dict):

| param | Δ (MHz) | source term |
|---|---|---|
| Origin | −0.006 | −Λ⁴D/c (hand version dropped it) |
| Be | −0.012 | −2Λ²D B-row |
| p+2q | −0.007 | +p2q_D Λ-doubling partner |
| ASO | +0.93 (on ~1.9×10⁷) | +A_D spin-orbit partner |

All shifts are far below the 1–2 MHz validation tolerances. The hand G-row was
an under-specified version of the R² convention; the converter applies the full
Table 7.2 transform self-consistently. Thesis reference: App A.3.1 "N and R
Formalisms" (p. 297) — R² is the backend; N² differs by the −K² electronic
origin offset plus the centrifugal corrections, exactly the converter's
+Be·Λ²/c G-row and Table 7.2 terms. Item-3 correctness is established by the
delta test above and the regression gate (§6), not by the citation.

## 3. Physics conventions (thesis-grounded)

Source: A. Jadbabaie, *Measuring Fundamental Symmetry Violation in Polyatomic
Molecules* (Springer Theses 2025), §3.2.4.1, pp. 107–115. Double prime =
ground (J″), single prime = excited (J′). `branching_ratios` and
`Calculate_TDMs` (Source Code) auto-diagonalize and auto-convert to case (a);
the module calls them and never reimplements a matrix element.

**Line strength (the spectrum stick height).** Thesis p. 109, definition of S:

> S = Σ_{M″,M′,p} |⟨g,i|T¹_p(d)|e,j⟩|²  — "the line strength of a transition."

- **M-resolved basis:** strength = S, summed over ground M″, excited M′, and
  polarization p ∈ {−1,0,+1}, built from `Calculate_TDMs` per q. The legacy
  overlapping-Lorentzian pile-up did the M′ sum implicitly; the module does it
  explicitly.
- **No-M basis (`M_sublevels='none'`):** strength = `branching_ratios` entry =
  |⟨g,J″|d|e,J′⟩|² = c(J″,J′)²·D_ge² (the J-dependent reduced ME squared).

**Frequency.** ν = (E_e − E_g) + c·`Origin`. Default `Origin` =
`Excited.parameters['Origin']`. Units MHz (default) or cm⁻¹.

**ω₀³ is deliberately omitted.** Thesis Eqs. 3.1/3.3/3.6 carry
ω₀³/(3πε₀ℏc³); p. 109 assumes ω_ij³ ≈ ω₀³ within a band. The factor cancels
identically in any branching ratio (Eq. 3.6: γ_ij = r_ij·γ) and varies <10⁻⁵
across one X–A band, so bare |d|² is the correct relative-intensity weight
in-band. It matters only across widely separated bands (footnote 20: YbOH
v″=1/0 gives an ω³ ratio ≈ 0.92). `line_list` exposes
`weight='dipole2'` (default) or `weight='rate'` (×ω₀³) for the cross-band case.

**Boltzmann population is removed; Doppler broadening is not.** The thesis
separates population weighting from lineshape. The legacy
`exp(-E_g·1.439/T)` population factor and the ad-hoc `J_adjust`/`P_adjust`
(an informal stand-in for the Eq. 3.11 (2J′+1)/(2J″+1) degeneracy ratio) are
both dropped. The explicit M-sum / reduced-ME path replaces them.

**Averaged branching ratio.** Thesis Eq. 3.4:

> r̃_{J″,J′} = |⟨g,J″‖d‖e,J′⟩|² / [(2J′+1)·D_ge²],  normalized Σ_{J″} r̃_{J′,J″} = 1.

The (2J′+1)⁻¹ is the average over initial (excited) orientations. `averaged_
branching` generalizes this to hyperfine: average the per-excited-sublevel
normalized BR over every excited F and M_F with equal weight per sublevel.

**Lineshape.** Thesis Eq. 3.12: σ(Δ) ∝ r_ij·D_ge²·g(Δ) — strength × normalized
lineshape. Doppler (Eqs. 3.14–3.15): Γ_std,D = (ω₀/c)·√(k_BT/M),
Γ_FWHM = 2√(2 ln 2)·Γ_std. Voigt = Doppler ⊗ Lorentzian when mechanisms are
comparable (p. 112). Thesis p. 111 uses two lineshape normalizations:
∫g dω = 1 (area) and ĝ(0) = 1 (peak). The module supports both.

## 4. Module API — `Jupyter Notebooks/RaX/xa_spectra.py`

Molecule-agnostic. Imports only from `Energy_Levels` (`branching_ratios`,
`Calculate_TDMs`). Caller passes initialized `MoleculeLevels` objects.

```python
line_list(Ground, Excited, g_idx, e_idx, *,
          origin=None, units='MHz', field=(0.0, 0.0),
          pol='all', weight='dipole2', thresh=1e-9) -> pandas.DataFrame
    # columns: freq, strength, plus dominant g/e labels (N,J,F[,M] / J,F,parity[,M])
    # strength: no-M -> branching_ratios entry; M -> S (Sec. 3), summed over q in pol
    # NO Boltzmann, NO J_adjust. Empty result -> empty DataFrame + warning.

doppler_fwhm(T, mass_amu, nu0) -> float                # thesis Eq. 3.15
    # FWHM returned in the same frequency unit as nu0 (Hz in -> Hz out);
    # Gamma_std = (nu0/c) sqrt(kB T / M), Gamma_FWHM = 2 sqrt(2 ln 2) Gamma_std

broaden(lines, *, shape='gaussian', fwhm=None,
        T=None, mass_amu=None, omega0_Hz=None,
        x=None, norm='area', cluster=True, cluster_thresh=None,
        n_per_fwhm=20, pad=5) -> (x, y_total, (comp_x, comp_y))
    # shape in {gaussian, lorentzian, voigt}; fwhm explicit OR from (T,mass,omega0)
    # norm in {area (int g dw = 1), peak (g_hat(0) = 1)}
    # cluster: merge sub-cluster_thresh lines (mean freq, summed strength)

plot_spectrum(lines=None, *, Ground=None, Excited=None,
              g_idx=None, e_idx=None, ax=None,
              sticks=True, broaden_kw=None, normalize=True,
              experimental=None, label='Simulation',
              **line_list_kw) -> matplotlib Axes
    # experimental = dict(freq, signal, err=None, offset=0, tweak=0, yscale=1)
    #   -> generalizes the BaF_spectrum_plot.py offset/tweak/yscale overlay

averaged_branching(Ground, Excited, e_select, *,
                    field=(0.0, 0.0), normalize=True) -> (avg_df, per_sublevel_df)
    # avg_df: BR per ground level, averaged over every excited F and M_F
    #   (equal weight per excited sublevel) == thesis Eq. 3.4, hyperfine-generalized
    # per_sublevel_df: the unaveraged per-excited-sublevel table
```

`pol`: `'all'` sums intensities over q ∈ {−1,0,+1}; `'z'`, `'+'`, `'-'`, `'x'`
select subsets. Summing intensities (|·|²), not amplitudes — a correctness fix
over the legacy `select_dipole('all')` amplitude sum.

## 5. Wiring

- **`BaF_spectrum_plot.py`**: build BaF X (N=0,1) and A (J=½) via
  `MoleculeLevels.initialize_state`, call `plot_spectrum` with the existing
  08/27/25 `fluor_data` as `experimental=`. No pasted `sim_freq`/`sim_data`.
- **`make_baf_xa_cooling_nb.py`**: the BR table/figure call `xa_spectra`;
  physics output matches the validated `baf_xa_validate.py` numbers within the
  §6 regression tolerances. Regenerate `baf_xa_cooling.ipynb`.
- **Snippet for the user**: a single import-and-call cell, delivered as
  `Jupyter Notebooks/RaX/xa_spectra_snippet.md`, for the user to paste into
  `RaF X-A.ipynb` / `BaF X-A.ipynb`. Those files are not edited here.

## 6. Verification

- **Regression (load-bearing).** `xa_spectra` no-M path reproduces
  `baf_xa_validate.py`: Table II (−)/N=0 frequencies ≤ 2 MHz; cooling-line BR
  A²Π₁/₂(J=½,+)→X N=1 parity-closed (Σ→N=1 ≥ 0.999, Σ→N=0 ≤ 1e-3).
- **Physics A.** `line_list` no-M strength equals the `branching_ratios` entry
  exactly (definitional).
- **Physics B.** `averaged_branching` on BaF A²Π₁/₂(J=½)→X N=1 satisfies
  Σ over ground = 1 and equals an equal-weight average of the validated
  per-sublevel cooling table.
- **Physics C.** M-resolved S, normalized per excited column, agrees with the
  no-M branching distribution for BaF X(N=1)←A(J=½,+) up to the c(J″,J′)²
  reduced-ME/degeneracy factor (Wigner-Eckart consistency between bases).
- **Item-3 gates.** `baf_xa_validate.py` exits 0 (stays green vs the
  converter path). The Source-Code change triggers the root CLAUDE.md gate:
  `conda run -n Structure jupyter execute "Jupyter Notebooks/RaF_Calcs_Tutorial.ipynb"`
  exits 0, no `--allow-errors`.
- **Notebook smoke.** `BaF_spectrum_plot.py` runs and renders the
  structure-derived simulation plus the experimental overlay.

Test file: `Jupyter Notebooks/RaX/test_xa_spectra.py` (hard-assert, exit 0),
patterned on `baf_xa_validate.py`.

## 7. Non-goals

- No thermal/Boltzmann population model (removed per the spectrum definition).
- No A²Π₁/₂–B²Σ⁺ perturbation or intensity modeling — the single-state
  idealization of the prior spec §3.6 stands.
- No edits to `RaF X-A.ipynb` / `BaF X-A.ipynb` (uncommitted user WIP).
- No `config_path` packaging; no constant refits; no sextic/H-order terms
  (the converter is quartic-truncated).

## 8. Files

| Action | Path |
|---|---|
| Modify | `Source Code/molecule_parameters.py` (BaF A0 + X0) |
| Create | `Jupyter Notebooks/RaX/xa_spectra.py` |
| Modify | `Jupyter Notebooks/RaX/BaF_spectrum_plot.py` |
| Modify | `Jupyter Notebooks/RaX/make_baf_xa_cooling_nb.py` (regenerate `baf_xa_cooling.ipynb`) |
| Create | `Jupyter Notebooks/RaX/test_xa_spectra.py` |
| Create | `Jupyter Notebooks/RaX/xa_spectra_snippet.md` (handed to user) |
