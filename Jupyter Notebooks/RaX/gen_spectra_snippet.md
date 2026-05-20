# gen_spectra — copy-paste cell for `RaF X-A.ipynb` / `BaF X-A.ipynb`

These two notebooks hold your own uncommitted work and were intentionally
not edited. Paste the cell below where the old `simulate_spectra` /
`simulate_spectra_noM` calls were. It replaces the drifting copy-pasted
helpers with the shared module.

```python
import gen_spectra as gs   # config_path.add_to_sys_path() must have run already

# spectrum: line strength S (= LIF excitation signal) by default; same
# observable whether the states are M-resolved or not. For the emission
# branching add initial='excited', initial_reduction='average'; for the
# per-molecule absorption cross section initial='ground', ='average'.
gidx = g.select_q({'N': 1})
eidx = e.select_q({'J': 0.5}, parity='+')
lines = gs.line_list(g, e, gidx, eidx, origin=e.parameters['Origin'])
gs.plot_spectrum(lines=lines, sticks=True,
                 broaden_kw=dict(shape='voigt', fwhm=40.0, lorentz_fwhm=10.0))

# excited-state-averaged branching (A hyperfine often unresolved)
avg_df, per_df = gs.averaged_branching(g, e, eidx)
display(avg_df)
```

Notes:
- FWHM units in `broaden_kw` match the `units` kwarg on `line_list`
  (default `units='MHz'`; pass `units='cm-1'` for the cm⁻¹ branch).
- Boltzmann population is off by default. Pass `boltzmann_T=<T_K>` to
  `line_list` to scale strengths by `exp(-E_g/(kB*T))`; T in Kelvin.
- For BR observable (matches legacy `simulate_spectra_noM`):
  `gs.line_list(..., initial='excited', initial_reduction='average')`.

Old → new map:
`simulate_spectra(..., T=)`        → `line_list(..., boltzmann_T=)` + `plot_spectrum`
`simulate_spectra_noM(...)`        → `line_list(..., initial='excited', initial_reduction='average')`
`plot_gaussian_spectrum(...)`      → `broaden(shape='gaussian')` / `plot_spectrum`
manual offset/tweak/yscale block   → `plot_spectrum(experimental=dict(...))`
