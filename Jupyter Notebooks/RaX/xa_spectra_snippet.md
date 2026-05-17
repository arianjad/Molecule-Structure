# xa_spectra — copy-paste cell for `RaF X-A.ipynb` / `BaF X-A.ipynb`

These two notebooks hold your own uncommitted work and were intentionally
not edited. Paste the cell below where the old `simulate_spectra` /
`simulate_spectra_noM` calls were. It replaces the drifting copy-pasted
helpers with the shared module.

```python
import xa_spectra as xs   # config_path.add_to_sys_path() must have run already

# spectrum (auto: no-M -> TDM^2; M-resolved -> line strength S; no Boltzmann)
gidx = g.select_q({'N': 1})
eidx = e.select_q({'J': 0.5}, parity='+')
lines = xs.line_list(g, e, gidx, eidx, origin=e.parameters['Origin'])
xs.plot_spectrum(lines=lines, sticks=True,
                 broaden_kw=dict(shape='voigt', fwhm=40.0, lorentz_fwhm=10.0))

# excited-state-averaged branching (A hyperfine often unresolved)
avg_df, per_df = xs.averaged_branching(g, e, eidx)
display(avg_df)
```

Old → new map:
`simulate_spectra(...)`           → `line_list` + `plot_spectrum`
`simulate_spectra_noM(...)`       → `line_list` (no-M auto-detected)
`plot_gaussian_spectrum(...)`     → `broaden(shape='gaussian')` / `plot_spectrum`
manual offset/tweak/yscale block  → `plot_spectrum(experimental=dict(...))`
