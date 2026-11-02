"""Zeeman-shift plotting helpers, extracted from RaF X-A.ipynb.

Caller must have run `state.ZeemanMap(Bz, ...)` first so `state.Bz` and
`state.evals_B` are populated. See plot_zeeman_shifts docstring for details.

Sibling module to gen_spectra; same import pattern (callers put "Source Code"
on sys.path via config_path first).
"""
import numpy as np
import matplotlib.pyplot as plt


def plot_zeeman_shifts(state, queries, labels=None, parity=None,
                       relative=True, baseline='mean',
                       kG=False, GHz=False,
                       derivative=0, g_factor=False,
                       ncols=None, figsize_per=(5, 4)):
    """Plot Zeeman shifts (or their derivatives) for groups of levels.

    Pulls eigenvalues from a `MoleculeLevels` object that has already had
    `ZeemanMap` run on it, slices out the levels matching each query, and
    draws one subplot per query in a grid.

    Parameters
    ----------
    state : MoleculeLevels
        Ground- or excited-state object (e.g. `g`, `e`). Must have run
        `state.ZeemanMap(Bz, ...)` so that `state.Bz` and `state.evals_B`
        are populated. Eigenvalues are in MHz, field in Gauss.
    queries : list
        One entry per panel. Each entry is either:
            - a dict passed straight to `state.select_q`, e.g.
              `{'N': [1]}` or `{'J': 0.5, 'M': 0}`, or
            - a precomputed integer index array (e.g. the return value
              of `state.select_q(...)`).
    labels : list of str, optional
        Subplot titles. Defaults to `repr(query)` for each entry.
    parity : {'+', '-', None}, optional
        Forwarded to `state.select_q` for every dict query. Ignored for
        index-array queries; if you need mixed parities, pass index
        arrays. Default `None`.
    relative : bool, optional
        Only used when `derivative == 0`. If True, subtract a common
        scalar baseline per panel so that the B = 0 splitting remains
        visible (subtracting per-line baselines would collapse every
        level to start at zero). Default True.
    baseline : {'mean', 'min', float, None}, optional
        How to pick the scalar baseline (only when `relative=True`):
            - 'mean' : mean of the group's zero-field eigenvalues.
            - 'min'  : lowest level in the group at B = 0.
            - float  : explicit reference energy in MHz, e.g.
                       `g.evals_B[0, g.select_q({'N':[0]})[0]]`.
            - None   : no shift (equivalent to `relative=False`).
        Default 'mean'.
    kG : bool, optional
        Plot field in kGauss instead of Gauss. Default False.
    GHz : bool, optional
        Plot energy in GHz instead of MHz. Default False.
    derivative : {0, 1, 2}, optional
        0 plots the energy itself; 1 plots dE/dB (slope); 2 plots
        d²E/dB² (curvature). Computed with `np.gradient` on the actual
        `Bz` grid, so non-uniform spacing is fine; expect edge artifacts
        at the first and last B point. `relative` / `baseline` are
        ignored when `derivative > 0`.
    g_factor : bool, optional
        Only meaningful with `derivative == 1`. Divides the slope by the
        Bohr magneton (1.39962 MHz/G) so the y-axis is the dimensionless
        effective g-factor. Default False.
    ncols : int, optional
        Subplot grid columns. Default `ceil(sqrt(len(queries)))` for an
        approximately square layout.
    figsize_per : (float, float), optional
        Figure size per subplot in inches. Default (5, 4); the total
        figure is `(ncols * w, nrows * h)`.

    Returns
    -------
    fig : matplotlib.figure.Figure
    axes : ndarray of matplotlib.axes.Axes, shape (nrows, ncols)
        Always 2-D (squeeze=False). Unused cells are hidden with
        `axis('off')`.
    panels : list of dict
        One entry per query, in the same order. Each dict contains the
        data actually drawn on that subplot, so it can be saved or
        replotted without re-running ZeemanMap:
            'x'     : 1-D ndarray, field axis as plotted (Gauss or kGauss)
            'y'     : 2-D ndarray, shape (N_B, n_selected), values as
                      plotted (after derivative, baseline, unit scaling)
            'idx'   : 1-D ndarray of int, indices into state.evals_B
            'query' : the original query dict / index array
            'label' : str, subplot title used
            'ax'    : the matplotlib Axes object for this panel
        Individual `Line2D` objects are also accessible via `ax.lines`
        if you need to tweak per-line styling after the fact.

    Notes
    -----
    - At avoided crossings the slope and curvature change rapidly; make
      `Bz` dense enough or those panels will look noisy.
    - Eigenvalues come from `evals_B` after the state-ordering step in
      `ZeemanMap`, so each column is a single adiabatic level vs B.
    - For an overlay (no subplots), pass a single query and ignore the
      grid, or call `plt.plot(state.Bz, state.evals_B[:, idx])` directly.

    Examples
    --------
    >>> # Energy levels, mean baseline, MHz vs Gauss
    >>> plot_zeeman_shifts(g, [{'N':[0]}, {'N':[1]}, {'N':[2]}])

    >>> # Slope dE/dB
    >>> plot_zeeman_shifts(g, [{'N':[1]}], derivative=1)

    >>> # Effective g-factor
    >>> plot_zeeman_shifts(g, [{'N':[1]}], derivative=1, g_factor=True)

    >>> # Reference everything in N=1 to the lowest M=0 level at B=0
    >>> ref = g.evals_B[0, g.select_q({'N':[1], 'M':0})[0]]
    >>> plot_zeeman_shifts(g, [{'N':[1]}], baseline=ref)

    >>> # Excited state, positive parity, J=1/2 and J=3/2
    >>> plot_zeeman_shifts(e, [{'J':0.5}, {'J':1.5}], parity='+')
    """
    Bz    = state.Bz
    evals = state.evals_B                       # (N_B, n_levels), MHz

    traces_full = evals.copy()
    for _ in range(derivative):
        traces_full = np.gradient(traces_full, Bz, axis=0)

    x       = Bz * (1e-3 if kG else 1)
    yscale  = 1e-3 if GHz else 1
    xlabel  = 'B (kGauss)' if kG else 'B (Gauss)'

    if derivative == 0:
        ylabel = ('ΔE' if relative else 'E') + (' (GHz)' if GHz else ' (MHz)')
    elif derivative == 1 and g_factor:
        traces_full = traces_full / 1.39962
        ylabel = r'$g_\mathrm{eff}$'
        yscale = 1                                 # dimensionless: don't rescale
    else:
        unit_num = 'GHz' if GHz else 'MHz'
        unit_den = 'kG'  if kG  else 'G'
        if kG:
            traces_full = traces_full * 1e3 ** derivative
        if derivative == 1:
            ylabel = rf'$dE/dB$ ({unit_num}/{unit_den})'
        else:
            ylabel = rf'$d^{{{derivative}}}E/dB^{{{derivative}}}$ ({unit_num}/{unit_den}$^{{{derivative}}}$)'

    idx_list = [state.select_q(q, parity=parity) if isinstance(q, dict)
                else np.asarray(q) for q in queries]

    n = len(idx_list)
    if ncols is None:
        ncols = int(np.ceil(np.sqrt(n)))
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, squeeze=False,
                             figsize=(figsize_per[0]*ncols,
                                      figsize_per[1]*nrows))

    panels = []
    for k, idx in enumerate(idx_list):
        ax = axes[k // ncols, k % ncols]
        title = labels[k] if labels else str(queries[k])
        if len(idx) == 0:
            ax.set_title(f'{queries[k]}  (no match)'); ax.axis('off')
            panels.append({'x': x, 'y': np.empty((len(x), 0)),
                           'idx': np.asarray(idx, dtype=int),
                           'query': queries[k], 'label': title, 'ax': ax})
            continue
        traces = traces_full[:, idx]
        if derivative == 0 and relative:
            ref = evals[0, idx]
            if baseline == 'mean':
                b = ref.mean()
            elif baseline == 'min':
                b = ref.min()
            elif np.isscalar(baseline):
                b = float(baseline)
            else:
                b = 0.0
            traces = traces - b
        y_plot = yscale * traces
        ax.plot(x, y_plot, lw=1)
        ax.set_title(title, fontsize=11)
        ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
        ax.axhline(0, color='k', lw=0.4, alpha=0.4)
        panels.append({'x': x, 'y': y_plot,
                       'idx': np.asarray(idx, dtype=int),
                       'query': queries[k], 'label': title, 'ax': ax})

    for k in range(n, nrows*ncols):
        axes[k // ncols, k % ncols].axis('off')

    if derivative == 1 and g_factor:
        title_str = r'effective $g$-factor'
    else:
        title_bits = {0: 'Zeeman shifts', 1: 'dE/dB', 2: r'd$^2$E/dB$^2$'}
        title_str = title_bits.get(derivative, f'd^{derivative}E/dB^{derivative}')
    fig.suptitle(f'{state.state_str} ' + title_str, fontsize=14)
    fig.tight_layout()
    return fig, axes, panels
