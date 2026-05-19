"""Molecule-agnostic spectrum / line-strength toolkit (gen_spectra).

One primitive -- the line strength
  S = sum_{M'',M',p} |<g,i|T^1_p(d)|e,j>|^2          (thesis Eq.3.1 / p.109)
computed identically whether the states were built M-resolved or not. Two
orthogonal knobs select the physical observable:

  initial            : 'ground' (default) | 'excited'   -- the initial state
  initial_reduction  : 'sum'     -> S                    line strength = LIF
                                    excitation signal (the thermal (2F''+1)
                                    ground population is already inside S)
                       'average' -> S / (2F_initial+1):
                          initial='excited' -> emission branching ratio
                                    (== Energy_Levels.branching_ratios)
                          initial='ground'  -> per-molecule absorption cross
                                    section (thesis Eq.3.11)

cross_section and branching are the SAME operation (divide by the initial-state
degeneracy) with the initial state swapped; line strength is the un-reduced
sum. M-resolution is only a basis choice, not an observable switch -- the same
(initial, initial_reduction) gives the same physics in either basis.

Within a band omega_0^3 is dropped (cancels in any ratio; <1e-5 across one
band) -- bare |d|^2; weight='rate' restores x nu^3 for the cross-band case.
NO Boltzmann factor (hyperfine << kT). F is good at zero field, so the
(2F+1) reductions are exact there; for finite field use an M-resolved build.

Callers put "Source Code" on sys.path first (notebooks: config_path).
This module never reimplements a matrix element -- it calls
Energy_Levels.branching_ratios / Calculate_TDMs.
"""
import warnings
import numpy as np
import pandas as pd
from Energy_Levels import branching_ratios, Calculate_TDMs

C_CM = 29979.2458                       # MHz per cm^-1 (== molecule_parameters.c)
_PARITY = {1: '+', -1: '-'}
_POL = {'all': [-1, 0, 1], 'z': [0], '+': [1], '-': [-1], 'x': [-1, 1]}
LINE_COLUMNS = ['freq', 'strength',
                'g_N', 'g_J', 'g_F', 'g_par', 'e_J', 'e_F', 'e_par',
                'g_idx', 'e_idx']


def _has_M(state):
    """True iff the state was built with resolved M sublevels."""
    M = getattr(state, 'q_numbers', {}).get('M', None)
    if M is None:
        return False
    return np.unique(np.asarray(M, dtype=float)).size > 1


def _dom(state, i, key):
    """Dominant basis label of quantum `key` for eigenstate i."""
    qd = int(np.argmax(state.evecs0[i] ** 2))
    return state.q_numbers[key][qd]


def _twoF1(state):
    """(size,) array of (2F+1) per eigenstate (dominant F; F is good at B=0)."""
    return np.array([2.0 * float(_dom(state, i, 'F')) + 1.0
                     for i in range(state.size)])


def _strength_matrix(Ground, Excited, Ez, Bz, pol,
                     initial='ground', initial_reduction='sum'):
    """(nG, nE) matrix whose cluster-sum is the line strength S, reduced over
    the chosen initial state.

    One primitive: S = sum_{M_g,M_e,p} |<g|T^1_p(d)|e>|^2  (thesis p.109),
    computed the same way in either basis:
      no-M  : branching_ratios = S/(2F_excited+1); rescaled by (2F_excited+1)
              -> S per (F'',F') level pair.
      M-res : sum_{p in pol} |Calculate_TDMs(p)|^2 per (M'',M') eigenstate
              pair; the M sums are realised by clustering downstream -> S.

    initial            : 'ground' | 'excited'
    initial_reduction  : 'sum'     -> S
                         'average' -> S / (2F_initial+1)
                             initial='excited' -> branching_ratios (emission)
                             initial='ground'  -> cross section (Eq.3.11)
    """
    if initial_reduction not in ('sum', 'average'):
        raise ValueError("initial_reduction must be 'sum' or 'average'")
    if initial not in ('ground', 'excited'):
        raise ValueError("initial must be 'ground' or 'excited'")
    eF1 = _twoF1(Excited)                                        # (nE,)
    if not _has_M(Ground):
        # branching_ratios[g,e] = S(g,e)/(2F'_e+1); undo -> S per level pair.
        S = branching_ratios(Ground, Excited, Ez, Bz) * eF1[None, :]
    else:
        S = None
        for pp in _POL[pol]:
            Tp = Calculate_TDMs(pp, Ground, Excited, Ez, Bz, q=[-1, 0, 1])  # (nE,nG)
            S = np.abs(Tp) ** 2 if S is None else S + np.abs(Tp) ** 2
        S = S.T                                                 # -> (nG, nE)
    if initial_reduction == 'sum':                              # line strength
        return S
    if initial == 'excited':                                    # / (2F'+1)
        return S / eF1[None, :]
    return S / _twoF1(Ground)[:, None]                          # initial='ground'


def line_list(Ground, Excited, g_idx, e_idx, *,
              origin=None, units='MHz', field=(0.0, 0.0),
              pol='all', weight='dipole2', thresh=1e-9,
              initial='ground', initial_reduction='sum'):
    """Tidy DataFrame of transitions between selected level sets.

    Parameters
    ----------
    Ground, Excited : MoleculeLevels (eigensystem run or run here via `field`)
    g_idx, e_idx     : index arrays (e.g. from state.select_q)
    origin           : band origin in cm^-1; default Excited.parameters['Origin']
    units            : 'MHz' (default) or 'cm-1' for the freq column
    field            : (Ez, Bz)
    pol              : 'all'|'z'|'+'|'-'|'x'  (sums |TDM|^2 over chosen p)
    weight           : 'dipole2' (bare |d|^2, in-band) or 'rate' (x nu^3)
    thresh           : drop rows with strength < thresh (absolute)
    initial          : 'ground' (default) | 'excited' -- the initial state
    initial_reduction: 'sum' (default) -> line strength S (LIF signal);
                       'average' -> S/(2F_initial+1):
                         initial='excited' -> emission branching ratio
                         initial='ground'  -> absorption cross section
                       (the strength column carries this AFTER clustering in
                       broaden(); in an M-resolved build each row is one
                       (M'',M') component and the cluster-sum realises S.)

    Returns a DataFrame with columns LINE_COLUMNS. Empty selection ->
    empty DataFrame (same columns) + UserWarning.
    """
    Ez, Bz = field
    Ground.eigensystem(Ez, Bz)
    Excited.eigensystem(Ez, Bz)
    if origin is None:
        origin = Excited.parameters['Origin']
    g_idx = np.asarray(g_idx, dtype=int)
    e_idx = np.asarray(e_idx, dtype=int)
    S = _strength_matrix(Ground, Excited, Ez, Bz, pol,
                         initial, initial_reduction)             # (nG, nE)

    rows = []
    for ig in g_idx:
        Eg = Ground.evals0[ig]
        for ie in e_idx:
            strength = float(S[ig, ie])
            if strength < thresh:
                continue
            nu_MHz = (Excited.evals0[ie] - Eg) + C_CM * origin
            if weight == 'rate':
                strength *= nu_MHz ** 3                          # relative; const dropped
            elif weight != 'dipole2':
                raise ValueError("weight must be 'dipole2' or 'rate'")
            freq = nu_MHz if units == 'MHz' else nu_MHz / C_CM
            rows.append({
                'freq': freq, 'strength': strength,
                'g_N': _dom(Ground, ig, 'N'), 'g_J': _dom(Ground, ig, 'J'),
                'g_F': _dom(Ground, ig, 'F'), 'g_par': _PARITY[Ground.parities[ig]],
                'e_J': _dom(Excited, ie, 'J'), 'e_F': _dom(Excited, ie, 'F'),
                'e_par': _PARITY[Excited.parities[ie]],
                'g_idx': int(ig), 'e_idx': int(ie)})
    if not rows:
        warnings.warn("line_list: no transitions above threshold for the "
                      "given selection; returning empty DataFrame.")
        return pd.DataFrame(columns=LINE_COLUMNS)
    return pd.DataFrame(rows, columns=LINE_COLUMNS).sort_values(
        'freq').reset_index(drop=True)


def averaged_branching(Ground, Excited, e_select, *,
                        field=(0.0, 0.0), normalize=True):
    """Branching ratio averaged over every excited M_F sublevel in `e_select`.

    Thesis Eq.3.4 orientation average, hyperfine-generalized: each excited
    level's normalized BR distribution is averaged with EQUAL weight per
    excited M_F sublevel. In an M-resolved build every M_F is its own
    eigenstate (weight 1 each); in a no-M build each selected eigenstate is
    a whole F' level carrying weight (2F'+1) (its M_F multiplicity), so F'
    enters (2F'+1)-weighted -- the (2J'+1)(2I+1)^-1 isotropic average, NOT
    equal-per-F'. branching_ratios columns are F'-independent in total
    (each ~Gamma), so the (2F'+1) must be applied explicitly here.

    Returns
    -------
    avg_df : DataFrame [g_N, g_J, g_F, g_par, BR, g_idx]  (sum BR == 1)
    per_df : DataFrame ['ground' label + one column per excited sublevel]
    """
    Ez, Bz = field
    BR = branching_ratios(Ground, Excited, Ez, Bz)              # (nG, nE)
    e_select = np.asarray(e_select, dtype=int)
    mres = _has_M(Excited)            # M-resolved: each M_F is its own eigenstate
    cols, labels, wts = [], [], []
    for ie in e_select:
        col = BR[:, ie].astype(float)
        tot = col.sum()
        if tot <= 0:
            raise ValueError(f"excited idx {ie} has zero total decay")
        cols.append(col / tot if normalize else col)
        Fp = float(_dom(Excited, ie, 'F'))
        wts.append(1.0 if mres else (2.0 * Fp + 1.0))           # M_F multiplicity
        labels.append("e[J=%s,F=%s,%s]" % (
            _dom(Excited, ie, 'J'), _dom(Excited, ie, 'F'),
            _PARITY[Excited.parities[ie]]))
    P = np.vstack(cols).T                                       # (nG, nSel)
    w = np.asarray(wts, float)
    w = w / w.sum()
    avg = P @ w                            # equal weight per M_F sublevel (Eq.3.4)

    g_all = np.arange(Ground.size)
    avg_df = pd.DataFrame({
        'g_N': [_dom(Ground, i, 'N') for i in g_all],
        'g_J': [_dom(Ground, i, 'J') for i in g_all],
        'g_F': [_dom(Ground, i, 'F') for i in g_all],
        'g_par': [_PARITY[Ground.parities[i]] for i in g_all],
        'BR': avg, 'g_idx': g_all})
    per_df = pd.DataFrame({'ground': [
        "X N=%s J=%s F=%s%s" % (_dom(Ground, i, 'N'), _dom(Ground, i, 'J'),
                                _dom(Ground, i, 'F'),
                                _PARITY[Ground.parities[i]]) for i in g_all]})
    for lab, col in zip(labels, cols):
        per_df[lab] = col
    return avg_df, per_df


_FWHM_SIG = 2.0 * np.sqrt(2.0 * np.log(2.0))     # FWHM = sigma * this (Gaussian)
_KB = 1.380649e-23
_AMU = 1.66053906660e-27
_C_SI = 299792458.0


def doppler_fwhm(T, mass_amu, nu0):
    """Doppler FWHM (thesis Eq.3.15). Returned in nu0's frequency unit."""
    sigma = (nu0 / _C_SI) * np.sqrt(_KB * T / (mass_amu * _AMU))
    return _FWHM_SIG * sigma


def _cluster(pos, s, thr):
    pos = np.asarray(pos, float); s = np.asarray(s, float)
    o = np.argsort(pos); pos, s = pos[o], s[o]
    if len(pos) == 0:
        return pos, s
    cut = np.where(np.diff(pos) >= thr)[0] + 1
    gp, gs = np.split(pos, cut), np.split(s, cut)
    return (np.array([p.mean() for p in gp]),
            np.array([q.sum() for q in gs]))


def broaden(lines, *, shape='gaussian', fwhm=None, lorentz_fwhm=None,
            T=None, mass_amu=None, nu0=None, x=None, norm='area',
            cluster=True, cluster_thresh=None, n_per_fwhm=20, pad=5):
    """Convolve a line list to a curve.

    shape : 'gaussian' | 'lorentzian' | 'voigt'
    fwhm  : explicit FWHM, OR omit and pass (T, mass_amu, nu0) for Doppler.
            For 'voigt', `fwhm` is the Gaussian (Doppler) FWHM and
            `lorentz_fwhm` the Lorentzian FWHM.
    norm  : 'area' (int g dnu = 1) | 'peak' (g_hat(0)=1)   (thesis p.111)
    cluster : merge lines closer than cluster_thresh (default fwhm/2).
    Returns (x, y_total, (component_pos, component_strength)).
    """
    if isinstance(lines, pd.DataFrame):
        pos = lines['freq'].to_numpy(float)
        s = lines['strength'].to_numpy(float)
    else:
        pos, s = np.asarray(lines[0], float), np.asarray(lines[1], float)
    if fwhm is None:
        if None in (T, mass_amu, nu0):
            raise ValueError("pass fwhm, or all of (T, mass_amu, nu0)")
        fwhm = doppler_fwhm(T, mass_amu, nu0)
    if len(pos) == 0:
        xx = np.array([]) if x is None else np.asarray(x, float)
        return xx, np.zeros_like(xx), (pos, s)
    if cluster:
        thr = fwhm / 2.0 if cluster_thresh is None else cluster_thresh
        pos, s = _cluster(pos, s, thr)
    if x is None:
        lo, hi = pos.min() - pad * fwhm, pos.max() + pad * fwhm
        n = max(int(np.ceil((hi - lo) / fwhm * n_per_fwhm)) + 1, 2)
        x = np.linspace(lo, hi, n)
    x = np.asarray(x, float)

    sig = fwhm / _FWHM_SIG
    if shape == 'gaussian':
        def k(d):
            g = np.exp(-0.5 * (d / sig) ** 2)
            return g / (sig * np.sqrt(2 * np.pi)) if norm == 'area' else g
    elif shape == 'lorentzian':
        hw = fwhm / 2.0
        def k(d):
            L = hw ** 2 / (d ** 2 + hw ** 2)             # peak = 1
            return (L / (np.pi * hw)) if norm == 'area' else L
    elif shape == 'voigt':
        from scipy.special import voigt_profile
        lf = fwhm if lorentz_fwhm is None else lorentz_fwhm
        gam = lf / 2.0
        v0 = voigt_profile(0.0, sig, gam)
        def k(d):
            v = voigt_profile(d, sig, gam)               # area-normalized
            return v if norm == 'area' else v / v0
    else:
        raise ValueError("shape must be gaussian|lorentzian|voigt")

    y = np.zeros_like(x)
    for p, a in zip(pos, s):
        y = y + a * k(x - p)
    return x, y, (pos, s)


def plot_spectrum(lines=None, *, Ground=None, Excited=None,
                  g_idx=None, e_idx=None, ax=None, sticks=True,
                  broaden_kw=None, normalize=True, experimental=None,
                  label='Simulation', **line_list_kw):
    """Plot sticks and/or a broadened curve, with an optional experiment.

    Forwards **line_list_kw to line_list (incl. initial / initial_reduction).
    experimental : dict(freq, signal, err=None, offset=0, tweak=0, yscale=1)
        generalizes the BaF_spectrum_plot.py offset/tweak/yscale overlay.
    """
    import matplotlib.pyplot as plt
    if lines is None:
        lines = line_list(Ground, Excited, g_idx, e_idx, **line_list_kw)
    if ax is None:
        ax = plt.gca()
    if len(lines):
        smax = lines['strength'].max()
        norm_s = (lambda v: v / smax) if (normalize and smax > 0) else (lambda v: v)
        if broaden_kw:
            x, y, _ = broaden(lines, **broaden_kw)
            if normalize and y.size and y.max() > 0:
                y = y / y.max()
            ax.plot(x, y, label=label)
        if sticks:
            ax.vlines(lines['freq'].to_numpy(),
                      0.0, norm_s(lines['strength'].to_numpy()),
                      color='0.4', alpha=0.7, linewidth=1)
    if experimental is not None:
        ex = experimental
        xd = np.asarray(ex['freq'], float) + ex.get('tweak', 0)
        sg = np.asarray(ex['signal'], float)
        yd = sg / sg.max() * ex.get('yscale', 1.0)
        if ex.get('err') is not None:
            er = np.asarray(ex['err'], float) / sg.max() * ex.get('yscale', 1.0)
            ax.errorbar(xd, yd, yerr=er, linestyle='none', marker='o',
                        label='Data')
        else:
            ax.plot(xd, yd, linestyle='none', marker='o', label='Data')
    ax.legend(loc='best')
    return ax
