"""Molecule-agnostic X-A spectrum / branching toolkit.

Generalizes the copy-pasted simulate_spectra / simulate_spectra_noM /
plot_gaussian_spectrum helpers (RaF X-A.ipynb, BaF X-A.ipynb,
BaF_spectrum_plot.py) into one physics-correct module.

Physics (thesis §3.2.4.1, pp.107-115; spec
docs/superpowers/specs/2026-05-16-xa-spectra-toolkit-design.md §3):
  line strength S = sum_{M'',M',p} |<g,i|T^1_p(d)|e,j>|^2   (thesis p.109)
  averaged BR  r~ = |<g,J''||d||e,J'>|^2 / [(2J'+1) D_ge^2]  (thesis Eq.3.4)
omega_0^3 is dropped (cancels in any branching ratio; <1e-5 across one band).
NO Boltzmann population weight, NO ad-hoc J_adjust.

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


def _strength_matrix(Ground, Excited, Ez, Bz, pol):
    """(nG, nE) line-strength matrix.

    no-M  : branching_ratios entry = |<g,J''|d|e,J'>|^2.
    M-res : S[g,e] = sum_{p in pol} |Calculate_TDMs(p)|^2  (thesis p.109 S,
            the M'' / M' sum is realized by clustering the degenerate
            M-components downstream in broaden()).
    """
    if not _has_M(Ground):
        return branching_ratios(Ground, Excited, Ez, Bz)        # (nG, nE)
    S = None
    for pp in _POL[pol]:
        Tp = Calculate_TDMs(pp, Ground, Excited, Ez, Bz, q=[-1, 0, 1])  # (nE,nG)
        S = np.abs(Tp) ** 2 if S is None else S + np.abs(Tp) ** 2
    return S.T                                                  # -> (nG, nE)


def line_list(Ground, Excited, g_idx, e_idx, *,
              origin=None, units='MHz', field=(0.0, 0.0),
              pol='all', weight='dipole2', thresh=1e-9):
    """Tidy DataFrame of X-A transitions between selected level sets.

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
    S = _strength_matrix(Ground, Excited, Ez, Bz, pol)          # (nG, nE)

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
