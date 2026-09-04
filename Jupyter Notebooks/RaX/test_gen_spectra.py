"""Hard-assert tests for gen_spectra (exit 0 iff all pass)."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'Source Code'))
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from Energy_Levels import MoleculeLevels, branching_ratios
import gen_spectra as xs

c = 29979.2458

def build(elec, N_list):
    return MoleculeLevels.initialize_state(
        molecule_name='BaF', elec_state=elec, vib_state=0,
        N_list=np.array(N_list), fermion_or_boson='boson',
        M_sublevels='none', I_nuclei=[0, 1/2], isotope=138,
        round=8, params=None, P_values=[1/2])

def build_M(elec, N_list):
    return MoleculeLevels.initialize_state(
        molecule_name='BaF', elec_state=elec, vib_state=0,
        N_list=np.array(N_list), fermion_or_boson='boson',
        M_sublevels='all', I_nuclei=[0, 1/2], isotope=138,
        round=8, params=None, P_values=[1/2])

g = build('X', [0, 1]); e = build('A', [1, 2])
g.eigensystem(0, 0); e.eigensystem(0, 0)

# Physics A: no-M emission strength == branching_ratios entry exactly
#   (initial='excited', initial_reduction='average' IS the branching observable)
BR = branching_ratios(g, e, 0, 0)                       # (nG, nE)
gx0 = g.select_q({'N': 0, 'J': 0.5})
exm = e.select_q({'J': 0.5}, parity='-')
ll = xs.line_list(g, e, gx0, exm, origin=e.parameters['Origin'],
                  initial='excited', initial_reduction='average')
for _, r in ll.iterrows():
    assert abs(r['strength'] - BR[int(r['g_idx']), int(r['e_idx'])]) < 1e-12, r
# Regression: the 3 Table II (-)/N=0 components within 2 MHz
TBL2 = {(1, 0): 348666402.6, (1, 1): 348666424.4, (0, 1): 348666490.0}
hits = 0
for _, r in ll.iterrows():
    key = (int(r['g_F']), int(r['e_F']))
    if key in TBL2:
        assert abs(r['freq'] - TBL2[key]) <= 2.0, (key, r['freq'], TBL2[key])
        hits += 1
assert hits == 3, hits

# ---- Tobs: modular observables, no-M <-> M='all' agreement ----
# Same molecule/params both ways (BaF DB, iso 138) -- only M_sublevels differs.
gN, eN = build('X', [0, 1]), build('A', [1, 2])
gMr, eMr = build_M('X', [0, 1]), build_M('A', [1, 2])
for _s in (gN, eN, gMr, eMr):
    _s.eigensystem(0, 0)
giN = gN.select_q({'N': 0}); eiN = eN.select_q({'J': 0.5}, parity='-')
giM = gMr.select_q({'N': 0}); eiM = eMr.select_q({'J': 0.5}, parity='-')
def _grp(L):
    return (L.groupby(['g_F', 'e_F'], as_index=False)
             .agg(s=('strength', 'sum'))
             .sort_values(['g_F', 'e_F']).reset_index(drop=True))
# line strength S (default sum): no-M == M='all' clustered
S_no = _grp(xs.line_list(gN, eN, giN, eiN, origin=eN.parameters['Origin']))
S_M = _grp(xs.line_list(gMr, eMr, giM, eiM, origin=eMr.parameters['Origin']))
assert np.allclose(S_no['s'].values, S_M['s'].values, rtol=1e-6, atol=1e-9), \
    np.c_[S_no['s'].values, S_M['s'].values]
# emission branching: no-M == M='all'
br_no = _grp(xs.line_list(gN, eN, giN, eiN, origin=eN.parameters['Origin'],
                          initial='excited', initial_reduction='average'))
br_M = _grp(xs.line_list(gMr, eMr, giM, eiM, origin=eMr.parameters['Origin'],
                         initial='excited', initial_reduction='average'))
assert np.allclose(br_no['s'].values, br_M['s'].values, rtol=1e-6, atol=1e-9)
# S == branching x (2F'+1)  (F' = excited dominant F)
Fp = S_no['e_F'].values * 2 + 1
assert np.allclose(S_no['s'].values, br_no['s'].values * Fp,
                   rtol=1e-6, atol=1e-9), \
    np.c_[S_no['s'].values, br_no['s'].values * Fp]
print("Tobs OK  S(no)==S(M); br(no)==br(M); S==br*(2F'+1)")
# Empty selection -> empty df + warning, no crash
import warnings
with warnings.catch_warnings(record=True) as w:
    warnings.simplefilter("always")
    empty = xs.line_list(g, e, np.array([], dtype=int), exm)
    assert len(empty) == 0 and list(empty.columns) == xs.LINE_COLUMNS
    assert any(issubclass(x.category, UserWarning) for x in w)
# ---- T3: averaged_branching == thesis Eq.3.4 generalization ----
exp = e.select_q({'J': 0.5}, parity='+')                # cooling excited state
avg_df, per_df = xs.averaged_branching(g, e, exp)
assert abs(avg_df['BR'].sum() - 1.0) < 1e-9, avg_df['BR'].sum()
# averaged BR == M_F-degeneracy-weighted mean of per-sublevel normalized cols
# (thesis Eq.3.4 orientation avg; no-M build -> weight (2F'+1) per F' eigenstate)
_BRm = branching_ratios(g, e, 0, 0)
_nc = np.array([_BRm[:, k] / _BRm[:, k].sum() for k in exp])
_Fv = np.array([float(xs._dom(e, k, 'F')) for k in exp])
_w = (2 * _Fv + 1); _w = _w / _w.sum()
_expected = _nc.T @ _w
assert np.allclose(avg_df['BR'].values, _expected, atol=1e-12), \
    np.abs(avg_df['BR'].values - _expected).max()
# guard against regression to the equal-per-F' bug
_plain = _nc.mean(axis=0)
assert not np.allclose(avg_df['BR'].values, _plain, atol=1e-6), \
    "averaged_branching regressed to equal-per-F' (must be (2F'+1)-weighted)"
# parity closure: X N=0 gets ~0
n0_mask = avg_df['g_N'].astype(float) == 0
assert avg_df.loc[n0_mask, 'BR'].sum() <= 1e-3, avg_df.loc[n0_mask, 'BR'].sum()
print("T3 OK  sum(BR)=%.6f  N0=%.1e" % (avg_df['BR'].sum(),
      avg_df.loc[n0_mask, 'BR'].sum()))
# ---- T4: broaden / doppler_fwhm / plot_spectrum ----
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fw = xs.doppler_fwhm(T=4.0, mass_amu=157.0, nu0=348.66e12)
ref = 2*np.sqrt(2*np.log(2)) * (348.66e12/299792458.0) * np.sqrt(
    1.380649e-23*4.0/(157.0*1.66053906660e-27))
assert abs(fw - ref)/ref < 1e-6, (fw, ref)

# area-normalized Gaussian: integral ~ total strength
lines = xs.line_list(g, e, g.select_q({'N': 1}), e.select_q({'J': 0.5}, parity='+'),
                     origin=e.parameters['Origin'])
x, y, (cp, cs) = xs.broaden(lines, shape='gaussian', fwhm=5.0,
                            norm='area', cluster=False)
area = np.trapz(y, x)
assert abs(area - lines['strength'].sum()) / lines['strength'].sum() < 0.01, area
# peak-normalized single line peaks at its strength
one = lines.iloc[[0]]
xx, yy, _ = xs.broaden(one, shape='gaussian', fwhm=5.0, norm='peak', cluster=False)
assert abs(yy.max() - one['strength'].iloc[0]) / one['strength'].iloc[0] < 1e-3
# clustering merges co-frequent lines
import pandas as pd
dup = pd.DataFrame({'freq': [100.0, 100.0, 300.0], 'strength': [1.0, 2.0, 4.0]})
_, _, (mp, ms) = xs.broaden(dup, shape='lorentzian', fwhm=10.0,
                            cluster=True, cluster_thresh=1.0)
assert len(mp) == 2 and abs(ms[0] - 3.0) < 1e-9 and abs(ms[1] - 4.0) < 1e-9, (mp, ms)
# plot_spectrum smoke (sticks + broadened + experimental)
fig, ax = plt.subplots()
expt = dict(freq=np.array([0.0, 50.0]), signal=np.array([1.0, 2.0]),
            err=np.array([0.1, 0.1]), offset=0, tweak=0, yscale=1.0)
ax2 = xs.plot_spectrum(lines=lines, ax=ax, sticks=True,
                       broaden_kw=dict(shape='voigt', fwhm=5.0, lorentz_fwhm=2.0),
                       experimental=expt, label='Sim')
assert ax2 is ax and len(ax.lines) >= 1
plt.close(fig)
print("T4 OK  doppler=%.3f MHz" % (fw/1e6))
print("T2 OK", len(ll), "rows; hits", hits)

# ---- Tclose: excited-state decay normalization over ALL E1 paths ----
#   Sum the emission observable over the COMPLETE E1-allowed ground manifold.
#   no-M : per-excited colsum == 1 (branching already normalized; thesis
#          Eq.3.6) -> a SUBSET of paths therefore sums to <1.
#   M=all: per-excited-M' colsum == 1/(2F'+1) (the average over the excited
#          (2F'+1)); the (2F'+1) M' of a given F' sum back to 1. Both bases
#          give unit closure after /colsum. Single vibronic band only
#          (vibrational/FCF branching is separate, thesis Eq.3.6).
def _cl(Msub):
    # X N=0..3 spans every E1-allowed final level of A(N=1, J=1/2&3/2):
    # truncating to N<=2 loses ~30% (-> sum 0.70) for the J=3/2->N=3 paths.
    gC = MoleculeLevels.initialize_state(molecule_name='BaF', elec_state='X',
        vib_state=0, N_list=np.array([0, 1, 2, 3]), fermion_or_boson='boson',
        M_sublevels=Msub, I_nuclei=[0, 1/2], isotope=138, round=8,
        params=None, P_values=[1/2])
    eC = MoleculeLevels.initialize_state(molecule_name='BaF', elec_state='A',
        vib_state=0, N_list=np.array([1]), fermion_or_boson='boson',
        M_sublevels=Msub, I_nuclei=[0, 1/2], isotope=138, round=8,
        params=None, P_values=[1/2, 3/2])
    gC.eigensystem(0, 0); eC.eigensystem(0, 0)
    Lc = xs.line_list(gC, eC, np.arange(gC.size), np.arange(eC.size),
                      origin=eC.parameters['Origin'], initial='excited',
                      initial_reduction='average', thresh=0.0)
    cs = Lc.groupby('e_idx')['strength'].sum()
    twoFp1 = np.array([2 * float(xs._dom(eC, k, 'F')) + 1 for k in cs.index])
    return cs.values, twoFp1
_csn, _ = _cl('none')
assert np.allclose(_csn, 1.0, atol=1e-6), ('no-M decay not normalized', _csn)
_csm, _twoFp1 = _cl('all')
assert np.allclose(_csm, 1.0 / _twoFp1, atol=1e-6), \
    ("M=all colsum != 1/(2F'+1)", np.c_[_csm, 1.0 / _twoFp1])
print("Tclose OK  no-M sum_paths=1 ; M=all sum_paths=1/(2F'+1)")

# ---- Tboltz: optional Boltzmann weighting on line strengths -----------
# boltzmann_T=None (default) returns the unweighted strength; boltzmann_T=T
# scales each row by exp(-E_g[MHz] * h/kB / T[K]).  Verify against an
# independent recompute and that boltzmann_T=None matches the default.
gB = build('X', [0, 1, 2]); eB = build('A', [1, 2])
gB.eigensystem(0, 0); eB.eigensystem(0, 0)
gB_idx = np.arange(gB.size); eB_idx = np.arange(eB.size)
L_none = xs.line_list(gB, eB, gB_idx, eB_idx,
                      origin=eB.parameters['Origin'], thresh=0.0)
T_test = 4.0
L_warm = xs.line_list(gB, eB, gB_idx, eB_idx,
                      origin=eB.parameters['Origin'], thresh=0.0,
                      boltzmann_T=T_test)
# Default of None must match boltzmann_T omitted
assert L_none['strength'].equals(
    xs.line_list(gB, eB, gB_idx, eB_idx,
                 origin=eB.parameters['Origin'], thresh=0.0,
                 boltzmann_T=None)['strength']), "boltzmann_T=None != default"
# Per-row check: L_warm.strength == L_none.strength * exp(-Eg*h/kB/T)
H_OVER_KB = 4.799243073e-5  # K/MHz
matched = 0
for (_, rn), (_, rw) in zip(L_none.iterrows(), L_warm.iterrows()):
    Eg = gB.evals0[int(rn['g_idx'])]
    expected = rn['strength'] * np.exp(-Eg * H_OVER_KB / T_test)
    assert abs(rw['strength'] - expected) < 1e-12 * max(abs(expected), 1e-30), \
        f"row {rn['g_idx']}: {rw['strength']} vs {expected}"
    matched += 1
assert matched > 0
# Negative T raises
try:
    xs.line_list(gB, eB, gB_idx, eB_idx, origin=eB.parameters['Origin'],
                 boltzmann_T=-1.0)
    raise AssertionError("boltzmann_T=-1.0 should have raised")
except ValueError:
    pass
print(f"Tboltz OK  T={T_test}K factor applied per-row ({matched} rows checked); "
      f"None=default; negative T raises")
