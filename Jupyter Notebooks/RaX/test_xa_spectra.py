"""Hard-assert tests for xa_spectra (exit 0 iff all pass)."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'Source Code'))
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from Energy_Levels import MoleculeLevels, branching_ratios
import xa_spectra as xs

c = 29979.2458

def build(elec, N_list):
    return MoleculeLevels.initialize_state(
        molecule_name='BaF', elec_state=elec, vib_state=0,
        N_list=np.array(N_list), fermion_or_boson='boson',
        M_sublevels='none', I_nuclei=[0, 1/2], isotope=138,
        round=8, params=None, P_values=[1/2])

g = build('X', [0, 1]); e = build('A', [1, 2])
g.eigensystem(0, 0); e.eigensystem(0, 0)

# Physics A: no-M strength == branching_ratios entry exactly
BR = branching_ratios(g, e, 0, 0)                       # (nG, nE)
gx0 = g.select_q({'N': 0, 'J': 0.5})
exm = e.select_q({'J': 0.5}, parity='-')
ll = xs.line_list(g, e, gx0, exm, origin=e.parameters['Origin'])
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
# Empty selection -> empty df + warning, no crash
import warnings
with warnings.catch_warnings(record=True) as w:
    warnings.simplefilter("always")
    empty = xs.line_list(g, e, np.array([], dtype=int), exm)
    assert len(empty) == 0 and list(empty.columns) == xs.LINE_COLUMNS
    assert any(issubclass(x.category, UserWarning) for x in w)
print("T2 OK", len(ll), "rows; hits", hits)
