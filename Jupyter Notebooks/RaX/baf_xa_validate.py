"""BaF X-A validation (spec §6; scope refined per user 2026-05-16: cooling line).

Frequencies: validated on the (-)-parity line A2Pi1/2(v=0,J=1/2,-) <- X2Sig+(v=0,
  N=0) vs arXiv:2511.06986 Table II (the paper's precision-measured components).
Branching:  the LASER-COOLING transition A2Pi1/2(v=0,J=1/2,+) -> X2Sig+(v=0,N=1),
  rotationally closed by parity -- A(+) decays E1 only to -parity X (odd N=1;
  N=0/N=2 are +parity, forbidden). This (not the (-)->N=0,1 of the old spec) is
  the physically meaningful cooling BR.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'Source Code'))
import numpy as np
from Energy_Levels import MoleculeLevels, branching_ratios

c = 29979.2458  # MHz per cm^-1 (== molecule_parameters.c)

def build(elec, N_list):
    return MoleculeLevels.initialize_state(
        molecule_name='BaF', elec_state=elec, vib_state=0,
        N_list=np.array(N_list), fermion_or_boson='boson',
        M_sublevels='none', I_nuclei=[0, 1/2], isotope=138,
        round=8, params=None, P_values=[1/2])

g = build('X', [0, 1])     # X2Sig+ v=0: N=0 (+par), N=1 (-par)
e = build('A', [1, 2])     # A2Pi1/2 v=0: J=1/2 (both Lambda-parities), J=3/2
g.eigensystem(0, 0)
e.eigensystem(0, 0)

def dom(s, i, q):
    """Dominant basis label of quantum q for eigenstate i (codebase convention)."""
    return s.q_numbers[q][np.argmax(s.evecs0[i] ** 2)]

# ---- (1) Hyperfine-splitting gates (mapping physics; Origin-independent) ----
gx0 = g.select_q({'N': 0, 'J': 0.5})                  # X N=0,J=1/2  (F=0,1)
x_split = float(abs(np.diff(np.sort(g.evals0[gx0]))[-1]))
exm = e.select_q({'J': 0.5}, parity='-')              # A 2Pi1/2 J=1/2 (-)  (F=0,1)
a_split = float(abs(np.diff(np.sort(e.evals0[exm]))[-1]))
print(f"X  N=0,J=1/2 F-split       = {x_split:8.3f} MHz  (Table II 65.6, tol 1.0)")
print(f"A  2Pi1/2 J=1/2(-) F-split = {a_split:8.3f} MHz  (Table II 21.8, tol 1.0)")

# ---- (2) Absolute freq: the 3 Table II (-)-line components, matched F_X->F_A ---
Org = e.parameters['Origin']
TBL2 = {(1, 0): 348666402.6, (1, 1): 348666424.4, (0, 1): 348666490.0}
freq_res = []
print("\n(-)-line vs Table II:  F_X->F_A     computed (MHz)      Table II       d")
for jx in gx0:
    Fx = int(dom(g, jx, 'F'))
    for ia in exm:
        Fa = int(dom(e, ia, 'F'))
        nu = (e.evals0[ia] - g.evals0[jx]) + c * Org
        ref = TBL2.get((Fx, Fa))
        if ref is None:
            continue                                   # F0->F0 is E1-forbidden
        freq_res.append(nu - ref)
        print(f"                        F{Fx}->F{Fa}   {nu:16.1f}  {ref:14.1f}  {nu-ref:+6.1f}")
max_freq_res = max(abs(d) for d in freq_res)

# ---- (3) Cooling-line BR: A2Pi1/2(J=1/2,+) -> X N=1 (parity-closed) ----
exp = e.select_q({'J': 0.5}, parity='+')               # cooling excited state
gN = np.array([dom(g, i, 'N') for i in range(g.size)]) # N per X eigenstate
BR = branching_ratios(Ground=g, Excited=e, Ez=0, Bz=0)
assert np.all(np.isfinite(BR)) and np.all(BR >= 0), "BR not finite/non-negative"
n1_idx = [i for i in range(g.size) if dom(g, i, 'N') == 1]
print("\ncooling line  A2Pi1/2(v=0,J=1/2,+) -> X2Sig+(v=0,N=1)  (parity-closed):")
hdr = "  ".join(f"N1[J={dom(g,i,'J')},F={int(dom(g,i,'F'))}]" for i in n1_idx)
print(f"  A(+) sublevel        to-X-N1  to-X-N0     {hdr}")
br_ok = True
for col in exp:
    tot = BR[:, col].sum()
    assert tot > 0, f"A(+) idx {col} has zero total decay"
    p = BR[:, col] / tot
    to_N0 = float(p[gN == 0].sum())
    to_N1 = float(p[gN == 1].sum())
    per = "  ".join(f"{p[i]:13.4f}" for i in n1_idx)
    npop = int(np.sum([p[i] > 1e-6 for i in n1_idx]))
    print(f"  idx {col} (F={int(dom(e,col,'F'))})  {to_N1:9.4f}  {to_N0:8.1e}     {per}")
    br_ok &= abs(p.sum() - 1.0) <= 1e-6     # normalized
    br_ok &= to_N1 >= 0.999                 # rotational closure into N=1
    br_ok &= to_N0 <= 1e-3                  # N=0 parity-forbidden (the defining property)
    br_ok &= npop >= 2                      # spreads over the N=1 hyperfine manifold

# ---- verdict ----
ok = True
ok &= abs(x_split - 65.6) <= 1.0            # X Fermi-contact mapping
ok &= abs(a_split - 21.8) <= 1.0            # A 2Pi1/2 hyperfine mapping (spec §3.3)
ok &= max_freq_res <= 2.0                   # paper abs acc ~1 MHz; fixed (non-refit) lit constants
ok &= br_ok                                 # parity-closed cooling BR
print(f"\nmax |freq residual| = {max_freq_res:.1f} MHz (tol 2.0; paper abs acc ~1 MHz)")
print("VALIDATION PASSED" if ok else "VALIDATION FAILED")
sys.exit(0 if ok else 1)
