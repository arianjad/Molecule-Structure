"""RaF A²Π eigenvalue-spacing regression: migrated-N² (via converter) vs the
explicit pre-migration post-ef06b86 R² constants. Spacings are Origin-free."""
import os, sys, numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from Energy_Levels import MoleculeLevels
from molecule_parameters import c, params_general

def _old_r2(spin):
    be0 = 5743.96 if spin == 'boson' else 5729.03
    d = {'Be': be0 - 2*1.4e-7*c, 'ASO': 1350*c, 'p+2q': -0.41071*c + 1.9e-7*c
         if spin == 'boson' else -0.4109*c + 1.9e-7*c, 'q': 0,
         'D': 1.4e-7*c, 'p2q_D': 1.9e-7*c, 'muE': 1*0.503412, 'g_S': 2.0023}
    if spin == 'boson':
        d.update({'h1/2': 0, 'a': 19/2, 'bF': 0, 'c': 0, 'd': -9,
                  'g_lp': -0.41071/(2*5743.96/c), 'g_l': -6000/(2*5743.96)})
    else:
        d.update({'h1/2Yb': -2852/2, 'dYb': -1*-0.076*c, 'h1/2H': 19/2,
                  'aH': 19/2, 'bFH': 0, 'cH': 0, 'dH': -9, 'e2Qq0': 0,
                  'g_lp': -0.4109/(2*0.191100)})
    md = dict(params_general); md.update(d)
    return md

def _evals(params, spin):
    # I_nuclei = [I_metal, I_ligand]. Boson RaF = even-A Ra (I_metal=0);
    # fermion RaF = 225Ra (I=1/2) + 19F (I=1/2), as RaF225 X-A.ipynb uses.
    # Omitting this defaults to [0,1/2], which feeds I_metal=0 into the
    # two-nuclear-spin _odd_aBJ path and yields an all-nan fermion spectrum.
    I_nuclei = [0, 1/2] if spin == 'boson' else [1/2, 1/2]
    m = MoleculeLevels.initialize_state(
        molecule_name='RaF', elec_state='A', vib_state='0',
        fermion_or_boson=spin, N_list=[1, 2, 3], round=8,
        params=params, P_values=[1/2, 3/2], I_nuclei=I_nuclei)
    return np.sort(m.eigensystem(0.0, 1e-9)[0])

if __name__ == "__main__":
    ok = True
    for spin in ('boson', 'fermion'):
        old = _evals(_old_r2(spin), spin)              # explicit R² (pre-migration)
        new = _evals(None, spin)                        # migrated N² via converter
        # Fail-loud on a non-finite spectrum: a nan/inf eigenvalue means the
        # build itself is broken (e.g. wrong I_nuclei), not a migration delta.
        # Surface that explicitly instead of letting it masquerade as a
        # spacing FAIL (nan < 1e-8 is False).
        for tag, ev in (("explicit-R²", old), ("migrated-N²", new)):
            bad = int((~np.isfinite(ev)).sum())
            if bad:
                print(f"  {spin}: {tag} spectrum has {bad}/{len(ev)} "
                      f"non-finite eigenvalues — build is broken (check "
                      f"I_nuclei), not a migration regression")
                print("FAIL build")
                sys.exit(2)
        n = min(len(old), len(new))
        d_old = np.diff(old[:n]); d_new = np.diff(new[:n])
        m = float(np.max(np.abs(d_old - d_new))) if n > 1 else 0.0
        print(f"  {spin}: max|Δspacing| = {m:.3e} MHz  ({n} levels)")
        ok = ok and (m < 1e-8)
    print("OK physics" if ok else "FAIL physics")
    sys.exit(0 if ok else 1)
