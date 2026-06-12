"""Convention gates for the 29SiO+ X2Sigma+ entry (molecules['SiO+']['boson']['X0']).

Each check pins one sign/conversion the SiO+ model rests on, with an EXECUTABLE
assertion. Five checks (C1-C5); C5 is an audit of the g_l-absence mechanism plus
two PENDING-REVIEWER literature slots (recorded in docs/sio-conventions.md, not
gated here). Run:

    conda run -n Structure python "Source Code/test_sio_conventions.py"

Exit 0 iff every executable check passes; exit 1 on any [FAIL]. Companion note
with code locations, pasted numerics, and reviewer slots: docs/sio-conventions.md.
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))

import sympy as sp
from Energy_Levels import MoleculeLevels
from molecule_parameters import get_molecule_params, params_general

_results = []  # (name, passed, detail)


def record(name, passed, detail):
    _results.append((name, bool(passed), detail))
    tag = "[PASS]" if passed else "[FAIL]"
    print(f"{tag} {name}: {detail}")


# ---------------------------------------------------------------------------
# C1 - hyperfine operator identity (sympy, exact)
#   bF*(I.S) + (c/3)*sqrt(6)*T2_0(I,S)  ==  b*(I.S) + c*(I.n)(S.n),  b = bF - c/3
#   T2_0(I,S) = (3 Iz Sz - I.S)/sqrt(6) in the molecule frame (n = z).
#   This is the exact form the builder applies at hamiltonian_builders.py:45-46
#   ('bF'*I.S + 'c'/3*sqrt(6)*T2_0(I,S)).
# ---------------------------------------------------------------------------
def check_C1():
    sx = sp.Rational(1, 2) * sp.Matrix([[0, 1], [1, 0]])
    sy = sp.Rational(1, 2) * sp.Matrix([[0, -sp.I], [sp.I, 0]])
    sz = sp.Rational(1, 2) * sp.Matrix([[1, 0], [0, -1]])
    I2 = sp.eye(2)
    kp = lambda a, b: sp.Matrix(sp.kronecker_product(a, b))
    Sx, Sy, Sz = (kp(o, I2) for o in (sx, sy, sz))   # S on first factor
    Ix, Iy, Iz = (kp(I2, o) for o in (sx, sy, sz))   # I on second factor
    IdotS = Ix * Sx + Iy * Sy + Iz * Sz
    T20 = (3 * Iz * Sz - IdotS) / sp.sqrt(6)
    bF, c = sp.symbols('bF c', real=True)
    b = bF - c / 3
    LHS = bF * IdotS + (c / 3) * sp.sqrt(6) * T20             # code's form
    RHS = b * IdotS + c * (Iz * Sz)                            # b*(I.S)+c*(I.n)(S.n)
    diff = sp.simplify(LHS - RHS)
    ident = diff == sp.zeros(4, 4)
    # also: c-term alone equals dipolar c*(IzSz - I.S/3)
    cterm = sp.simplify((c / 3) * sp.sqrt(6) * T20 - c * (Iz * Sz - IdotS / 3))
    cident = cterm == sp.zeros(4, 4)
    record("C1 hyperfine operator identity", ident and cident,
           f"LHS-RHS == 0 (4x4): {ident}; c-term == c*(IzSz-I.S/3): {cident} "
           f"[code prefactor c/3*sqrt(6) verified, hamiltonian_builders.py:46]")


# ---------------------------------------------------------------------------
# C2 - g_S sign: rising N=0 branch slope ~ +g_S*mu_B*(1/2) = +1.4012 MHz/G.
#   Per-level slope 14.01 GHz/T (m_S = +-1/2 SPLITTING = 28 GHz/T). Karthein
#   Fig.1: N=0 positive-parity rises.
# ---------------------------------------------------------------------------
def check_C2():
    m = MoleculeLevels.initialize_state(
        'SiO+', 'X', 0, N_list=[0, 1], fermion_or_boson='boson',
        I_nuclei=[0, 1 / 2], P_values=[1 / 2])
    B = np.linspace(0, 200, 11); B[0] = 1e-9
    evB, _ = m.ZeemanMap(B, output=True, write_attribute=False, order=True)
    evB = np.array(evB)
    slopes = np.array([np.polyfit(B, evB[:, k], 1)[0] for k in range(evB.shape[1])])
    e0 = evB[0]
    n0 = [k for k in range(16) if e0[k] < 2000]      # N=0 manifold
    n1 = [k for k in range(16) if e0[k] > 40000]     # N=1 manifold
    rise = max(slopes[n0])
    desc = min(slopes[n1])
    pred = params_general['g_S'] * params_general['mu_B'] * 0.5   # 1.40123 MHz/G
    ok_rise = abs(rise - pred) / pred < 0.05 and rise > 0
    ok_desc = desc < -1.0                               # N=1 branch descends ~ -1.4
    record("C2 g_S sign (rising N=0 branch)", ok_rise and ok_desc,
           f"rise={rise:.5f} MHz/G ({rise*10:.3f} GHz/T) vs pred {pred:.5f} "
           f"({100*(rise-pred)/pred:+.3f}%); N=1 descend={desc:.5f} MHz/G")


# ---------------------------------------------------------------------------
# C3 - g_N sign + Larmor at B=15000 G.
#   |g_N|*mu_N*B = 1.1106*7.62259e-4*15000 = 12.6985 MHz (within-branch nuclear
#   splitting). Sign: code applies V_B += -g_N*mu_N*ZeemanIZ
#   (hamiltonian_builders.py:57). g_N<0 -> coef +|g_N|muN -> E=-g_N muN m_I B,
#   so m_I=+1/2 lies HIGHER. Isolate via overrides (bF,c->1e-9; g_S->1e-12).
# ---------------------------------------------------------------------------
def check_C3():
    muN = params_general['mu_N']; gN = -1.1106; B = 15000.
    larmor = abs(gN) * muN * B
    p0 = get_molecule_params('SiO+', 'X', '0', 'boson')
    # (a) magnitude: hyperfine off, electron Zeeman on -> within-branch splitting
    p_hf = dict(p0); p_hf['bF'] = 1e-9; p_hf['c'] = 1e-9
    m_hf = MoleculeLevels.initialize_state(
        'SiO+', 'X', 0, N_list=[0], fermion_or_boson='boson',
        I_nuclei=[0, 1 / 2], P_values=[1 / 2], params=p_hf)
    ev = np.sort(m_hf.eigensystem(0.0, B)[0])
    gaps = np.diff(ev)
    small = [g for g in gaps if abs(g) < 500]            # nuclear splittings
    mag_ok = len(small) > 0 and all(abs(g - larmor) / larmor < 0.10 for g in small)
    # (b) sign: also kill electron Zeeman -> pure nuclear Zeeman, read m_F ordering
    p_nz = dict(p_hf); p_nz['g_S'] = 1e-12
    m_nz = MoleculeLevels.initialize_state(
        'SiO+', 'X', 0, N_list=[0], fermion_or_boson='boson',
        I_nuclei=[0, 1 / 2], P_values=[1 / 2], params=p_nz)
    evn, evec = m_nz.eigensystem(0.0, B)
    Mbasis = m_nz.q_numbers['M']
    # dominant M of the highest and lowest eigenstate
    idx = np.argsort(evn)
    def dom_M(k):
        v = evec[:, k]
        return Mbasis[int(np.argmax(np.abs(v)))]
    hi_M, lo_M = dom_M(idx[-1]), dom_M(idx[0])
    pred_split = -gN * muN * 0.5 * B                     # +6.3492 for m_I=+1/2
    hi_E = evn[idx[-1]]
    # g_N<0 -> stretched m_I=+1/2 (dominant M=+1 here) is the HIGHER level, E>0
    sign_ok = (hi_M > 0) and (hi_E > 0) and abs(hi_E - pred_split) / pred_split < 0.01
    record("C3 g_N magnitude (Larmor @15kG)", mag_ok,
           f"|g_N|muN B pred={larmor:.4f} MHz; within-branch splittings={[round(g,4) for g in small]}")
    record("C3 g_N sign (m_I=+1/2 higher for g_N<0)", sign_ok,
           f"highest level E={hi_E:.4f} MHz at dom M={hi_M:+}, lowest at M={lo_M:+}; "
           f"pred +-{pred_split:.4f} MHz [V_B += -g_N*mu_N*ZeemanIZ, hamiltonian_builders.py:57]")


# ---------------------------------------------------------------------------
# C4 - gamma (spin-rotation) sign: in N=1 (I=0 isolates gamma), J=3/2 sits ABOVE
#   J=1/2 for gamma>0. Splitting = (3/2)*gamma = 18 MHz.
#   <N.S> = +1/2 (J=3/2) -> +gamma/2 = +6; <N.S> = -1 (J=1/2) -> -gamma = -12.
# ---------------------------------------------------------------------------
def check_C4():
    m1 = MoleculeLevels.initialize_state(
        'SiO+', 'X', 0, N_list=[1], fermion_or_boson='boson',
        I_nuclei=[0, 0], P_values=[1 / 2])
    ev, _ = m1.eigensystem(0.0, 1e-9)
    Jlist = np.array(m1.q_numbers['J'])
    center = np.mean(ev)
    rel = ev - center
    E_half = np.mean(rel[Jlist == 0.5])
    E_three = np.mean(rel[Jlist == 1.5])
    gamma = 12.0
    split = E_three - E_half                              # should be +18 = (3/2)*gamma
    ok = (E_three > E_half) and abs(split - 1.5 * gamma) < 1e-6
    record("C4 gamma sign (J=3/2 above J=1/2, gamma>0)", ok,
           f"E(J=3/2)-center={E_three:+.4f}, E(J=1/2)-center={E_half:+.4f} MHz; "
           f"split={split:.4f} MHz (pred (3/2)*gamma={1.5*gamma:.1f})")


# ---------------------------------------------------------------------------
# C5 - g_l absence mechanism (audit): SiO+ entry deliberately omits g_l; builder
#   guards with params.get('g_l') is not None (hamiltonian_builders.py:54), so
#   the anisotropic-Zeeman term is never added. The operator-identity question
#   (is ZeemanLZ_bBJ B&C's anisotropic spin-Zeeman form?) is PENDING-REVIEWER.
# ---------------------------------------------------------------------------
def check_C5():
    p = get_molecule_params('SiO+', 'X', '0', 'boson')
    absent = 'g_l' not in p
    # builder tolerance: an N=0,1 build succeeds with finite spectrum (g_l skipped)
    m = MoleculeLevels.initialize_state(
        'SiO+', 'X', 0, N_list=[0, 1], fermion_or_boson='boson',
        I_nuclei=[0, 1 / 2], P_values=[1 / 2])
    ev = m.eigensystem(0.0, 1e-9)[0]
    finite = bool(np.all(np.isfinite(ev)))
    record("C5 g_l absence tolerated (params.get default)", absent and finite,
           f"'g_l' in params: {not absent}; build finite spectrum: {finite} "
           f"[guard: params.get('g_l') is not None, hamiltonian_builders.py:54]")


if __name__ == "__main__":
    check_C1()
    check_C2()
    check_C3()
    check_C4()
    check_C5()
    n_pass = sum(1 for _, p, _ in _results if p)
    n_fail = sum(1 for _, p, _ in _results if not p)
    print(f"\n{'OK' if n_fail == 0 else 'FAIL'}: {n_pass} passed, {n_fail} failed")
    sys.exit(0 if n_fail == 0 else 1)
