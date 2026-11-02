"""Stark-X (p=+-1) + X-X transition-dipole (TDM) gates for 29SiO+ X2Sigma+.

Companion to test_sio_crossing.py / test_sio_conventions.py. Validates the
lab-frame TRANSVERSE electric-dipole element StarkX_bBJ (matrix_elements.py) and
the hyperfine-resolved X-X transition-dipole matrices d_p (p=-1,0,+1) between the
N=0 and N=1 manifolds.

    conda run -n Structure python "Source Code/test_sio_stark_tdm.py"

Exit 0 iff every gate passes; exit 1 on any [FAIL].

Registration route (deliverable #2, recorded here for the downstream consumers).
    H_even_X (hamiltonian_builders.py) has NO transverse-field hook -- neither E
    nor B: it wires only the axial p=0 operators (ZeemanZ/StarkZ/ZeemanLZ/ZeemanIZ)
    into V_B/V_E and exposes H_function(E,B) over a scalar (z) E and B only.
    Retrofitting a transverse field-map into that shared builder would change the
    H_function signature and risk every aBJ/bBS path + the tutorial regression,
    and the downstream consumers (master-equation export, two-level driver) need
    the operator MATRICES, not a field-map API. So StarkX (and the pre-existing,
    never-registered ZeemanX) are registered in the bBJ_even_X ext_fields dict
    (molecule_library_class.py) -- exposed through the SAME element-name -> matrix
    -element machinery as StarkZ -- and the d_p / d_x matrices are obtained as
    retrievable matrices via hamiltonian_builders.build_operator over the model's
    own q_numbers basis (the build_d helper below). This is the same operator the
    spectrum uses: -muE*StarkZ matrix == H_function(1,0)-H_function(0,0) to 0.0
    (gate G2 re-checks this), so d_z is THE trusted axial dipole, and d_x is built
    from StarkZ's reduced parts by the codebase's own axial->transverse transform.

X-X TDM access (deliverable #3). The codebase's TDM registration paths
(collect_p_TDM / collect_TDM in molecule_library_class.py) are keyed by excited
-state labels and exist only for the case-a (aBJ) A-X channel; there is no bBJ
X-X TDM entry for any molecule. For a 2Sigma+ ground state the X-X electric-dipole
matrices d_p ARE the Stark operators: d_0 == StarkZ matrix exactly, and
d_{+-1} == the p=+-1 spherical components whose real combination (d_{-1}-d_{+1})
/sqrt2 == StarkX. The build_d helper below returns those d_p matrices (in units of
muE; multiply by muE for MHz/(V/cm)) -- they are the deliverable.

PHASE note. StarkX_bBJ follows the codebase's existing transverse-operator phase
convention exactly (ZeemanZ_bBJ -> ZeemanX_bBJ): the (-1)**(p/2+1/2) weighting
(+T1_{-1} for p=-1, -T1_{+1} for p=+1) and the 1/sqrt2. This coincides with the
canonical d_x = (T1_{-1}-T1_{+1})/sqrt2. No GATE depends on the undetermined
overall phase: G1 uses |element|^2, G2 uses |<A|d|B>|, G4 uses Hermiticity and the
parity-block structure -- all overall-sign-free.
"""
import os
import sys

import numpy as np
from functools import partial

sys.path.insert(0, os.path.dirname(__file__))

from Energy_Levels import MoleculeLevels
import hamiltonian_builders as ham
import matrix_elements as me
from sympy.physics.wigner import wigner_3j, wigner_6j

E_DRIVE = 6.0                 # V/cm, Karthein axial drive field
ORACLE_KHZ = 3.3319           # kHz, trusted crossing-drive value (test_sio_crossing G5)

_results = []


def record(name, passed, detail):
    _results.append((name, bool(passed), detail))
    tag = "[PASS]" if passed else "[FAIL]"
    print(f"{tag} {name}: {detail}")


def _kron(a, b):
    return 1 if a == b else 0


# ---------------------------------------------------------------------------
# d_p builder. Raw spherical dipole component T^1_p(d), unit muE (the
# direction-cosine reduced element). Same reduced parts as StarkZ_bBJ; only the
# lab F,M Wigner factor carries the projection p. d_p[0] reproduces StarkZ exactly
# (checked in G1), and (d_p[-1]-d_p[+1])/sqrt2 reproduces StarkX exactly (G1).
# ---------------------------------------------------------------------------
def _dipole_p_me(p, K0, N0, J0, F0, M0, K1, N1, J1, F1, M1, S=1 / 2, I=1 / 2):
    if not _kron(K0, K1):
        return 0.0
    if int(round(M0 - M1)) != int(round(p)):
        return 0.0
    return (-1) ** (F0 - M0 + J0 + I + F1 + 1 + N0 + S + J1 + 1 + N0 - K0) * \
        np.sqrt((2 * F0 + 1) * (2 * F1 + 1) * (2 * J0 + 1) * (2 * J1 + 1) *
                (2 * N0 + 1) * (2 * N1 + 1)) * \
        float(wigner_3j(F0, 1, F1, -M0, p, M1)) * \
        float(wigner_6j(J1, F1, I, F0, J0, 1)) * \
        float(wigner_6j(N1, J1, S, J0, N0, 1)) * \
        float(wigner_3j(N0, 1, N1, -K0, 0, K1))


def build_d(M):
    """Return the hyperfine-resolved X-X dipole matrices over M's own basis.

    Returns dict with keys:
        'z'  : d_z matrix      == StarkZ matrix (the trusted axial dipole),
        'x'  : d_x matrix      == StarkX matrix (transverse, p=+-1),
        -1, 0, +1 : d_p spherical-component matrices (unit muE).
    Multiply by M.parameters['muE'] for MHz/(V/cm)."""
    qn = M.q_numbers
    out = {
        'z': np.array(ham.build_operator(qn, qn, M.matrix_elements['StarkZ'])),
        'x': np.array(ham.build_operator(qn, qn, M.matrix_elements['StarkX'])),
    }
    for p in (-1, 0, 1):
        out[p] = np.array(ham.build_operator(qn, qn, partial(_dipole_p_me, p)))
    return out


def build_model(N_list):
    return MoleculeLevels.initialize_state(
        'SiO+', 'X', 0, N_list=N_list, fermion_or_boson='boson',
        I_nuclei=[0, 1 / 2], P_values=[1 / 2])


# ===========================================================================
# G1 -- rotational invariance / closed-shell N=0->1 sum rule.
#   For each zero-field N=0 eigenstate |i>,
#       S_i = sum_{f in N=1} sum_{p in {-1,0,+1}} |<f|d_p|i>|^2
#   must be IDENTICAL across all M-sublevels and across F=0/F=1 (m- and
#   F-independence = the rotational-invariance test), and equal muE^2 * 1 (the
#   2Sigma N=0->1 line-strength factor, verified numerically not asserted: d_p
#   built unit-muE => S_i should be 1). Consistency of the new p=+-1 elements with
#   the trusted p=0 element: S_i computed with ONLY p=0 (over all M') must equal
#   1/3 of the total (isotropy). Also checks d_p[0] == StarkZ and
#   (d_p[-1]-d_p[+1])/sqrt2 == StarkX exactly.
# ===========================================================================
def gate1_sum_rule():
    M = build_model([0, 1])
    qn = M.q_numbers
    d = build_d(M)
    Nq = np.array(qn['N'])
    n1 = np.where(Nq == 1)[0]              # N=1 final-manifold basis indices

    # Internal consistency of the new elements vs the trusted p=0 / StarkX.
    dev_dz = float(np.max(np.abs(d[0] - d['z'])))
    dev_dx = float(np.max(np.abs(d['x'] - (d[-1] - d[1]) / np.sqrt(2))))

    ev, vec = M.eigensystem(0.0, 1e-9)     # vec[k] = kth eigenvector (rows)
    domN = np.array([Nq[int(np.argmax(vec[k] ** 2))] for k in range(M.size)])
    i0 = np.where(domN == 0)[0]

    def Si(vi, ps):
        tot = 0.0
        for p in ps:
            col = d[p] @ vi                 # <f|d_p|i> over all basis f
            tot += float(np.sum(col[n1] ** 2))
        return tot

    rows = []
    for k in i0:
        vi = vec[k]
        idx = int(np.argmax(vi ** 2))
        rows.append((qn['F'][idx], qn['J'][idx], qn['M'][idx],
                     Si(vi, (-1, 0, 1)), Si(vi, (0,))))

    S_all = np.array([r[3] for r in rows])
    S_p0 = np.array([r[4] for r in rows])
    spread = float(S_all.max() - S_all.min())          # m- & F-independence
    frac = S_p0 / S_all                                # isotropy: must be 1/3
    muE2 = M.parameters['muE'] ** 2

    print("\n  --- G1 sum rule (d_p built unit-muE; S in muE^2 units) ---")
    print(f"  {'F':>2} {'J':>4} {'M':>3}  {'S_total':>12} {'S_p0':>12} {'S_p0/S':>9}")
    for F, J, Mm, Sa, S0 in sorted(rows):
        print(f"  {F:2.0f} {J:4.1f} {Mm:+3.0f}  {Sa:12.8f} {S0:12.8f} {S0 / Sa:9.6f}")
    print(f"  S_total spread (max-min) = {spread:.2e}  (m- & F-INDEPENDENCE)")
    print(f"  S_total mean = {S_all.mean():.10f}  -> 2Sigma N=0->1 line-strength "
          f"factor = {S_all.mean():.6f} (x muE^2 = {S_all.mean() * muE2:.6f} MHz^2/(V/cm)^2)")
    print(f"  d_p[0]==StarkZ dev={dev_dz:.2e} ; (d_-1 - d_+1)/sqrt2==StarkX dev={dev_dx:.2e}")

    ok = (spread < 1e-12 and                              # m- & F-independent
          np.all(np.abs(frac - 1 / 3) < 1e-12) and        # isotropy (p=0 = 1/3)
          abs(S_all.mean() - 1.0) < 1e-9 and              # factor = muE^2 * 1
          dev_dz < 1e-12 and dev_dx < 1e-12)              # element consistency
    record("G1 rotational invariance (m- & F-independent sum rule = muE^2*1; "
           "p=0 = 1/3 by isotropy)", ok,
           f"S_total={S_all.mean():.8f} (spread {spread:.1e}, 4 N=0 states); "
           f"S_p0/S_total={frac.mean():.6f} (target 1/3); "
           f"d_p[0]=StarkZ ({dev_dz:.0e}), StarkX combo ({dev_dx:.0e})")


# ===========================================================================
# G2 -- cross-consistency with the trusted crossing oracle. Re-derive the
#   N=0(+)/N=1(-) drive at the Karthein crossing from the d_z matrix (V_E =
#   -muE*E*d_z) and reproduce 3.332 kHz at 6 V/cm. First checks -muE*d_z equals
#   the builder's own field operator H(1,0)-H(0,0) exactly, tying the new d_p
#   construction to the operator the spectrum uses.
# ===========================================================================
def gate2_cross_consistency():
    M = build_model([0, 1, 2])
    d = build_d(M)
    muE = M.parameters['muE']
    VE_mine = -muE * d['z']
    VE_builder = M.H_function(1.0, 1e-9) - M.H_function(0.0, 1e-9)
    dev_op = float(np.max(np.abs(VE_mine - VE_builder)))

    def zmap(Bv):
        ev, vec = M.ZeemanMap(Bv, output=True, write_attribute=False, order=True)
        return np.array(ev), np.array(vec)

    def domN(v):
        return M.q_numbers['N'][int(np.argmax(v ** 2))]

    def domM(v):
        return M.q_numbers['M'][int(np.argmax(v ** 2))]

    def par(v):
        return int(round(v @ M.Parity_mat @ v))

    B = np.linspace(1e-6, 16000, 2401)
    EV, VEC = zmap(B)
    dN0 = np.array([domN(VEC[0, k]) for k in range(M.size)])
    n0 = [k for k in range(M.size) if dN0[k] == 0]
    n1 = [k for k in range(M.size) if dN0[k] == 1]
    cands = []
    for a in n0:
        for b in n1:
            gap = EV[:, b] - EV[:, a]
            for i in np.where(np.diff(np.sign(gap)) != 0)[0]:
                f = gap[i] / (gap[i] - gap[i + 1])
                Bc = B[i] + f * (B[i + 1] - B[i])
                if (par(VEC[i, a]) != par(VEC[i, b]) and
                        domM(VEC[i, a]) == domM(VEC[i, b]) and
                        15050 <= Bc <= 15300):
                    cands.append(dict(Bc=Bc, a=a, b=b, idx=i))
    if not cands:
        record("G2 cross-consistency (reproduce 3.332 kHz drive via d_z)", False,
               f"no canonical crossing found; -muE*d_z==builder op dev={dev_op:.1e}")
        return
    c = min(cands, key=lambda r: abs(r['Bc'] - 15170.0))

    def evoff(Bc, refv, a, b, off):
        _, vo = M.eigensystem(0.0, Bc - off)
        ov = np.abs(vo @ refv.T)
        return vo[int(np.argmax(ov[:, a]))], vo[int(np.argmax(ov[:, b]))]

    va, vb = evoff(c['Bc'], VEC[c['idx']], c['a'], c['b'], 15.0)
    va2, vb2 = evoff(c['Bc'], VEC[c['idx']], c['a'], c['b'], -15.0)
    om_m = abs(va @ (VE_mine * E_DRIVE) @ vb) * 1e3
    om_p = abs(va2 @ (VE_mine * E_DRIVE) @ vb2) * 1e3
    om = 0.5 * (om_m + om_p)
    ratio = om / ORACLE_KHZ
    ok = (dev_op < 1e-9) and (0.95 <= ratio <= 1.05)
    record("G2 cross-consistency (reproduce 3.332 kHz drive via d_z)", ok,
           f"Omega/2pi={om:.4f} kHz at Bc={c['Bc']:.1f} G, E={E_DRIVE} V/cm; "
           f"ratio to oracle {ORACLE_KHZ} kHz = {ratio:.4f}x; "
           f"-muE*d_z == builder op dev={dev_op:.1e}")


# ===========================================================================
# G3 -- the one-line bug fix: TDM_p_odd_aBJ now returns a finite number.
#   (The fix -- adding `return TDM_p` -- ships in a SEPARATE commit; this gate
#   guards against regression.)
# ===========================================================================
def gate3_tdm_return():
    vals = []
    for p in (-1, 0, 1):
        for M0, M1 in ((0.5, 0.5), (0.5, 1.5), (0.5, -0.5)):
            v = me.TDM_p_odd_aBJ(p, 0, 0.5, 0.5, 0.5, 3.0, 3.5, M0,
                                 1, 0.5, 0.5, 1.5, 3.0, 3.5, M1,
                                 S=0.5, I=2.5, iH=0.5)
            vals.append(v)
    not_none = all(v is not None for v in vals)
    finite = not_none and all(np.isfinite(float(v)) for v in vals)
    nonzero = finite and any(abs(float(v)) > 0 for v in vals)
    record("G3 TDM_p_odd_aBJ returns a finite number (missing-return bug fixed)",
           not_none and finite and nonzero,
           f"returns (not None)={not_none}, finite={finite}, "
           f"nonzero-present={nonzero}; sample p=0 value="
           f"{float(vals[3]):+.6f}")


# ===========================================================================
# G4 -- Hermiticity / reality / parity structure.
#   d_z Hermitian; d_x Hermitian; d_x (and d_z) connect ONLY opposite parity
#   (zero matrix elements within same-parity blocks of the zero-field parity
#   eigenbasis).
# ===========================================================================
def gate4_hermiticity_parity():
    M = build_model([0, 1])
    d = build_d(M)
    dz, dx = d['z'], d['x']
    asym_z = float(np.max(np.abs(dz - dz.T)))
    asym_x = float(np.max(np.abs(dx - dx.T)))
    real = bool(np.isrealobj(dz) and np.isrealobj(dx))

    P = np.array(M.Parity_mat)
    pev, pvec = np.linalg.eigh(P)
    U = pvec
    dz_p = U.T @ dz @ U
    dx_p = U.T @ dx @ U
    sgn = np.sign(np.round(pev))
    samep = np.outer(sgn, sgn) > 0
    oppp = ~samep
    samepar_z = float(np.max(np.abs(dz_p[samep])))
    samepar_x = float(np.max(np.abs(dx_p[samep])))
    opp_x = float(np.max(np.abs(dx_p[oppp])))   # must be nonzero (operator nontrivial)

    ok = (asym_z < 1e-12 and asym_x < 1e-12 and real and
          samepar_z < 1e-12 and samepar_x < 1e-12 and opp_x > 1e-6)
    record("G4 Hermiticity + reality + parity structure (d_x, d_z connect "
           "opposite parity only)", ok,
           f"asym d_z={asym_z:.1e}, d_x={asym_x:.1e}, real={real}; "
           f"same-parity |d_z|={samepar_z:.1e}, |d_x|={samepar_x:.1e}; "
           f"opposite-parity |d_x|max={opp_x:.4f}")


if __name__ == "__main__":
    print("=" * 78)
    print("29SiO+ StarkX (p=+-1) + X-X TDM gates")
    print("=" * 78)
    gate1_sum_rule()
    gate2_cross_consistency()
    gate3_tdm_return()
    gate4_hermiticity_parity()

    n_pass = sum(1 for _, p, _ in _results if p)
    n_fail = sum(1 for _, p, _ in _results if not p)
    print("\n" + "=" * 78)
    print(f"{'OK' if n_fail == 0 else 'FAIL'}: {n_pass} passed, {n_fail} failed")
    sys.exit(0 if n_fail == 0 else 1)
