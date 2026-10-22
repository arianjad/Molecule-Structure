"""NSD-PV anapole operator C = (n_hat x S).I/I gates for 29SiO+ X2Sigma+ (Task 4).

Companion to test_sio_crossing.py / test_sio_conventions.py / test_sio_stark_tdm.py.
Validates the case-b(beta-J) NSD-PV operator:
    - matrix element  NSDPV_bBJ  (matrix_elements.py), and
    - builder         build_PTV_NSDPV  (hamiltonian_builders.py, complex dtype), exposed on
      the SiO+ object as MoleculeLevels.NSDPV_operator() -> self.H_NSDPV
against the reviewed, signed-off derivation note docs/sio-nsdpv-derivation.md (commit 7a8a50e).
The note is the spec; values here are NOT tuned -- they are the note's verified gate numbers.

    conda run -n Structure python "Source Code/test_nsdpv_operator.py"

Exit 0 iff every gate passes; exit 1 on any [FAIL].

Gates
-----
G1  Hermiticity (C = C^dagger, machine precision) + P-oddness: C vanishes between
    SAME-parity zero-field eigenstates (via Parity_mat), and the diagonal eigen-expectation
    v.C.v = 0 at 5 field points across 0-16000 G for all real eigenvectors (P-odd diagonal
    null, note 4a).
G2  Decoupled-channel values (note 4b/4c) on directly-constructed decoupled kets, tol 1e-9:
    flip-flop +i/sqrt3; stretched -i/sqrt6; zero-field coupled F0<->F0 = -1.0i;
    F1(J1/2)<->F1(J1/2) = +i/3 (all M); F1(J1/2)<->F1(J3/2) = -i*sqrt2/3.
G3  Karthein Eq. (B1) sign pin (note 4d flag 2) -- PRIMARY GATE. For (bra = even-N |+>
    canonical state, ket = odd-N |-> state) the field-dressed element at B_c is +i*(positive
    real)*(scale), i.e. <+|C|->/i > 0. Uses the +-15 G off-crossing eigenvector method
    (test_sio_crossing.py). The dressed |<C>| must equal 1/sqrt3 = 0.5774 within 3%.
    RECORDED FINDING (2026-06-11, note 4b addendum): dressed value is 0.59015, a stable
    +2.22% above 1/sqrt3 -- stretched-channel admixture in partner B (0.105% weight,
    c-driven, hyperfine-scale intra-manifold gap; offset-independent 0.5-50 G). The note's
    original ~0.5% estimate assumed Zeeman-scale suppression of the channel mixing and is
    superseded; tolerance recalibrated 2% -> 3% (decision A, Arian 2026-06-11, who re-ran
    the gate independently). W/2pi = kprime*W_A*|<C>| at kprime=0.05, W_A/2pi=16 Hz
    = 0.4721 Hz.
G4  Transpose structure: pair-subspace H(W) = H0(B_c) + kprime*W_A*C restricted to the
    crossing pair satisfies H(-W) = H(W)^T in the real basis (manuscript app:transpose), and
    the 2x2 kprime*W_A*C block has the Eq. (B1) form [[.,iW],[-iW,.]] with W real.
G5  PV-active crossing map (deliverable): <C> (field-dressed) at all 7 crossings of the
    Task-3 table -- the five with M_F(A)!=M_F(B) ZERO to machine precision (DM_F=0 kill, note
    2.2); the canonical (0,0) -> +i/sqrt3-class; the (+1,+1) -> stretched -i/sqrt6-class with
    its sign RECORDED. Then finds the SECOND (0,0) crossing (partner
    |N=1,m_N=+1,m_S=-1/2,m_I=-1/2>, also M_F=0) and records its B_c, dressed <C>, sign. Prints
    the full PV-active table (B_c, channel, <C>, sign, W/2pi at kprime=0.05).
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))

from Energy_Levels import MoleculeLevels
from molecule_parameters import get_molecule_params
import matrix_elements as me

# ---------------------------------------------------------------------------
# Constants (note 4b/4c; Karthein assumed kprime=0.05, W_A/2pi=16 Hz)
# ---------------------------------------------------------------------------
KPRIME = 0.05            # assumed mixing coefficient (Karthein Fig. 3)
WA_HZ = 16.0             # W_A/2pi in Hz (Karthein assumed electronic m.e.)
SCAN_HI = 16000.0        # G
WINDOW_LO, WINDOW_HI = 15050.0, 15300.0   # G, Karthein window (Task-3 crossing table)
INV_SQRT3 = 1.0 / np.sqrt(3.0)            # 0.57735, flip-flop |<C>|
INV_SQRT6 = 1.0 / np.sqrt(6.0)            # 0.40825, stretched |<C>|

_results = []  # (name, passed, detail)


def record(name, passed, detail):
    _results.append((name, bool(passed), detail))
    tag = "[PASS]" if passed else "[FAIL]"
    print(f"{tag} {name}: {detail}")


# ---------------------------------------------------------------------------
# Shared model + crossing machinery (reused from test_sio_crossing.py)
# ---------------------------------------------------------------------------
def build_model():
    p = dict(get_molecule_params('SiO+', 'X', '0', 'boson'))
    return MoleculeLevels.initialize_state(
        'SiO+', 'X', 0, N_list=[0, 1, 2], fermion_or_boson='boson',
        I_nuclei=[0, 1 / 2], P_values=[1 / 2], params=p)


def zeeman_map(m, B):
    ev, vec = m.ZeemanMap(B, output=True, write_attribute=False, order=True)
    return np.array(ev), np.array(vec)


def dom_N(m, vec):
    return m.q_numbers['N'][int(np.argmax(vec ** 2))]


def dom_M(m, vec):
    return m.q_numbers['M'][int(np.argmax(vec ** 2))]


def parity_of(m, vec):
    return int(np.round(vec @ m.Parity_mat @ vec))


def evec_offcrossing(m, Bc, ref_vecs, a, b, offset=15.0):
    """Clean partner eigenvectors evaluated OFF the crossing (Bc - offset), re-identified by
    maximum overlap with the tracked-map evecs at the crossing index. Verbatim the method in
    test_sio_crossing.py (the AT-degeneracy evecs are an ill-conditioned arbitrary rotation)."""
    _, vec_off = m.eigensystem(0.0, Bc - offset)
    ov = np.abs(vec_off @ ref_vecs.T)
    ka = int(np.argmax(ov[:, a]))
    kb = int(np.argmax(ov[:, b]))
    return vec_off[ka], vec_off[kb]


def decoupled_components(m, vec, thresh=0.04):
    """Dominant (N, m_N, m_S, m_I, weight) decoupled components."""
    dq = m.alt_q_numbers['decoupled']
    dv = m.convert_evecs('decoupled', evecs=np.array([vec]), verbose=False)[0]
    comps = []
    for idx in np.argsort(-dv ** 2):
        w = dv[idx] ** 2
        if w < thresh:
            break
        comps.append((int(dq['N'][idx]), dq['M_N'][idx], dq['M_S'][idx],
                      dq['M_I'][idx], w))
    return comps


def find_opp_crossings(m, B, EV, VEC, N0_idx, N1_idx):
    out = []
    for a in N0_idx:
        for b in N1_idx:
            gap = EV[:, b] - EV[:, a]
            sc = np.where(np.diff(np.sign(gap)) != 0)[0]
            for i in sc:
                f = gap[i] / (gap[i] - gap[i + 1])
                Bc = B[i] + f * (B[i + 1] - B[i])
                out.append(dict(Bc=Bc, a=a, b=b, idx=i,
                                pa=parity_of(m, VEC[i, a]), pb=parity_of(m, VEC[i, b]),
                                Ma=dom_M(m, VEC[i, a]), Mb=dom_M(m, VEC[i, b])))
    out.sort(key=lambda r: r['Bc'])
    return [c for c in out if c['pa'] != c['pb']]


# ===========================================================================
# Module-level shared setup
# ===========================================================================
print("=" * 78)
print("29SiO+ NSD-PV anapole operator C = (n_hat x S).I/I -- gate suite")
print("=" * 78)

M = build_model()
SIZE = M.size
C = M.NSDPV_operator()                # complex (i/I)*NSDPV matrix, per unit kprime*W_A

B_SCAN = np.linspace(1e-6, SCAN_HI, 2401)
EV, VEC = zeeman_map(M, B_SCAN)
DOMN0 = np.array([dom_N(M, VEC[0, k]) for k in range(SIZE)])
N0_IDX = [k for k in range(SIZE) if DOMN0[k] == 0]
N1_IDX = [k for k in range(SIZE) if DOMN0[k] == 1]
OPP = find_opp_crossings(M, B_SCAN, EV, VEC, N0_IDX, N1_IDX)
WIN = [c for c in OPP if WINDOW_LO <= c['Bc'] <= WINDOW_HI]


def oriented_dressed(m, c, offset=15.0):
    """Field-dressed <C> for crossing c, oriented bra = even-N |+>, ket = odd-N |->.
    Returns (value, va_even, vb_odd)."""
    va, vb = evec_offcrossing(m, c['Bc'], VEC[c['idx']], c['a'], c['b'], offset=offset)
    if dom_N(m, va) != 0:           # ensure bra is the even-N (positive-parity) partner
        va, vb = vb, va
    return va @ C @ vb, va, vb


# ===========================================================================
# G1 -- Hermiticity + P-oddness
# ===========================================================================
def gate1_hermiticity_poddness():
    herm = np.max(np.abs(C - C.conj().T))
    imag_only = np.max(np.abs(C.real))
    herm_ok = herm < 1e-12 and imag_only < 1e-12

    # P-oddness 1: C vanishes between SAME-parity zero-field eigenstates.
    evecs0 = M.evecs0
    par0 = np.array([parity_of(M, v) for v in evecs0])
    same_par_max = 0.0
    for i in range(SIZE):
        for j in range(SIZE):
            if par0[i] == par0[j]:
                same_par_max = max(same_par_max, abs(evecs0[i] @ C @ evecs0[j]))
    same_par_ok = same_par_max < 1e-9

    # P-oddness 2: diagonal eigen-expectation v.C.v = 0 at 5 field points (real evecs).
    diag_max = 0.0
    for Bpt in np.linspace(0.0, SCAN_HI, 5):
        _, vecs = M.eigensystem(0.0, max(Bpt, 1e-6))
        for v in vecs:
            diag_max = max(diag_max, abs((v @ C @ v)))
    diag_ok = diag_max < 1e-9

    ok = herm_ok and same_par_ok and diag_ok
    record("G1 Hermiticity + P-oddness", ok,
           f"max|C-Cdag|={herm:.2e}, max|Re C|={imag_only:.2e} (imaginary-Hermitian); "
           f"same-parity zero-field max|<C>|={same_par_max:.2e} (expect 0); "
           f"diagonal v.C.v max over 5 fields={diag_max:.2e} (P-odd null, expect 0)")


# ===========================================================================
# G2 -- decoupled-channel values on directly-constructed kets
# ===========================================================================
def gate2_decoupled_values():
    I = 0.5
    # Build C in bBJ directly from NSDPV_bBJ, transform to decoupled basis with the model's
    # own audited b_decoupled change-of-basis (decoupled <- bBJ).
    q = M.q_numbers
    qstr = list(q)
    n = len(q[qstr[0]])
    Cbbj = np.zeros((n, n), dtype=complex)
    for i in range(n):
        for j in range(n):
            so = {x + '0': q[x][i] for x in qstr}
            si = {x + '1': q[x][j] for x in qstr}
            Cbbj[i, j] = (1j / I) * me.NSDPV_bBJ(**so, **si, I=I)
    dq = M.alt_q_numbers['decoupled']
    dstr = list(dq)
    nd = len(dq[dstr[0]])
    U = M.library.basis_changers['b_decoupled'](q, dq)   # (nd, n): rows decoupled, cols bBJ
    Cdec = U @ Cbbj @ U.conj().T

    def didx(N, mN, mS, mI):
        for k in range(nd):
            if (dq['N'][k] == N and dq['M_N'][k] == mN and
                    dq['M_S'][k] == mS and dq['M_I'][k] == mI):
                return k
        return None

    tol = 1e-9
    checks = []

    # flip-flop +i/sqrt3
    a = didx(0, 0, 0.5, -0.5); b = didx(1, 0, -0.5, 0.5)
    v = Cdec[a, b]
    checks.append(("flip-flop <0,0,+1/2,-1/2|C|1,0,-1/2,+1/2>", v, 1j * INV_SQRT3))
    # stretched -i/sqrt6
    a = didx(0, 0, 0.5, -0.5); b = didx(1, 1, -0.5, -0.5)
    v = Cdec[a, b]
    checks.append(("stretched <0,0,+1/2,-1/2|C|1,+1,-1/2,-1/2>", v, -1j * INV_SQRT6))

    # zero-field coupled gates -- evaluate NSDPV_bBJ directly in the (N,J,F,M) coupled basis
    def Ccoup(bra, ket):   # bra,ket = (N,J,F,M)
        return (1j / I) * me.NSDPV_bBJ(0, *bra, 0, *ket, I=I)
    checks.append(("F0<->F0 (N0J1/2F0)<->(N1J1/2F0)",
                   Ccoup((0, 0.5, 0, 0), (1, 0.5, 0, 0)), -1.0j))
    for Mm in (-1, 0, 1):
        checks.append((f"F1(J1/2)<->F1(J1/2) M={Mm:+d}",
                       Ccoup((0, 0.5, 1, Mm), (1, 0.5, 1, Mm)), 1j / 3.0))
    checks.append(("F1(J1/2)<->F1(J3/2) (M=0)",
                   Ccoup((0, 0.5, 1, 0), (1, 1.5, 1, 0)), -1j * np.sqrt(2) / 3.0))

    all_ok = True
    lines = []
    for label, got, exp in checks:
        ok = abs(got - exp) < tol
        all_ok = all_ok and ok
        lines.append(f"    {'OK ' if ok else 'BAD'} {label}: <C>={got:+.6f} (exp {exp:+.6f})")
    print("\n".join(lines))
    record("G2 decoupled-channel values (tol 1e-9)", all_ok,
           f"{sum(1 for _,g,e in checks if abs(g-e)<tol)}/{len(checks)} channel "
           f"matrix elements match note 4b/4c exactly")


# ===========================================================================
# G3 -- Karthein Eq. (B1) sign pin (PRIMARY GATE)
# ===========================================================================
CANON = None


def gate3_b1_signpin():
    global CANON
    canon = [c for c in WIN if c['Ma'] == 0 and c['Mb'] == 0]
    canon = min(canon, key=lambda c: abs(c['Bc'] - 15167.2)) if canon else None
    if canon is None:
        record("G3 Karthein Eq.(B1) sign pin (PRIMARY)", False,
               "no canonical (0,0) crossing found in Karthein window")
        return
    CANON = canon
    val_p, va, vb = oriented_dressed(M, canon, offset=15.0)
    val_m, _, _ = oriented_dressed(M, canon, offset=-15.0)
    mag = 0.5 * (abs(val_p) + abs(val_m))
    coeff = (val_p / 1j).real            # <+|C|->/i ; B1 demands > 0
    sign_ok = coeff > 0
    rel = abs(mag - INV_SQRT3) / INV_SQRT3
    mag_ok = rel < 0.03   # recalibrated from 2% (decision A, 2026-06-11); see docstring
    W_hz = KPRIME * WA_HZ * mag
    ok = sign_ok and mag_ok
    record("G3 Karthein Eq.(B1) sign pin (PRIMARY)", ok,
           f"Bc={canon['Bc']:.1f} G; <+|C|->/i={coeff:+.5f} (B1: >0 -> {sign_ok}); "
           f"dressed |<C>|={mag:.5f} vs 1/sqrt3={INV_SQRT3:.5f} ({100*rel:.2f}%, tol 3%); "
           f"dressed value +{100*rel:.2f}% above decoupled 1/sqrt3 -- stretched-channel "
           f"admixture (0.105% weight, c-driven, hyperfine-scale gap); note 4b estimate of "
           f"~0.5% superseded; "
           f"+/-offset agree to {abs(abs(val_p)-abs(val_m))/mag:.2%}; "
           f"W/2pi = kprime*W_A*|<C>| = {W_hz:.4f} Hz (kprime={KPRIME}, W_A/2pi={WA_HZ} Hz)")


# ===========================================================================
# G4 -- transpose structure of the pair subspace
# ===========================================================================
def gate4_transpose():
    if CANON is None:
        record("G4 transpose structure (B1 2x2)", False, "no canonical crossing (G3 failed)")
        return
    val, va, vb = oriented_dressed(M, CANON, offset=15.0)
    basis = [va, vb]                     # bra = even-N |+>, ket = odd-N |->
    H0 = M.H_function(0.0, CANON['Bc'])

    def proj(Op):
        return np.array([[b1.conj() @ Op @ b2 for b2 in basis] for b1 in basis])

    Cpair = proj(C)                      # the bare C block (per unit kprime*W_A)
    H0p = proj(H0)
    HW = H0p + KPRIME * WA_HZ * Cpair    # H(W)
    HmW = H0p - KPRIME * WA_HZ * Cpair   # H(-W) (flip the PV term sign)

    transpose_ok = np.allclose(HmW, HW.T, atol=1e-9)

    # B1 form [[.,iW],[-iW,.]] with W real: the kprime*W_A*C block.
    Cblk = KPRIME * WA_HZ * Cpair
    W_up = Cblk[0, 1] / 1j               # should be +W (real, positive)
    W_lo = Cblk[1, 0] / 1j               # should be -W
    W = W_up.real
    b1_ok = (abs(W_up.imag) < 1e-12 and abs(W_lo.imag) < 1e-12 and
             abs(Cblk[0, 0]) < 1e-12 and abs(Cblk[1, 1]) < 1e-12 and
             abs(W_lo.real + W) < 1e-12 and W > 0)

    ok = transpose_ok and b1_ok
    record("G4 transpose structure (B1 2x2)", ok,
           f"H(-W) == H(W)^T: {transpose_ok}; kprime*W_A*C block = "
           f"[[{Cblk[0,0]:.3g}, {Cblk[0,1]:.4f}], [{Cblk[1,0]:.4f}, {Cblk[1,1]:.3g}]]; "
           f"upper/i=+W={W_up.real:+.5f}, lower/i={W_lo.real:+.5f} (W real, >0: {b1_ok})")


# ===========================================================================
# G5 -- PV-active crossing map (deliverable)
# ===========================================================================
def gate5_crossing_map():
    # All 7 window crossings (Task-3 table): five DM_F!=0 must be ZERO, two DM_F=0 nonzero.
    print("\n  --- All opposite-parity N0xN1 crossings in Karthein window (Task-3 table) ---")
    print(f"  {'Bc[G]':>9} {'Ma':>4} {'Mb':>4} {'DM_F=0':>7} {'|<C>|dressed':>13} "
          f"{'<C>/i (sign)':>13}")
    five_zero_max = 0.0
    canon_row = None
    stretched_row = None
    for c in WIN:
        dm0 = (c['Ma'] == c['Mb'])
        if not dm0:
            va, vb = evec_offcrossing(M, c['Bc'], VEC[c['idx']], c['a'], c['b'])
            val = va @ C @ vb
            five_zero_max = max(five_zero_max, abs(val))
            print(f"  {c['Bc']:9.1f} {c['Ma']:+4.0f} {c['Mb']:+4.0f} {'no':>7} "
                  f"{abs(val):13.2e} {'0 (DM_F kill)':>13}")
        else:
            val, va, vb = oriented_dressed(M, c)
            coeff = (val / 1j).real
            chan = "(0,0) flip-flop" if (c['Ma'] == 0) else "(+1,+1) stretched"
            print(f"  {c['Bc']:9.1f} {c['Ma']:+4.0f} {c['Mb']:+4.0f} {'YES':>7} "
                  f"{abs(val):13.5f} {coeff:+13.5f}")
            if c['Ma'] == 0:
                canon_row = (c['Bc'], "(0,0) flip-flop", abs(val), coeff)
            else:
                stretched_row = (c['Bc'], "(+1,+1) stretched", abs(val), coeff)

    five_zero_ok = five_zero_max < 1e-12

    # Second (0,0) crossing: partner |N=1, m_N=+1, m_S=-1/2, m_I=-1/2> (stretched-channel
    # partner of A; also M_F=0), hyperfine-shifted from 15167.2 G. Scan wider window.
    Bw = np.linspace(14600.0, 15800.0, 4001)
    EVw, VECw = zeeman_map(M, Bw)
    domn0w = np.array([dom_N(M, VECw[0, k]) for k in range(SIZE)])
    n0w = [k for k in range(SIZE) if domn0w[k] == 0]
    n1w = [k for k in range(SIZE) if domn0w[k] == 1]
    oppw = find_opp_crossings(M, Bw, EVw, VECw, n0w, n1w)
    zero_zero_w = [c for c in oppw if c['Ma'] == 0 and c['Mb'] == 0]
    second = None
    for c in zero_zero_w:
        if abs(c['Bc'] - 15167.2) < 5.0:    # skip the canonical one
            continue
        va, vb = evec_offcrossing(M, c['Bc'], VECw[c['idx']], c['a'], c['b'])
        if dom_N(M, va) != 0:
            va, vb = vb, va
        comps = decoupled_components(M, vb)
        # confirm the stretched partner character
        is_stretched = any(abs(n[1] - 1.0) < 1e-9 and abs(n[2] + 0.5) < 1e-9 and
                           abs(n[3] + 0.5) < 1e-9 for n in comps)
        val = va @ C @ vb
        second = (c['Bc'], "(0,0)' stretched-partner", abs(val), (val / 1j).real, is_stretched)

    print("\n  --- PV-ACTIVE crossing table (DM_F=0 only; the deliverable) ---")
    print(f"  {'Bc[G]':>9} {'channel':>26} {'|<C>|':>9} {'<C>/i sign':>11} "
          f"{'W/2pi[Hz]':>10}")
    pv_rows = []
    if canon_row:
        pv_rows.append(canon_row)
    if second:
        pv_rows.append(second[:4])
    if stretched_row:
        pv_rows.append(stretched_row)
    for r in pv_rows:
        Bc, chan, mag, coeff = r[0], r[1], r[2], r[3]
        sign = "+" if coeff > 0 else "-"
        W = KPRIME * WA_HZ * mag
        print(f"  {Bc:9.1f} {chan:>26} {mag:9.5f} {sign:>11} {W:10.4f}")

    second_ok = second is not None and second[4]   # found + confirmed stretched character
    map_ok = (five_zero_ok and canon_row is not None and stretched_row is not None
              and second_ok)
    record("G5 PV-active crossing map (deliverable)", map_ok,
           f"5 DM_F!=0 crossings: max|<C>|={five_zero_max:.2e} (expect 0, DM_F kill); "
           f"canonical (0,0) at {canon_row[0]:.1f}G |<C>|={canon_row[2]:.5f} "
           f"({'+' if canon_row[3]>0 else '-'}); "
           f"(+1,+1) at {stretched_row[0]:.1f}G |<C>|={stretched_row[2]:.5f} "
           f"({'+' if stretched_row[3]>0 else '-'}); "
           f"2nd (0,0) at {second[0]:.1f}G |<C>|={second[2]:.5f} "
           f"({'+' if second[3]>0 else '-'}, stretched-partner confirmed={second[4]})")


if __name__ == "__main__":
    gate1_hermiticity_poddness()
    gate2_decoupled_values()
    gate3_b1_signpin()
    gate4_transpose()
    gate5_crossing_map()

    n_pass = sum(1 for _, p, _ in _results if p)
    n_fail = sum(1 for _, p, _ in _results if not p)
    print("\n" + "=" * 78)
    print(f"{'OK' if n_fail == 0 else 'FAIL'}: {n_pass} passed, {n_fail} failed")
    sys.exit(0 if n_fail == 0 else 1)
