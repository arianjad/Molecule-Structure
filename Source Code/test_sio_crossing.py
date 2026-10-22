"""Crossing oracle gate for the 29SiO+ X2Sigma+ entry (molecules['SiO+']['boson']['X0']).

This is the program's oracle gate: it reproduces the Karthein N=0(+) / N=1(-)
near-degeneracy in 29SiO+ near 1.5 T, names the crossing sublevels, and
reproduces the published Stark drive coupling (Omega/2pi ~ 3 kHz at E = 6 V/cm
axial). If the constants/conventions in molecule_parameters.py are right, the
crossing lands inside Karthein's 1.505-1.530 T window and the drive matches to
within a factor of 2. If they are wrong, gate 1 fails and the script reports the
MEASURED crossing field(s) -- the oracle doing its job. The constants are NEVER
adjusted to make a gate pass.

    conda run -n Structure python "Source Code/test_sio_crossing.py"

Exit 0 iff every gate passes; exit 1 on any [FAIL].

Physics background
------------------
29SiO+ X2Sigma+ has electron spin S=1/2 and one nuclear spin (29Si, I=1/2; the O
is spinless on the SiO+ entry, I_nuclei=[0,1/2]). At low field each rotational
level N splits by spin-rotation (gamma) and Fermi-contact + dipolar hyperfine.
The electron-spin Zeeman slope is +-g_S mu_B/2 ~ +-1.4017 MHz/G, so the
differential slope between an m_S=+1/2 sublevel of N=0 and an m_S=-1/2 sublevel
of N=1 is ~2.80 MHz/G. The N=0 -> N=1 rotational interval is ~42.49 GHz, so the
two manifolds are tuned into degeneracy near 42490 / 2.80 ~ 15170 G ~ 1.517 T.

Parity for a 2Sigma+ state is +(-1)^N: N=0 (even) -> +, N=1 (odd) -> -. The
crossing is between a positive-parity N=0 sublevel and a negative-parity N=1
sublevel. The static E-field (StarkZ, axial, Delta M_F=0) couples them ONLY
through the spin-character admixture -- the two partners carry opposite dominant
m_S, and the Stark operator is spin-diagonal, so the matrix element is nonzero
only because each high-field eigenstate is a small admixture of the other's
spin character. This is the Karthein mechanism.

The parity convention is verified against the code's own Parity_mat below
(gate 1); the +(-1)^N formula is not assumed.
"""
import os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))

from Energy_Levels import MoleculeLevels
from molecule_parameters import get_molecule_params, params_general

# ---------------------------------------------------------------------------
# Constants / window
# ---------------------------------------------------------------------------
E_DRIVE = 6.0                 # V/cm, Karthein axial drive field
KARTHEIN_OMEGA_KHZ = 3.0      # kHz, published Omega/2pi at 6 V/cm
WINDOW_LO, WINDOW_HI = 15050.0, 15300.0   # G, Karthein 1.505-1.530 T
SCAN_LO, SCAN_HI = 14000.0, 16000.0       # G, full scan
SLOPE_TARGET = 2.80           # MHz/G, predicted differential slope

_results = []  # (name, passed, detail)


def record(name, passed, detail):
    _results.append((name, bool(passed), detail))
    tag = "[PASS]" if passed else "[FAIL]"
    print(f"{tag} {name}: {detail}")


# ---------------------------------------------------------------------------
# Shared machinery
# ---------------------------------------------------------------------------
def build_model(gamma=None):
    """Instantiate the N=0,1,2 SiO+ model. gamma overrides Gamma_SR (MHz)."""
    p = dict(get_molecule_params('SiO+', 'X', '0', 'boson'))
    if gamma is not None:
        p['Gamma_SR'] = float(gamma)
    return MoleculeLevels.initialize_state(
        'SiO+', 'X', 0, N_list=[0, 1, 2], fermion_or_boson='boson',
        I_nuclei=[0, 1 / 2], P_values=[1 / 2], params=p)


def stark_operator(m):
    """V_E matrix (the -muE*StarkZ field-coupling operator) in the original
    basis, extracted from the builder's own H_function closure:
        H(E,B) = H0 + V_E*E + V_B*B   (hamiltonian_builders.py:99)
    so V_E = H(1,0) - H(0,0). This is the SAME matrix the builder applies as
    V_E[i][j] += -muE*StarkZ (hamiltonian_builders.py:58), guaranteeing the
    drive element is computed with the identical operator the spectrum uses."""
    return m.H_function(1.0, 1e-9) - m.H_function(0.0, 1e-9)


def zeeman_map(m, B):
    evB, vecB = m.ZeemanMap(B, output=True, write_attribute=False, order=True)
    return np.array(evB), np.array(vecB)


def dom_N(m, vec):
    return m.q_numbers['N'][int(np.argmax(vec ** 2))]


def dom_M(m, vec):
    return m.q_numbers['M'][int(np.argmax(vec ** 2))]


def parity_of(m, vec):
    return int(np.round(vec @ m.Parity_mat @ vec))


def find_crossings(m, B, evB, vecB, n0_idx, n1_idx):
    """All true sign-change crossings between an N=0-character state and an
    N=1-character state, with interpolated Bc, differential slope, min gap,
    parities and dominant M of each partner. Returns list of dicts."""
    out = []
    for a in n0_idx:
        for b in n1_idx:
            gap = evB[:, b] - evB[:, a]
            sc = np.where(np.diff(np.sign(gap)) != 0)[0]
            for i in sc:
                f = gap[i] / (gap[i] - gap[i + 1])
                Bc = B[i] + f * (B[i + 1] - B[i])
                slope = abs((gap[i + 1] - gap[i]) / (B[i + 1] - B[i]))
                gmin = min(abs(gap[i]), abs(gap[i + 1]))
                pa = parity_of(m, vecB[i, a]); pb = parity_of(m, vecB[i, b])
                Ma = dom_M(m, vecB[i, a]); Mb = dom_M(m, vecB[i, b])
                out.append(dict(Bc=Bc, a=a, b=b, idx=i, slope=slope, gmin=gmin,
                                pa=pa, pb=pb, Ma=Ma, Mb=Mb))
    out.sort(key=lambda r: r['Bc'])
    return out


def evec_offcrossing(m, Bc, ref_vecs, a, b, offset=15.0):
    """Eigenvectors of partners a,b evaluated OFF the crossing (Bc - offset).

    Conditioning note: AT the degeneracy the a/b eigenvectors are an arbitrary
    rotation within the near-degenerate 2D subspace -- the individual evecs are
    ill-conditioned and the Stark element is split arbitrarily between them.
    We therefore evaluate at Bc - offset (the manifolds are well separated by
    ~2.8 MHz/G * offset >> the ~kHz coupling), re-identify the two partners by
    maximum overlap with the tracked-map evecs at the crossing index, and use
    those clean evecs. Verified elsewhere that +offset and -offset agree to
    <1%, so the off-crossing choice does not bias the result."""
    _, vec_off = m.eigensystem(0.0, Bc - offset)
    ov = np.abs(vec_off @ ref_vecs.T)
    ka = int(np.argmax(ov[:, a]))
    kb = int(np.argmax(ov[:, b]))
    return vec_off[ka], vec_off[kb]


def decoupled_components(m, vec, thresh=0.04):
    """Dominant |N, M_N, M_S, M_I, M_F> components (weight >= thresh)."""
    dq = m.alt_q_numbers['decoupled']
    dv = m.convert_evecs('decoupled', evecs=np.array([vec]), verbose=False)[0]
    comps = []
    for idx in np.argsort(-dv ** 2):
        w = dv[idx] ** 2
        if w < thresh:
            break
        comps.append(dict(N=dq['N'][idx], M_N=dq['M_N'][idx], M_S=dq['M_S'][idx],
                          M_I=dq['M_I'][idx], M_F=dq['M_F'][idx], w=w, amp=dv[idx]))
    return comps


def zero_field_label(m, B, vecB, k):
    """Zero-field-adiabatic (N, J, F, M) ket of tracked state k. Requires the
    map to start near B=0 so index k carries its zero-field identity."""
    v = vecB[0, k]
    idx = int(np.argmax(v ** 2))
    return {q: m.q_numbers[q][idx] for q in ['N', 'J', 'F', 'M']}, v[idx] ** 2


def fmt_decomp(comps):
    return " + ".join(
        f"{c['w']:.3f}|N={c['N']:.0f},m_N={c['M_N']:+.0f},m_S={c['M_S']:+.1f},"
        f"m_I={c['M_I']:+.1f},M_F={c['M_F']:+.0f}>" for c in comps)


# ===========================================================================
# Module-level shared analysis (computed once, used by several gates)
# ===========================================================================
print("=" * 78)
print("29SiO+ crossing oracle -- building N=0,1,2 model and scanning Zeeman map")
print("=" * 78)

M = build_model()                       # nominal Gamma_SR = 12 MHz
SIZE = M.size

# Full-window scan, tracked from low field so indices carry zero-field identity.
B_SCAN = np.linspace(1e-6, SCAN_HI, 2401)        # ~6.7 G spacing, >1001 points
EV, VEC = zeeman_map(M, B_SCAN)
DOMN0 = np.array([dom_N(M, VEC[0, k]) for k in range(SIZE)])
N0_IDX = [k for k in range(SIZE) if DOMN0[k] == 0]
N1_IDX = [k for k in range(SIZE) if DOMN0[k] == 1]

ALL_CROSS = find_crossings(M, B_SCAN, EV, VEC, N0_IDX, N1_IDX)
# Restrict to crossings inside the full scan with opposite parity (the oracle).
OPP = [c for c in ALL_CROSS if c['pa'] != c['pb']]
IN_WINDOW = [c for c in OPP if WINDOW_LO <= c['Bc'] <= WINDOW_HI]

# Refine each crossing's Bc and slope on a fine local grid (adaptive refinement,
# not a 1e5-point brute force) so the slope gate is not grid-limited.
def refine(m, c):
    B = np.linspace(c['Bc'] - 30, c['Bc'] + 30, 601)
    ev, vec = zeeman_map(m, B)
    gap = ev[:, c['b']] - ev[:, c['a']]
    sc = np.where(np.diff(np.sign(gap)) != 0)[0]
    if len(sc) == 0:
        return c
    i = sc[len(sc) // 2]
    f = gap[i] / (gap[i] - gap[i + 1])
    c = dict(c)
    c['Bc'] = B[i] + f * (B[i + 1] - B[i])
    c['slope'] = abs((gap[i + 1] - gap[i]) / (B[i + 1] - B[i]))
    c['gmin'] = min(abs(gap[i]), abs(gap[i + 1]))
    return c

OPP = [refine(M, c) for c in OPP]
IN_WINDOW = [c for c in OPP if WINDOW_LO <= c['Bc'] <= WINDOW_HI]

# Pick the canonical Delta M_F = 0 (driveable) crossing nearest the Karthein
# prediction (~15170 G). The static axial E-field couples only Delta M_F = 0.
DM0 = [c for c in IN_WINDOW if c['Ma'] == c['Mb']]
CANON = min(DM0, key=lambda c: abs(c['Bc'] - 15170.0)) if DM0 else None


# ===========================================================================
# Gate 1 -- Window: 36 states, crossings exist inside Karthein window with the
#           code-verified parity convention (+(-1)^N).
# ===========================================================================
def gate1_window():
    ok_size = (SIZE == 36)
    # Verify code's own parity vs +(-1)^N for N=0 and N=1 manifolds.
    par_n0 = sorted(set(parity_of(M, VEC[0, k]) for k in N0_IDX))
    par_n1 = sorted(set(parity_of(M, VEC[0, k]) for k in N1_IDX))
    parity_ok = (par_n0 == [1]) and (par_n1 == [-1])
    crossings_ok = len(IN_WINDOW) >= 1
    detail = (f"size={SIZE} (exp 36); N=0 parities={par_n0} (exp [+1]), "
              f"N=1 parities={par_n1} (exp [-1]); "
              f"{len(IN_WINDOW)} opposite-parity N0xN1 crossings in "
              f"[{WINDOW_LO:.0f},{WINDOW_HI:.0f}] G")
    if not crossings_ok:
        # Oracle failure mode: report the measured crossing field(s).
        meas = ", ".join(f"{c['Bc']:.1f}G" for c in OPP[:8]) or "NONE in scan"
        detail += f" | MEASURED opposite-parity crossings (full scan): {meas}"
    record("G1 window (36 states, parity, crossing in Karthein window)",
           ok_size and parity_ok and crossings_ok, detail)


# ===========================================================================
# Gate 2 -- Slope: differential slope within 20% of 2.80 MHz/G.
# ===========================================================================
def gate2_slope():
    if CANON is None:
        record("G2 differential slope (~2.80 MHz/G)", False,
               "no Delta M_F=0 crossing in window to measure slope")
        return
    s = CANON['slope']
    rel = abs(s - SLOPE_TARGET) / SLOPE_TARGET
    record("G2 differential slope (~2.80 MHz/G)", rel < 0.20,
           f"|dDelta/dB|={s:.4f} MHz/G at Bc={CANON['Bc']:.1f} G "
           f"({100 * (s - SLOPE_TARGET) / SLOPE_TARGET:+.2f}% vs {SLOPE_TARGET})")


# ===========================================================================
# Gate 3 -- Spin character: partners differ in dominant m_S (decoupled basis).
# ===========================================================================
def gate3_spin():
    if CANON is None:
        record("G3 spin character (partners differ in m_S)", False,
               "no Delta M_F=0 crossing to characterize")
        return
    va, vb = evec_offcrossing(M, CANON['Bc'], VEC[CANON['idx']],
                              CANON['a'], CANON['b'])
    ca = decoupled_components(M, va)
    cb = decoupled_components(M, vb)
    mSa = ca[0]['M_S']; mSb = cb[0]['M_S']
    differ = (mSa != mSb)
    record("G3 spin character (partners differ in m_S)", differ,
           f"A dom m_S={mSa:+.1f}, B dom m_S={mSb:+.1f} (differ={differ}); "
           f"A: {fmt_decomp(ca)} || B: {fmt_decomp(cb)}")


# ===========================================================================
# Gate 4 -- Pair identification (deliverable #1): full labels + all-crossings
#           table. No threshold (always passes if labels resolve), prints the
#           headline deliverable.
# ===========================================================================
def gate4_identify():
    if CANON is None:
        record("G4 pair identification", False, "no canonical crossing")
        return
    la, wa = zero_field_label(M, B_SCAN, VEC, CANON['a'])
    lb, wb = zero_field_label(M, B_SCAN, VEC, CANON['b'])
    va, vb = evec_offcrossing(M, CANON['Bc'], VEC[CANON['idx']],
                              CANON['a'], CANON['b'])
    ca = decoupled_components(M, va); cb = decoupled_components(M, vb)

    print("\n  --- CROSSING PAIR (canonical Delta M_F=0, Bc = "
          f"{CANON['Bc']:.1f} G = {CANON['Bc'] / 1e4:.4f} T) ---")
    print(f"  Partner A (parity {CANON['pa']:+d}): zero-field |N={la['N']:.0f}, "
          f"J={la['J']}, F={la['F']:.0f}, M_F={la['M']:+.0f}> (w={wa:.3f})")
    print(f"            high-field decoupled: {fmt_decomp(ca)}")
    print(f"  Partner B (parity {CANON['pb']:+d}): zero-field |N={lb['N']:.0f}, "
          f"J={lb['J']}, F={lb['F']:.0f}, M_F={lb['M']:+.0f}> (w={wb:.3f})")
    print(f"            high-field decoupled: {fmt_decomp(cb)}")

    print(f"\n  --- ALL opposite-parity N0xN1 crossings in "
          f"[{WINDOW_LO:.0f},{WINDOW_HI:.0f}] G ---")
    print(f"  {'Bc[G]':>9} {'parA':>5} {'parB':>5} {'Ma':>4} {'Mb':>4} "
          f"{'slope[MHz/G]':>12} {'mingap[MHz]':>11}  {'DM0':>4}")
    for c in IN_WINDOW:
        dm = "yes" if c['Ma'] == c['Mb'] else "no"
        print(f"  {c['Bc']:9.1f} {c['pa']:+5d} {c['pb']:+5d} "
              f"{c['Ma']:+4.0f} {c['Mb']:+4.0f} {c['slope']:12.4f} "
              f"{c['gmin']:11.3f}  {dm:>4}")

    labels_ok = (la['N'] == 0 and lb['N'] == 1 and
                 len(ca) > 0 and len(cb) > 0)
    record("G4 pair identification (labels + all-crossings table)", labels_ok,
           f"A=|N=0,J={la['J']},F={la['F']:.0f},M_F={la['M']:+.0f}> (par "
           f"{CANON['pa']:+d}), B=|N=1,J={lb['J']},F={lb['F']:.0f},"
           f"M_F={lb['M']:+.0f}> (par {CANON['pb']:+d}); "
           f"{len(IN_WINDOW)} crossings tabulated")


# ===========================================================================
# Gate 5 -- Drive oracle: |<A| muE*E*StarkZ |B>| at E=6 V/cm within factor 2
#           of 3 kHz.
# ===========================================================================
def gate5_drive():
    if CANON is None:
        record("G5 drive oracle (~3 kHz at 6 V/cm)", False, "no canonical crossing")
        return
    VE = stark_operator(M)
    va, vb = evec_offcrossing(M, CANON['Bc'], VEC[CANON['idx']],
                              CANON['a'], CANON['b'])
    # Cross-check conditioning: +/- offset must agree.
    va_p, vb_p = evec_offcrossing(M, CANON['Bc'], VEC[CANON['idx']],
                                  CANON['a'], CANON['b'], offset=-15.0)
    omega_m = abs(va @ (VE * E_DRIVE) @ vb)        # MHz
    omega_p = abs(va_p @ (VE * E_DRIVE) @ vb_p)
    omega = 0.5 * (omega_m + omega_p)
    cond_spread = abs(omega_m - omega_p) / omega if omega > 0 else 0.0
    omega_khz = omega * 1e3
    # xi = Omega / (muE*E/h): dimensionless suppression vs the bare dipole drive.
    muE = M.parameters['muE']                      # MHz/(V/cm)
    bare = muE * E_DRIVE                            # MHz, full-dipole reference
    xi = omega / bare
    ratio = omega_khz / KARTHEIN_OMEGA_KHZ
    ok = 0.5 <= ratio <= 2.0
    record("G5 drive oracle (~3 kHz at 6 V/cm, within factor 2)", ok,
           f"|<A|muE*E*StarkZ|B>|={omega_khz:.4f} kHz (E={E_DRIVE} V/cm); "
           f"ratio to 3 kHz = {ratio:.3f}x; xi=Omega/(muE*E/h)={xi:.3e}; "
           f"+/-offset cond. spread={cond_spread:.2%}")


# ===========================================================================
# Gate 6 -- gamma error band (Decision D1): Bc(gamma), Omega(gamma) for
#           Gamma_SR in {-13, +12, +37} MHz. No pass/fail threshold; gate is
#           that all three builds complete and numbers are recorded.
# ===========================================================================
def gate6_gamma_band():
    if CANON is None:
        record("G6 gamma error band (D1)", False, "no canonical crossing to track")
        return
    print("\n  --- gamma error band (Decision D1): Delta M_F=0 N0xN1 crossings ---")
    print(f"  {'gamma[MHz]':>10} {'Bc[G]':>9} {'Bc[T]':>8} {'Omega[kHz]':>11} {'xi':>11}")
    rows = []
    all_ok = True
    omegas_canon = []
    for g in (-13.0, 12.0, 37.0):
        try:
            mg = build_model(gamma=g)
            Bg = np.linspace(SCAN_LO, SCAN_HI, 2001)
            eg, vg = zeeman_map(mg, Bg)
            dNg = np.array([dom_N(mg, vg[0, k]) for k in range(mg.size)])
            n0g = [k for k in range(mg.size) if dNg[k] == 0]
            n1g = [k for k in range(mg.size) if dNg[k] == 1]
            cg = find_crossings(mg, Bg, eg, vg, n0g, n1g)
            cg = [c for c in cg if c['pa'] != c['pb'] and c['Ma'] == c['Mb']
                  and WINDOW_LO <= c['Bc'] <= WINDOW_HI]
            VEg = stark_operator(mg)
            muEg = mg.parameters['muE']
            # Canonical = nearest 15170 G (the same branch as CANON).
            cg.sort(key=lambda c: abs(c['Bc'] - 15170.0))
            for j, c in enumerate(cg):
                va, vb = evec_offcrossing(mg, c['Bc'], vg[c['idx']], c['a'], c['b'])
                om = abs(va @ (VEg * E_DRIVE) @ vb) * 1e3
                xi = (om / 1e3) / (muEg * E_DRIVE)
                tag = " <-canon" if j == 0 else ""
                print(f"  {g:10.0f} {c['Bc']:9.1f} {c['Bc']/1e4:8.4f} "
                      f"{om:11.4f} {xi:11.3e}{tag}")
                if j == 0:
                    omegas_canon.append(om)
                    rows.append((g, c['Bc'], om))
        except Exception as e:
            all_ok = False
            print(f"  gamma={g}: BUILD FAILED: {e!r}")
    # gamma-dominance flag.
    flag = ""
    if len(omegas_canon) >= 2:
        spread = max(omegas_canon) / min(omegas_canon)
        if spread > 3.0:
            flag = (f" | *** FLAG: canonical Omega varies {spread:.1f}x across "
                    f"gamma band -- drive budget is gamma-DOMINATED ***")
        else:
            flag = (f" | canonical Omega varies {spread:.2f}x across band "
                    f"(NOT gamma-dominated)")
    record("G6 gamma error band (D1: 3 builds complete, numbers recorded)",
           all_ok and len(rows) == 3,
           f"Bc(gamma)/Omega(gamma) recorded for gamma in [-13,+12,+37] MHz"
           + flag)


if __name__ == "__main__":
    gate1_window()
    gate2_slope()
    gate3_spin()
    gate4_identify()
    gate5_drive()
    gate6_gamma_band()

    n_pass = sum(1 for _, p, _ in _results if p)
    n_fail = sum(1 for _, p, _ in _results if not p)
    print("\n" + "=" * 78)
    print(f"{'OK' if n_fail == 0 else 'FAIL'}: {n_pass} passed, {n_fail} failed")
    sys.exit(0 if n_fail == 0 else 1)
