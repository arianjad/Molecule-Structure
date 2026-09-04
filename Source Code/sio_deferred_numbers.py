"""Deferred-number package I for 29SiO+ X2Sigma+ (Task 6): xi(B)/d_eff, Protocol-I
E-field budget, the B_perp x E_z cross term, and the B_perp leakage table.

Companion to test_sio_crossing.py / test_sio_conventions.py / test_sio_stark_tdm.py /
test_nsdpv_operator.py. This script SHIPS numbers for the manuscript and the
reversal-suite budget; it does NOT introduce new physics. Every operator it uses
is the codebase's own:

    d_z  == StarkZ matrix   (the trusted axial dipole; -muE*d_z == H_function(1,0)-H_function(0,0))
    zx   == ZeemanX matrix  (transverse electron-spin Zeeman, p=+-1, registered in
            bBJ_even_X ext_fields), applied with the SAME coefficient the builder
            uses for the axial Zeeman: V_B += g_S*mu_B*ZeemanZ (hamiltonian_builders.py:53),
            so the transverse partner is  H_Bx = g_S*mu_B*B_x*zx.

    conda run -n Structure python "Source Code/sio_deferred_numbers.py"

Exit 0 iff every gate passes; exit 1 on any [FAIL]. Gates are SELF-CONSISTENCY
assertions (arithmetic / parity / reality / coverage), never literature numbers.

The canonical pair (test_sio_crossing.py G4/G5, test_nsdpv_operator.py G3):
    A = |N=0, J=1/2, F=0, M_F=0> (parity +, rises +1.40166 MHz/G)   ~ |0,0,+1/2,-1/2> decoupled
    B = |N=1, J=3/2, F=2, M_F=0> (parity -, falls -1.40166 MHz/G)   ~ |1,0,-1/2,+1/2> decoupled
B_c = 15167.2 G; xi(B_c) = |<A|d_z|B>| = 2.660e-4 (unit-muE dipole fraction);
Omega/2pi = 3.3319 kHz at E_z = 6 V/cm.

RECORDED FINDING (Task-6 deliverable 3, 2026-06-11). The B_perp x E_z bilinear
coefficient d2M/dB_x dE_z is ZERO to machine precision, NOT a small finite number.
This is the strong, quantitative form of the static-fields-cannot-fake-sigma_y
theorem: VBx = g_S*mu_B*zx changes M_F by +-1, the Stark drive VE = -muE*d_z keeps
M_F=0, and the M_F=0 pair (A,B) therefore has NO linear-in-B_x coupling at any order
combined with the axial drive -- both the literal operator element <A|(E_z*VE +
B_x*VBx)|B> and the second-order van Vleck cross term vanish by M_F selection. The
leading genuine B_x systematic on the drive is QUADRATIC in B_x (deliverable 3b ships
that coefficient + the relative drive shift at 10 mG as the reversal-suite budget
number). The pure-B_x pair element at E=0 is exactly zero (parity/M_F exact, 3c).
The drive coupling stays REAL (sigma_x-class) to machine precision for any combination
of static real fields (3a) -- no imaginary sigma_y can be faked.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))

from Energy_Levels import MoleculeLevels
from molecule_parameters import get_molecule_params
import hamiltonian_builders as ham

# ---------------------------------------------------------------------------
# Operating point / constants (manuscript Protocol I; Karthein crossing)
# ---------------------------------------------------------------------------
E_DRIVE = 6.0                 # V/cm, Karthein axial drive field
XI_TARGET = 2.660e-4          # trusted unit-muE dipole fraction at B_c (test_sio_crossing G5)
DEBYE_PER_UNIT_MUE = 4.147    # muE = 4.147 D * 0.503412 MHz/(V/cm); 1 unit-muE = 4.147 D
WINDOW_LO, WINDOW_HI = 15050.0, 15300.0   # G, Karthein window for crossing search
SWEEP_LO, SWEEP_HI, SWEEP_N = 14950.0, 15350.0, 21   # G, deliverable-1 xi(B) sweep
GAMMA_BAND = (-13.0, 12.0, 37.0)          # MHz, Gamma_SR error band (Decision D1)
# Protocol-I manuscript operating point (recorded, used in deliverable 2 header):
TAU_MS, N_CYC, T_TOTAL_S = 0.78, 256, 0.2
TP_LIST_US = (1.0, 5.0, 10.0, 50.0)       # pi_x pulse durations, microseconds
EPI_WATERMARK = 100.0         # V/cm, trap-feasibility flag (comment, not gate)
BX_LIST_MG = (1.0, 5.0, 10.0) # mG, transverse-B probes for the cross term

_results = []  # (name, passed, detail)


def record(name, passed, detail):
    _results.append((name, bool(passed), detail))
    tag = "[PASS]" if passed else "[FAIL]"
    print(f"{tag} {name}: {detail}")


# ---------------------------------------------------------------------------
# Shared machinery (verbatim methods from the companion test files)
# ---------------------------------------------------------------------------
def build_model(gamma=None):
    """N=0,1,2 SiO+ boson X0 model. gamma overrides Gamma_SR (MHz)."""
    p = dict(get_molecule_params('SiO+', 'X', '0', 'boson'))
    if gamma is not None:
        p['Gamma_SR'] = float(gamma)
    return MoleculeLevels.initialize_state(
        'SiO+', 'X', 0, N_list=[0, 1, 2], fermion_or_boson='boson',
        I_nuclei=[0, 1 / 2], P_values=[1 / 2], params=p)


def op_matrix(m, name):
    """Retrieve a registered ext-field operator as a matrix over the model basis,
    via the same hamiltonian_builders.build_operator path used by test_sio_stark_tdm
    build_d. 'StarkZ' -> d_z; 'ZeemanX' -> zx (transverse electron-spin Zeeman)."""
    qn = m.q_numbers
    return np.array(ham.build_operator(qn, qn, m.matrix_elements[name]))


def zeeman_map(m, B):
    ev, vec = m.ZeemanMap(B, output=True, write_attribute=False, order=True)
    return np.array(ev), np.array(vec)


def dom_N(m, vec):
    return m.q_numbers['N'][int(np.argmax(vec ** 2))]


def dom_M(m, vec):
    return m.q_numbers['M'][int(np.argmax(vec ** 2))]


def dom_q(m, vec, q):
    return m.q_numbers[q][int(np.argmax(vec ** 2))]


def parity_of(m, vec):
    return int(np.round(vec @ m.Parity_mat @ vec))


def find_canonical(m):
    """The canonical Delta M_F=0, opposite-parity N=0xN=1 crossing nearest 15167.2 G,
    tracked from low field so the partner indices a,b carry their zero-field identity.
    Returns (crossing_dict, B_scan, EV, VEC)."""
    BS = np.linspace(1e-6, 16000.0, 2401)
    EV, VEC = zeeman_map(m, BS)
    domn0 = np.array([dom_N(m, VEC[0, k]) for k in range(m.size)])
    n0 = [k for k in range(m.size) if domn0[k] == 0]
    n1 = [k for k in range(m.size) if domn0[k] == 1]
    cands = []
    for a in n0:
        for b in n1:
            gap = EV[:, b] - EV[:, a]
            for i in np.where(np.diff(np.sign(gap)) != 0)[0]:
                f = gap[i] / (gap[i] - gap[i + 1])
                Bc = BS[i] + f * (BS[i + 1] - BS[i])
                if (parity_of(m, VEC[i, a]) != parity_of(m, VEC[i, b]) and
                        dom_M(m, VEC[i, a]) == 0 and dom_M(m, VEC[i, b]) == 0 and
                        WINDOW_LO <= Bc <= WINDOW_HI):
                    cands.append(dict(Bc=Bc, a=a, b=b, idx=i))
    cands.sort(key=lambda r: abs(r['Bc'] - 15167.2))
    return cands[0], BS, EV, VEC


def evec_offcrossing(m, Bc, ref_vecs, a, b, offset=15.0):
    """Clean partner eigenvectors evaluated OFF the crossing (Bc - offset), re-identified by
    maximum overlap with the tracked-map evecs at the crossing index. AT degeneracy the a/b
    evecs are an ill-conditioned arbitrary rotation; +-offset agree to <1% (verified in the
    companion files). Verbatim the test_sio_crossing.py method."""
    _, vec_off = m.eigensystem(0.0, Bc - offset)
    ov = np.abs(vec_off @ ref_vecs.T)
    return vec_off[int(np.argmax(ov[:, a]))], vec_off[int(np.argmax(ov[:, b]))]


def xi_at(m, c, ref_vecs, dz):
    """xi(B_c) = |<A|d_z|B>| from the +-offset off-crossing eigenvectors (averaged)."""
    va, vb = evec_offcrossing(m, c['Bc'], ref_vecs, c['a'], c['b'], +15.0)
    va2, vb2 = evec_offcrossing(m, c['Bc'], ref_vecs, c['a'], c['b'], -15.0)
    return 0.5 * (abs(va @ dz @ vb) + abs(va2 @ dz @ vb2))


# ===========================================================================
# Module-level shared setup
# ===========================================================================
print("=" * 78)
print("29SiO+ deferred-number package I -- xi(B)/d_eff, E budget, B_perp x E "
      "cross term, leakage table")
print("=" * 78)

M = build_model()
MUE = M.parameters['muE']            # MHz/(V/cm) = 4.147 D * 0.503412
G_S = M.parameters['g_S']
MU_B = M.parameters['mu_B']          # MHz/G
DZ = op_matrix(M, 'StarkZ')          # d_z (unit-muE axial dipole)
ZX = op_matrix(M, 'ZeemanX')         # transverse electron-spin Zeeman operator
CANON, B_SCAN, EV, VEC = find_canonical(M)
BC = CANON['Bc']
XI_C = xi_at(M, CANON, VEC[CANON['idx']], DZ)
D_EFF_MHZ = XI_C * MUE               # MHz/(V/cm)
D_EFF_DEBYE = XI_C * DEBYE_PER_UNIT_MUE


# ===========================================================================
# Deliverable 1 -- xi(B) and d_eff over the crossing window + gamma error band.
#   GATE: xi(B_c) reproduces 2.660e-4 within 0.5%.
# ===========================================================================
def deliverable1_xi():
    print("\n  --- D1: xi(B) = |<A|d_z|B>| and d_eff over the crossing window ---")
    print(f"  canonical pair A=|N=0,J=1/2,F=0,M_F=0> (+), B=|N=1,J=3/2,F=2,M_F=0> (-); "
          f"B_c = {BC:.1f} G")
    print(f"  muE = {MUE:.9f} MHz/(V/cm) = {DEBYE_PER_UNIT_MUE} D")
    print(f"  {'B[G]':>9} {'xi':>13} {'d_eff[MHz/(V/cm)]':>18} {'d_eff[Debye]':>13}")
    # Tracked-pair sweep (indices a,b carry identity from the low-field map).
    Bwin = np.linspace(SWEEP_LO, SWEEP_HI, SWEEP_N)
    EVw, VECw = zeeman_map(M, Bwin)
    for j, B in enumerate(Bwin):
        va = VECw[j, CANON['a']]
        vb = VECw[j, CANON['b']]
        xi = abs(va @ DZ @ vb)
        print(f"  {B:9.1f} {xi:13.6e} {xi * MUE:18.6e} {xi * DEBYE_PER_UNIT_MUE:13.6e}")
    print(f"  [off-crossing +-15 G at B_c] xi(B_c) = {XI_C:.6e}  "
          f"d_eff = {D_EFF_MHZ:.6e} MHz/(V/cm) = {D_EFF_DEBYE:.6e} Debye")

    # gamma error band: recompute xi(B_c) for the Gamma_SR band.
    print(f"\n  --- D1 gamma error band (Decision D1): xi(B_c) vs Gamma_SR ---")
    print(f"  {'Gamma_SR[MHz]':>13} {'B_c[G]':>9} {'xi(B_c)':>13} {'rel-to-nominal':>15}")
    band = []
    for g in GAMMA_BAND:
        mg = build_model(gamma=g)
        dzg = op_matrix(mg, 'StarkZ')
        cg, _, _, vecg = find_canonical(mg)
        xig = xi_at(mg, cg, vecg[cg['idx']], dzg)
        rel = abs(xig - XI_C) / XI_C
        band.append((g, cg['Bc'], xig, rel))
        print(f"  {g:13.0f} {cg['Bc']:9.1f} {xig:13.6e} {rel:15.3e}")
    band_lo = min(r[2] for r in band)
    band_hi = max(r[2] for r in band)

    rel_err = abs(XI_C - XI_TARGET) / XI_TARGET
    ok = rel_err < 0.005
    record("D1 xi(B_c) reproduces 2.660e-4 within 0.5%", ok,
           f"xi(B_c)={XI_C:.6e} (target {XI_TARGET:.3e}, {100 * rel_err:.4f}%); "
           f"d_eff={D_EFF_MHZ:.4e} MHz/(V/cm) = {D_EFF_DEBYE:.4e} D; "
           f"gamma band xi in [{band_lo:.4e}, {band_hi:.4e}] over Gamma_SR{GAMMA_BAND}")


# ===========================================================================
# Deliverable 2 -- E-field budget for Protocol I pi_x pulses.
#   Omega/2pi = 1/(2 t_p);  E_pi = Omega/(xi*muE).
#   GATE: E_pi(10us) * (xi*muE) * 2 * 10us == 1 within 1e-6 (arithmetic self-check).
# ===========================================================================
def deliverable2_budget():
    print(f"\n  --- D2: Protocol-I E-field budget (tau={TAU_MS} ms, N={N_CYC}, "
          f"T={T_TOTAL_S} s) ---")
    print(f"  pi_x pulse: Omega/2pi = 1/(2 t_p);  E_pi = Omega / (xi*muE);  "
          f"xi*muE = d_eff = {D_EFF_MHZ:.6e} MHz/(V/cm)")
    print(f"  {'t_p[us]':>8} {'Omega/2pi[kHz]':>15} {'E_pi[V/cm]':>13} {'flag':>22}")
    deff = D_EFF_MHZ  # MHz/(V/cm)
    for tp_us in TP_LIST_US:
        tp_s = tp_us * 1e-6
        omega_MHz = 1.0 / (2.0 * tp_s) * 1e-6     # 1/(2 t_p) in Hz -> MHz
        omega_kHz = omega_MHz * 1e3
        E_pi = omega_MHz / deff                    # V/cm
        flag = ""
        if E_pi > EPI_WATERMARK:
            flag = f">{EPI_WATERMARK:.0f} V/cm (trap watermark)"
        print(f"  {tp_us:8.1f} {omega_kHz:15.4f} {E_pi:13.4f} {flag:>22}")

    # Arithmetic self-check at t_p = 10 us: E_pi * (xi*muE) * 2 * t_p == 1.
    tp_s = 10e-6
    omega_MHz = 1.0 / (2.0 * tp_s) * 1e-6
    E_pi = omega_MHz / deff
    # E_pi[V/cm] * deff[MHz/(V/cm)] gives Omega[MHz]; *2*t_p[s] uses Omega in Hz.
    check = (E_pi * deff * 1e6) * 2.0 * tp_s       # = Omega[Hz] * 2 * t_p
    ok = abs(check - 1.0) < 1e-6
    record("D2 E_pi(10us) arithmetic self-check (E_pi*xi*muE*2*10us == 1)", ok,
           f"E_pi(10us)={E_pi:.4f} V/cm; closure E_pi*d_eff*2*t_p = {check:.10f} "
           f"(target 1, dev {abs(check - 1.0):.2e})")


# ===========================================================================
# Deliverable 3 -- B_perp x E_z cross term.
#   VE  = -muE*d_z  (axial Stark operator, MHz/(V/cm))
#   VBx =  g_S*mu_B*zx  (transverse electron-spin Zeeman operator, MHz/G)
#   (a) M(B_x) = <A| (E_z*VE + B_x*VBx) |B> is REAL to machine precision (sigma_x-class). GATE.
#   (b) bilinear d2M/dB_x dE_z (Hz/(mG V/cm)): ZERO to machine precision (M_F-forbidden).
#       Ship the QUADRATIC budget number: relative drive shift at B_x=10 mG and the
#       coefficient d2(drive)/dB_x^2 at E_z=6 V/cm. (reversal-suite budget)
#   (c) pure-B_x pair element at E=0 is zero (parity/M_F exact). GATE, machine precision.
# ===========================================================================
def deliverable3_cross():
    VE = -MUE * DZ                    # MHz/(V/cm)
    VBx = G_S * MU_B * ZX             # MHz/G
    H0 = M.H_function(0.0, BC)        # axial Hamiltonian at B_c
    real_ops = (np.isrealobj(H0) and np.isrealobj(VE) and np.isrealobj(VBx))
    # bare pair vectors (off-crossing, real, gauge-fixed)
    va, vb = evec_offcrossing(M, BC, VEC[CANON['idx']], CANON['a'], CANON['b'], 15.0)

    print(f"\n  --- D3: B_perp x E_z cross term (E_z={E_DRIVE} V/cm; static B_x) ---")
    print(f"  operators real: H0={np.isrealobj(H0)}, VE={np.isrealobj(VE)}, "
          f"VBx=g_S*mu_B*zx={np.isrealobj(VBx)} (any combination -> real H -> no sigma_y)")

    # (a) M(B_x) = <A|(E_z*VE + B_x*VBx)|B> in the bare pair basis; REAL to machine precision.
    print(f"  {'B_x[mG]':>8} {'Re M[MHz]':>15} {'Im M[MHz]':>12} {'Im/|M|':>11}")
    max_imag_ratio = 0.0
    for bx_mG in BX_LIST_MG:
        bx = bx_mG * 1e-3            # G
        Mel = complex(va.conj() @ (E_DRIVE * VE + bx * VBx) @ vb)
        ratio = abs(Mel.imag) / abs(Mel) if abs(Mel) > 0 else 0.0
        max_imag_ratio = max(max_imag_ratio, ratio)
        print(f"  {bx_mG:8.1f} {Mel.real:15.8e} {Mel.imag:12.2e} {ratio:11.2e}")
    a_ok = max_imag_ratio < 1e-10

    # (c) pure-B_x pair element at E=0: <A|VBx|B> exactly zero (parity/M_F).
    pure_bx = abs(va.conj() @ VBx @ vb)
    c_ok = pure_bx < 1e-10
    print(f"  pure-B_x pair element at E=0: |<A|g_S*mu_B*zx|B>| = {pure_bx:.2e} "
          f"(parity/M_F exact zero)")

    # (b) bilinear d2M/dB_x dE_z. Literal operator: linear in each field, no E_x*B_x
    #     product -> identically 0. Confirm by central finite differences.
    hE, hB = 1e-3, 1e-6           # V/cm, G
    def Mop(Ez, Bx):
        return float(va.conj() @ (Ez * VE + Bx * VBx) @ vb)
    bilin_op = (Mop(E_DRIVE + hE, hB) - Mop(E_DRIVE + hE, -hB)
                - Mop(E_DRIVE - hE, hB) + Mop(E_DRIVE - hE, -hB)) / (4 * hE * hB)
    bilin_op_Hz = bilin_op * 1e6 * 1e-3   # MHz/(V/cm.G) -> Hz/(mG.V/cm)

    # Dressed-pair second-order perturbation theory (CONDITIONING-CLEAN). B_x admixes
    # M_F=+-1 spectators into A,B at 1st order; the drive correction comes from those
    # admixtures. Done in the H0(B_c) eigenbasis with EXPLICIT energy denominators,
    # EXCLUDING the near-degenerate pair partner (which is M_F=0 -> VBx element zero by
    # M_F anyway, so the exclusion changes nothing physical but removes the 1/(3.6e-7 MHz)
    # blowup). A full-36-dim re-diagonalization with overlap-tracking is ILL-CONDITIONED
    # near the ~1e-7 MHz pair degeneracy and is NOT used for the shipped number.
    ev0, U0 = np.linalg.eigh(H0)            # U0[:,k] = k-th eigenvector
    ov0 = np.abs(U0.T @ np.array([va, vb]).T)
    ia = int(np.argmax(ov0[:, 0])); ib = int(np.argmax(ov0[:, 1]))
    A0 = U0[:, ia].copy(); B0 = U0[:, ib].copy()
    if A0 @ va < 0: A0 = -A0
    if B0 @ vb < 0: B0 = -B0
    bare = A0 @ VE @ B0                      # bare drive amplitude (unit E)

    def first_order_corr(psi0, E0, skip):
        corr = np.zeros_like(psi0)
        for s in range(M.size):
            if s == skip:
                continue
            denom = E0 - ev0[s]
            if abs(denom) < 1e-6:
                continue                     # guards the degenerate partner (M_F-null anyway)
            corr = corr + ((U0[:, s] @ VBx @ psi0) / denom) * U0[:, s]
        return corr
    Acorr = first_order_corr(A0, ev0[ia], ia)   # per Gauss
    Bcorr = first_order_corr(B0, ev0[ib], ib)
    lin_coef = (Acorr @ VE @ B0) + (A0 @ VE @ Bcorr)    # drive coeff of B_x (per G); M_F-null
    quad_coef = Acorr @ VE @ Bcorr                       # drive coeff of B_x^2 (per G^2)
    # bilinear d2(drive)/dE_z dB_x: drive is linear in E_z, so = lin_coef (per G) -> Hz/(mG.V/cm)
    bilin_dressed_Hz = lin_coef * 1e6 * 1e-3
    rel_shift_10mG = (lin_coef * 10e-3 + quad_coef * (10e-3) ** 2) / bare
    quad_Hz_per_mG2_per_Vcm = quad_coef * 1e6 * 1e-6     # Hz/(V/cm.mG^2)

    print(f"\n  bilinear d2M/dB_x dE_z (literal operator)  = {bilin_op_Hz:.3e} Hz/(mG.V/cm)  "
          f"(M_F-forbidden -> ~0 noise)")
    print(f"  bilinear d2(drive)/dB_x dE_z (dressed PT2)  = {bilin_dressed_Hz:.3e} "
          f"Hz/(mG.V/cm)  (M_F-forbidden -> ~0 noise)")
    print(f"  *** SHIPPED budget number (leading B_x systematic; PT2 admixture, QUADRATIC) ***")
    print(f"      relative drive shift at B_x=10 mG = {rel_shift_10mG:+.3e}")
    print(f"      d2(drive)/dB_x^2 / E_z = {quad_Hz_per_mG2_per_Vcm:.3e} Hz/(V/cm.mG^2)")

    # GATE 3a: reality; GATE 3c: pure-B_x null. (3b is a recorded finding, not a pass/fail
    # number -- the bilinear being a machine-precision null IS the result.) Gate the two
    # machine-precision claims, plus require the bilinear to actually be a null (<1 Hz/(mG.V/cm)).
    bilin_null = abs(bilin_op_Hz) < 1.0 and abs(bilin_dressed_Hz) < 1.0
    ok = a_ok and c_ok and real_ops and bilin_null
    record("D3 cross term: M real (sigma_x), pure-B_x null at E=0, bilinear "
           "machine-precision null", ok,
           f"(a) max Im/|M| = {max_imag_ratio:.2e} (<1e-10: {a_ok}); "
           f"(c) |<A|VBx|B>|@E=0 = {pure_bx:.2e} (<1e-10: {c_ok}); "
           f"(b) bilinear d2M/dB_x dE_z = {bilin_op_Hz:.2e} (op) / {bilin_dressed_Hz:.2e} "
           f"(dressed) Hz/(mG.V/cm) = machine-precision null (M_F-forbidden); "
           f"SHIPPED: rel drive shift at 10 mG = {rel_shift_10mG:+.3e}, "
           f"quad coeff {quad_Hz_per_mG2_per_Vcm:.3e} Hz/(V/cm.mG^2)")
    return dict(rel_shift_10mG=rel_shift_10mG,
                quad=quad_Hz_per_mG2_per_Vcm,
                bilin_op_Hz=bilin_op_Hz, bilin_dressed_Hz=bilin_dressed_Hz,
                pure_bx=pure_bx)


# ===========================================================================
# Deliverable 3b -- the SURVIVING stray cross channel: B_x x E_x (both transverse).
#   D3 established that B_x x E_z is M_F-FORBIDDEN at linear order (axial drive
#   keeps M_F, transverse Zeeman changes it by +-1). The manuscript's fakes-table
#   row budgets the STRAY-field cross term, where the stray E has transverse
#   components too: B_x changes M_F by +1 and E_x by -1 (or vice versa), so the
#   bilinear pair coupling between the two M_F=0 members SURVIVES M_F selection.
#   Effective 2x2 pair coupling at second order (degenerate PT, pair excluded):
#     M_bilin = sum_s [<A|VEx|s><s|VBx|B> + <A|VBx|s><s|VEx|B>] / (E_pair - E_s)
#   per (V/cm x G). GATE: the coupling is REAL (sigma_x-class) -- static fields
#   keep H real, so this is a drive/contrast systematic, never a W-fake.
# ===========================================================================
def deliverable3b_cross_transverse():
    VEx = -MUE * op_matrix(M, 'StarkX')   # MHz/(V/cm), transverse Stark (Task 5)
    VBx = G_S * MU_B * ZX                 # MHz/G
    H0 = M.H_function(0.0, BC)
    ev0, U0 = np.linalg.eigh(H0)
    va, vb = evec_offcrossing(M, BC, VEC[CANON['idx']], CANON['a'], CANON['b'], 15.0)
    ov0 = np.abs(U0.T @ np.array([va, vb]).T)
    ia = int(np.argmax(ov0[:, 0])); ib = int(np.argmax(ov0[:, 1]))
    A0 = U0[:, ia].copy(); B0 = U0[:, ib].copy()
    if A0 @ va < 0: A0 = -A0
    if B0 @ vb < 0: B0 = -B0
    E_pair = 0.5 * (ev0[ia] + ev0[ib])    # pair degenerate at B_c to ~1e-7 MHz

    coef = 0.0 + 0.0j                      # MHz per (V/cm * G)
    for s in range(M.size):
        if s in (ia, ib):
            continue
        denom = E_pair - ev0[s]
        if abs(denom) < 1e-6:
            continue
        S = U0[:, s]
        coef += ((A0 @ VEx @ S) * (S @ VBx @ B0)
                 + (A0 @ VBx @ S) * (S @ VEx @ B0)) / denom
    coef_real = float(np.real(coef))
    imag_ratio = abs(np.imag(coef)) / abs(coef) if abs(coef) > 0 else 0.0
    coef_Hz_mG_Vcm = coef_real * 1e6 * 1e-3        # MHz/(V/cm.G) -> Hz/(mG.V/cm)

    bx, ex = 10.0, 0.1                              # mG, V/cm (stray budget point)
    budget_Hz = coef_Hz_mG_Vcm * bx * ex
    drive_Hz = 3331.9                               # deliberate drive at 6 V/cm (Hz)
    W_HZ = 0.4721

    print(f"\n  --- D3b: surviving stray cross channel B_x x E_x (both transverse) ---")
    print(f"  bilinear pair-coupling coeff d2M/dB_x dE_x = {coef_Hz_mG_Vcm:.4e} Hz/(mG.V/cm)")
    print(f"  Im/|coef| = {imag_ratio:.2e}  (real sigma_x-class -> drive systematic, NOT a W-fake)")
    print(f"  budget at B_x=10 mG, stray E_x=0.1 V/cm: {budget_Hz:.4e} Hz")
    print(f"    vs deliberate drive 3331.9 Hz (6 V/cm): {budget_Hz/drive_Hz:.2e}")
    print(f"    vs W = 0.4721 Hz:                       {budget_Hz/W_HZ:.2e}")

    ok = (imag_ratio < 1e-10) and np.isfinite(coef_Hz_mG_Vcm)
    record("D3b stray B_x x E_x channel (survives M_F selection; REAL sigma_x-class)",
           ok,
           f"coeff = {coef_Hz_mG_Vcm:.3e} Hz/(mG.V/cm); Im/|coef| = {imag_ratio:.1e} "
           f"(<1e-10: real); budget(10 mG x 0.1 V/cm) = {budget_Hz:.3e} Hz "
           f"= {budget_Hz/W_HZ:.1e} x W")
    return dict(coef_Hz_mG_Vcm=coef_Hz_mG_Vcm, budget_Hz=budget_Hz,
                imag_ratio=imag_ratio)


# ===========================================================================
# Deliverable 4 -- B_perp leakage table.
#   At B_c, for both pair members, every spectator in N=0,1 with its spin-resonance
#   matrix element |<spectator| g_S*mu_B*zx |pair member>| (per mG of B_x, kHz/mG)
#   and its detuning (MHz). Sort by |ME|/|detuning|.
#   GATE: the largest leakage-to-detuning ratio among nonzero-ME channels is finite,
#         and the table covers all 16 N=0,1 states.
# ===========================================================================
def deliverable4_leakage():
    VBx = G_S * MU_B * ZX            # MHz/G; per mG -> *1e-3 G; in kHz -> *1e3 -> numerically MHz/G
    # Eigenstates at B_c via DIRECT eigensystem (energy-sorted). A single-point ZeemanMap does
    # NOT preserve the low-field tracked index identity (it re-sorts by energy at that one
    # field), so spectators are identified here by dominant N from the direct B_c eigenstates
    # -- well-conditioned because ONLY the crossing PAIR is degenerate (~1e-7 MHz gap). The two
    # PAIR MEMBERS are taken from the OFF-CROSSING (B_c-15 G) eigenvectors -- the conditioning-
    # clean method of test_sio_crossing.py / test_nsdpv (at B_c the pair is an arbitrary
    # rotation; member A would otherwise pick up ~1.40 kHz/mG of B's electron-spin character).
    ev_c, vec_c = M.eigensystem(0.0, BC)     # vec_c[k] = k-th eigenvector (row), energy-sorted
    domn_bc = np.array([dom_N(M, vec_c[k]) for k in range(M.size)])
    n01 = [k for k in range(M.size) if domn_bc[k] in (0, 1)]
    pa_idx, pb_idx = CANON['a'], CANON['b']
    # Clean pair vectors (off-crossing) and their degenerate energy at B_c.
    va, vb = evec_offcrossing(M, BC, VEC[CANON['idx']], pa_idx, pb_idx, 15.0)
    # Pair eigenstate indices among the direct B_c eigenstates, by overlap with va,vb.
    ip_a = int(np.argmax(np.abs(vec_c @ va)))
    ip_b = int(np.argmax(np.abs(vec_c @ vb)))
    E_pair = 0.5 * (ev_c[ip_a] + ev_c[ip_b])     # true degenerate pair energy (~21449.42 MHz)
    # Zero-field-adiabatic identity of the canonical pair (low-field index carries its label).
    def zf_lab(k):
        v = VEC[0, k]; idx = int(np.argmax(v ** 2))
        return (int(M.q_numbers['N'][idx]), M.q_numbers['J'][idx],
                M.q_numbers['F'][idx], M.q_numbers['M'][idx])
    labA = zf_lab(pa_idx); labB = zf_lab(pb_idx)

    def lab(k):
        v = vec_c[k]
        return (int(dom_N(M, v)), dom_q(M, v, 'J'), dom_q(M, v, 'F'),
                dom_M(M, v), parity_of(M, v))

    rows = []
    for tag, vpm in (('A', va), ('B', vb)):
        for k in n01:
            if k in (ip_a, ip_b):           # exclude both pair members as spectators
                continue
            me_kHz_per_mG = abs(vec_c[k] @ VBx @ vpm)   # MHz/G == kHz/mG numerically
            det_MHz = ev_c[k] - E_pair
            rows.append(dict(pm=tag, k=k, lab=lab(k),
                             me=me_kHz_per_mG, det=det_MHz))

    # Rank by |ME|/|detuning|, considering only nonzero-ME channels (a zero ME carries no
    # leakage regardless of detuning; a 0/0 row is not a real channel).
    def ratio(r):
        return r['me'] / abs(r['det']) if abs(r['det']) > 1e-12 else np.inf
    nz = [r for r in rows if r['me'] > 1e-9]
    nz.sort(key=ratio, reverse=True)

    print(f"\n  --- D4: B_perp leakage table at B_c = {BC:.1f} G "
          f"(operator g_S*mu_B*zx, per mG of B_x) ---")
    print(f"  pair members (zero-field-adiabatic): "
          f"A=|N={labA[0]},J={labA[1]},F={labA[2]:.0f},M_F={labA[3]:+.0f}> (+), "
          f"B=|N={labB[0]},J={labB[1]},F={labB[2]:.0f},M_F={labB[3]:+.0f}> (-); "
          f"clean off-crossing vectors; E_pair={E_pair:.3f} MHz")
    print(f"  {'PM':>2} {'spectator (N,J,F,M_F,par)':>26} {'|ME|[kHz/mG]':>13} "
          f"{'detuning[MHz]':>14} {'|ME|/|det|':>12}")
    for r in nz:
        N, J, F, Mf, p = r['lab']
        ls = f"({N},{J:.1f},{F:.0f},{Mf:+.0f},{p:+d})"
        print(f"  {r['pm']:>2} {ls:>26} {r['me']:13.5f} {r['det']:14.4f} "
              f"{ratio(r):12.4e}")

    # Coverage: union of spectators in the table + the two pair members must be all 16 N=0,1.
    covered = set(r['k'] for r in rows) | {ip_a, ip_b}
    cover_ok = (len(covered) == 16) and (len(n01) == 16)
    max_ratio = max(ratio(r) for r in nz) if nz else np.inf
    ratio_finite = np.isfinite(max_ratio)

    ok = cover_ok and ratio_finite and len(nz) >= 1
    record("D4 leakage table (16 N=0,1 states covered, max finite |ME|/|det|)", ok,
           f"{len(nz)} nonzero-ME leakage channels of {len(rows)} "
           f"(2 members x 14 N=0,1 spectators); covers {len(covered)}/16 N=0,1 states "
           f"(cover_ok={cover_ok}); largest |ME|/|det| = {max_ratio:.4e} "
           f"kHz/mG/MHz (finite={ratio_finite}); top channel: "
           f"PM {nz[0]['pm']} -> {nz[0]['lab']}, |ME|={nz[0]['me']:.4f} kHz/mG, "
           f"det={nz[0]['det']:.3f} MHz")
    return nz


# ===========================================================================
# Deliverable 5 -- final summary table of every shipped number with units.
# ===========================================================================
def deliverable5_summary(d3, leak_nz):
    print("\n" + "=" * 78)
    print("  SHIPPED-NUMBER SUMMARY (29SiO+ deferred-number package I)")
    print("=" * 78)
    tp10 = 10e-6
    E_pi_10 = (1.0 / (2.0 * tp10) * 1e-6) / D_EFF_MHZ
    top3 = leak_nz[:3]
    lines = [
        ("B_c (canonical Delta M_F=0 crossing)", f"{BC:.1f}", "G"),
        ("xi(B_c) = |<A|d_z|B>| (unit-muE dipole fraction)", f"{XI_C:.6e}", "dimensionless"),
        ("d_eff = xi*muE", f"{D_EFF_MHZ:.6e}", "MHz/(V/cm)"),
        ("d_eff = xi*muE", f"{D_EFF_DEBYE:.6e}", "Debye"),
        ("Omega/2pi for pi_x at t_p=10 us", f"{1.0 / (2.0 * tp10) / 1e3:.3f}", "kHz"),
        ("E_pi (t_p=10 us)", f"{E_pi_10:.4f}", "V/cm"),
        ("B_perp x E_z bilinear d2M/dB_x dE_z (operator)", f"{d3['bilin_op_Hz']:.3e}",
         "Hz/(mG.V/cm) [machine-precision null]"),
        ("B_perp x E_z bilinear (dressed-state)", f"{d3['bilin_dressed_Hz']:.3e}",
         "Hz/(mG.V/cm) [machine-precision null]"),
        ("relative drive shift at B_x=10 mG (quadratic budget)",
         f"{d3['rel_shift_10mG']:+.3e}", "dimensionless"),
        ("d2(drive)/dB_x^2 / E_z (quadratic coeff)", f"{d3['quad']:.3e}",
         "Hz/(V/cm.mG^2)"),
        ("pure-B_x pair element at E=0", f"{d3['pure_bx']:.2e}",
         "(parity/M_F exact null)"),
    ]
    for name, val, unit in lines:
        print(f"  {name:<52} {val:>16}  {unit}")
    print(f"  top-3 B_perp leakage channels (|ME| kHz/mG, detuning MHz, ratio):")
    for r in top3:
        N, J, F, Mf, p = r['lab']
        print(f"    PM {r['pm']} -> (N={N},J={J:.1f},F={F:.0f},M_F={Mf:+.0f},par {p:+d}): "
              f"|ME|={r['me']:.5f} kHz/mG, det={r['det']:.3f} MHz, "
              f"ratio={r['me'] / abs(r['det']) if abs(r['det']) > 1e-12 else float('inf'):.4e}")


if __name__ == "__main__":
    deliverable1_xi()
    deliverable2_budget()
    d3 = deliverable3_cross()
    d3b = deliverable3b_cross_transverse()
    leak_nz = deliverable4_leakage()
    deliverable5_summary(d3, leak_nz)
    print(f"  B_x x E_x stray coeff (surviving channel)            "
          f"{d3b['coef_Hz_mG_Vcm']:>16.3e}  Hz/(mG.V/cm)")
    print(f"  B_x x E_x budget (10 mG x 0.1 V/cm)                  "
          f"{d3b['budget_Hz']:>16.3e}  Hz  ({d3b['budget_Hz']/0.4721:.1e} x W)")

    n_pass = sum(1 for _, p, _ in _results if p)
    n_fail = sum(1 for _, p, _ in _results if not p)
    print("\n" + "=" * 78)
    print(f"{'OK' if n_fail == 0 else 'FAIL'}: {n_pass} passed, {n_fail} failed")
    sys.exit(0 if n_fail == 0 else 1)
