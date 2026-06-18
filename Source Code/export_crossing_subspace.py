"""Frozen-operator export of the 29SiO+ X2Sigma+ crossing subspace (Task 9).

Builds the 36-state 29SiO+ model (N_list=[0,1,2], I_nuclei=[0,1/2], boson,
P_values=[1/2]) and freezes the operators the downstream two-level / three-level
toy model needs into ONE .npz file. The toy repo never imports Molecule-Structure;
it loads this npz and re-derives the pair reduction from the frozen matrices.

    conda run -n Structure python "Source Code/export_crossing_subspace.py"

Exit 0 iff every exporter gate passes; exit 1 on any [FAIL]. The npz is written
only after all gates pass.

What is exported (all energies in MHz, B in Gauss, E in V/cm; see the `units`
JSON manifest inside the npz):

  H0     (36x36 complex)  zero-field Hamiltonian = H_function(0, 0).
  V_B    (36x36 complex)  axial Zeeman coupling per Gauss, EXTRACTED from the
                          builder's own field map as (H_function(0,B0)-H0)/B0.
                          This is exactly the operator the spectrum uses
                          (g_S*mu_B*ZeemanZ - g_N*mu_N*ZeemanIZ for SiO+, which
                          has no g_l term). Verified field-linear to machine
                          precision over 14000-16000 G, so the single-slope
                          freeze is exact, not a finite-difference approximation.
  V_E    (36x36 complex)  axial Stark per V/cm = H_function(1,0) - H0
                          = -muE*StarkZ matrix (the trusted axial dipole drive).
  V_Ex   (36x36 complex)  transverse (lab-x, p=+-1) Stark per V/cm
                          = -muE * StarkX matrix (analogue of V_E for the
                          transverse field; StarkX registered in bBJ_even_X).
  V_Bx   (36x36 complex)  transverse (lab-x, p=+-1) Zeeman per Gauss
                          = g_S*mu_B * ZeemanX matrix (the transverse spin
                          Zeeman; SiO+ has no g_l, so this is the full transverse
                          magnetic coupling at the level the axial term is wired).
  C      (36x36 complex)  the NSD-PV anapole operator C = (i/I)(n_hat x S).I,
                          dimensionless. H_PV = kappaPrime * W_A * C. Phase pinned
                          to Karthein Eq. (B1): <+|C|->/i > 0 for bra = even-N.
  d_p    (3 x 36 x 36)    spherical dipole components p = -1, 0, +1 (unit muE;
                          multiply by muE_MHz_per_Vcm for MHz/(V/cm)). d_p[1]
                          (p=0) == StarkZ, and (d_p[0]-d_p[2])/sqrt2 == StarkX.
  parity (36,)           diagonal of Parity_mat (+-1 per zero-field eigenstate;
                          here the basis is the zero-field eigenbasis so parity
                          is a good quantum number on the diagonal).
  qn     structured       N, J, F, M per basis state (the coupled (N,J,F,M) basis
                          the operators above are written in).
  pair_idx (2,)          [A, B] basis indices of the canonical Delta M_F=0 (0,0)
                          crossing pair nearest B_c, in the SAME basis ordering as
                          all operators. A is the even-N (+) partner, B the odd-N.
                          NOTE: these are basis indices, not tracked-map indices;
                          the toy repo reconstructs the field-dressed pair from
                          the operators (the +-15 G eigenvector method), so the
                          pair_idx is provided only as a labeling convenience /
                          cross-check anchor.
  Bc_G   scalar          15167.2 (the canonical crossing field, Gauss).
  ref_evals (3 x 36)     eigenvalues (MHz, ascending) at the three reference
                          points (E,B) = (0,15000), (6,15167.2), (0,15300), for
                          the toy repo's round-trip gate.
  units  JSON string     the unit / convention manifest.

Exporter gates (run before writing; the npz is the deliverable, the gates guard
it):

  GATE A  reconstruction: H0 + E*V_E + B*V_B reproduces H_function(E,B) at the 3
          reference points to 1e-9 RELATIVE (max|rebuilt-true| / max|true|).
          Catches any nonlinearity, basis-ordering slip, or operator mismatch.
  GATE B  V_B linearity: (H_function(0,B)-H0)/B is B-independent over 14000-16000
          G to machine precision (this model's Zeeman is exactly linear).
  GATE C  d_p / StarkX consistency: d_p[1]==StarkZ and (d_p[0]-d_p[2])/sqrt2
          ==StarkX to machine precision; V_E == -muE*d_p[1] to 1e-9.
  GATE D  C phase pin: at the canonical (0,0) crossing, field-dressed <+|C|->/i
          > 0 and |<+|C|->| within 3% of 1/sqrt3 (Karthein Eq. B1; the +-15 G
          off-crossing eigenvector method, verbatim from test_nsdpv_operator.py).
"""
import json
import os
import sys
from functools import partial

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))

from Energy_Levels import MoleculeLevels
from molecule_parameters import get_molecule_params
import hamiltonian_builders as ham
import matrix_elements as me
from sympy.physics.wigner import wigner_3j, wigner_6j

# ---------------------------------------------------------------------------
# Constants / reference points
# ---------------------------------------------------------------------------
def _resolve_out_path():
    """Cross-machine: write to whichever three-level-toy-model/data exists."""
    cands = [
        "/Users/arianjadbabaie/Code-Local/three-level-toy-model/data/sio_crossing_subspace.npz",
        os.path.expanduser("~/Code-Local/three-level-toy-model/data/sio_crossing_subspace.npz"),
        os.path.expanduser("~/Code/three-level-toy-model/data/sio_crossing_subspace.npz"),
        "C:/Users/Arian/Code/three-level-toy-model/data/sio_crossing_subspace.npz",
    ]
    for c in cands:
        if os.path.isdir(os.path.dirname(c)):
            return c
    return cands[0]


OUT_PATH = _resolve_out_path()
BC_G = 15167.2                # canonical (0,0) crossing field, Gauss
E_DRIVE = 6.0                 # V/cm, Karthein axial drive field
KPRIME = 0.05                 # assumed mixing coefficient (Karthein Fig. 3)
WA_HZ = 16.0                  # W_A/2pi in Hz (Karthein assumed electronic m.e.)
INV_SQRT3 = 1.0 / np.sqrt(3.0)
# Three reference points (E [V/cm], B [G]) for the round-trip gate.
REF_POINTS = [(0.0, 15000.0), (6.0, 15167.2), (0.0, 15300.0)]

_results = []


def record(name, passed, detail):
    _results.append((name, bool(passed), detail))
    tag = "[PASS]" if passed else "[FAIL]"
    print(f"{tag} {name}: {detail}")


# ---------------------------------------------------------------------------
# Model + helpers (patterns lifted verbatim from the test_sio_* gate suite so
# the frozen operators are the SAME ones those gates validated).
# ---------------------------------------------------------------------------
def build_model():
    p = dict(get_molecule_params('SiO+', 'X', '0', 'boson'))
    return MoleculeLevels.initialize_state(
        'SiO+', 'X', 0, N_list=[0, 1, 2], fermion_or_boson='boson',
        I_nuclei=[0, 1 / 2], P_values=[1 / 2], params=p)


def _dipole_p_me(p, K0, N0, J0, F0, M0, K1, N1, J1, F1, M1, S=1 / 2, I=1 / 2):
    """Raw spherical dipole component T^1_p(d), unit muE. Verbatim from
    test_sio_stark_tdm.py: d_p[0]==StarkZ, (d_p[-1]-d_p[+1])/sqrt2==StarkX."""
    if int(round(K0)) != int(round(K1)):
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
    """d_p matrices over M's own (N,J,F,M) basis, unit muE."""
    qn = M.q_numbers
    out = {}
    for p in (-1, 0, 1):
        out[p] = np.array(ham.build_operator(qn, qn, partial(_dipole_p_me, p)))
    out['z'] = np.array(ham.build_operator(qn, qn, M.matrix_elements['StarkZ']))
    out['x'] = np.array(ham.build_operator(qn, qn, M.matrix_elements['StarkX']))
    return out


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
    """Clean partner eigenvectors OFF the crossing (Bc - offset), re-identified
    by max overlap with the tracked-map evecs. Verbatim from test_sio_crossing.py
    / test_nsdpv_operator.py."""
    _, vec_off = m.eigensystem(0.0, Bc - offset)
    ov = np.abs(vec_off @ ref_vecs.T)
    ka = int(np.argmax(ov[:, a]))
    kb = int(np.argmax(ov[:, b]))
    return vec_off[ka], vec_off[kb]


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
                                pa=parity_of(m, VEC[i, a]),
                                pb=parity_of(m, VEC[i, b]),
                                Ma=dom_M(m, VEC[i, a]),
                                Mb=dom_M(m, VEC[i, b])))
    out.sort(key=lambda r: r['Bc'])
    return [c for c in out if c['pa'] != c['pb']]


# ===========================================================================
# Build operators
# ===========================================================================
print("=" * 78)
print("29SiO+ crossing-subspace exporter (Task 9) -- building 36-state model")
print("=" * 78)

M = build_model()
SIZE = M.size
muE = M.parameters['muE']                  # MHz/(V/cm)
gS = M.parameters['g_S']
mu_B = M.parameters['mu_B']                 # MHz/G
gS_muB = gS * mu_B                          # MHz/G

# --- axial operators, extracted from the builder's own field map ---
H0 = M.H_function(0.0, 0.0)                 # zero-field Hamiltonian, MHz
V_E = M.H_function(1.0, 0.0) - H0           # axial Stark per V/cm  (= -muE*StarkZ)
B0_lin = 10000.0                            # any B in the linear regime
V_B = (M.H_function(0.0, B0_lin) - H0) / B0_lin   # axial Zeeman per Gauss

# --- transverse operators, built from the registered StarkX / ZeemanX MEs ---
d = build_d(M)
StarkX_mat = d['x']                         # unit-muE transverse dipole
ZeemanX_mat = np.array(ham.build_operator(M.q_numbers, M.q_numbers,
                                          M.matrix_elements['ZeemanX']))
V_Ex = -muE * StarkX_mat                    # transverse Stark per V/cm
V_Bx = gS_muB * ZeemanX_mat                 # transverse Zeeman per Gauss
# Cast transverse ops to complex to match the others (they are real-valued).
V_Ex = V_Ex.astype(complex)
V_Bx = V_Bx.astype(complex)

# --- NSD-PV operator (dimensionless; H_PV = kappaPrime * W_A * C) ---
C = M.NSDPV_operator()                      # complex 36x36

# --- d_p stack (p = -1, 0, +1), unit muE ---
d_p = np.array([d[-1], d[0], d[1]])         # (3, 36, 36)

# --- bare nuclear spin operators I_x, I_y, I_z (29Si, I=1/2; dimensionless) ---
# Built in the DECOUPLED product basis (1_rot (x) 1_S (x) I^(1/2), trivial) and
# transformed to the coupled (N,J,F,M) basis via the audited b_decoupled change of
# basis. The I_z so built equals ZeemanIZ_bBJ to machine precision (the C3b proof
# in test_sio_conventions.py); I_x, I_y are the identical construction. These give
# the toy repo the nuclear-spin vector the modulated-anapole (29Si NMR) route needs.
_dq = M.alt_q_numbers['decoupled']
_U = np.array(M.library.basis_changers['b_decoupled'](M.q_numbers, _dq), dtype=float)
_K = np.array(_dq['K']); _N = np.array(_dq['N']); _MN = np.array(_dq['M_N'])
_MS = np.array(_dq['M_S']); _MI = np.array(_dq['M_I'])
_Izd = np.diag(_MI.astype(complex))
_Ixd = np.zeros((SIZE, SIZE), dtype=complex)
_Iyd = np.zeros((SIZE, SIZE), dtype=complex)
for _i in range(SIZE):
    for _j in range(SIZE):
        if (_K[_i] == _K[_j] and _N[_i] == _N[_j] and _MN[_i] == _MN[_j]
                and _MS[_i] == _MS[_j] and abs(_MI[_i] - _MI[_j]) == 1):
            _Ixd[_i, _j] = 0.5
            _Iyd[_i, _j] = -0.5j if _MI[_i] > _MI[_j] else 0.5j
I_x = _U.T @ _Ixd @ _U
I_y = _U.T @ _Iyd @ _U
I_z = _U.T @ _Izd @ _U

# --- parity diagonal ---
parity = np.real(np.round(np.diag(np.array(M.Parity_mat)))).astype(int)

# --- quantum numbers (coupled basis the operators are written in) ---
qn = np.array(
    list(zip(M.q_numbers['N'], M.q_numbers['J'], M.q_numbers['F'],
             M.q_numbers['M'])),
    dtype=[('N', 'f8'), ('J', 'f8'), ('F', 'f8'), ('M', 'f8')])

# --- pair indices: canonical (0,0) crossing nearest B_c, as BASIS indices ---
B_SCAN = np.linspace(1e-6, 16000.0, 2401)
EV, VEC = zeeman_map(M, B_SCAN)
DOMN0 = np.array([dom_N(M, VEC[0, k]) for k in range(SIZE)])
N0_IDX = [k for k in range(SIZE) if DOMN0[k] == 0]
N1_IDX = [k for k in range(SIZE) if DOMN0[k] == 1]
OPP = find_opp_crossings(M, B_SCAN, EV, VEC, N0_IDX, N1_IDX)
WIN = [c for c in OPP if 15050.0 <= c['Bc'] <= 15300.0 and c['Ma'] == c['Mb']]
CANON = min([c for c in WIN if c['Ma'] == 0],
            key=lambda c: abs(c['Bc'] - BC_G)) if WIN else None
if CANON is None:
    print("[FATAL] no canonical (0,0) crossing found -- aborting export")
    sys.exit(1)

# Orient and resolve the field-dressed pair, then map back to dominant basis
# indices of A (even-N, +) and B (odd-N, -) for the pair_idx anchor.
va, vb = evec_offcrossing(M, CANON['Bc'], VEC[CANON['idx']], CANON['a'],
                          CANON['b'], offset=15.0)
if dom_N(M, va) != 0:
    va, vb = vb, va
A_basis = int(np.argmax(va ** 2))
B_basis = int(np.argmax(vb ** 2))
pair_idx = np.array([A_basis, B_basis], dtype=int)

# --- reference eigenvalues at the 3 round-trip reference points ---
ref_evals = np.array([np.linalg.eigvalsh(M.H_function(E, B))
                      for (E, B) in REF_POINTS])   # (3, 36), ascending

# --- units / convention manifest ---
units = json.dumps({
    "energy": "MHz",
    "B": "G",
    "E": "V/cm",
    "muE_MHz_per_Vcm": round(float(muE), 4),
    "gS_muB_MHz_per_G": round(float(gS_muB), 4),
    "phase": "Karthein Eq. B1: <+|C|->/i > 0, bra=even-N",
    "HPV": ("kappaPrime*W_A*C with W_A/2pi=16 Hz for 29SiO+ "
            "(kappaPrime=0.05 assumed); C dimensionless"),
    "I_ops": ("I_x,I_y,I_z bare 29Si nuclear spin (I=1/2), dimensionless, in the "
              "coupled (N,J,F,M) basis; I_z == ZeemanIZ_bBJ to machine precision"),
    "Bc_G": BC_G,
    "ref_points_E_B": REF_POINTS,
})


# ===========================================================================
# Exporter gates
# ===========================================================================
def gate_A_reconstruction():
    worst = 0.0
    for (E, B) in REF_POINTS:
        true = M.H_function(E, B)
        rebuilt = H0 + E * V_E + B * V_B
        denom = np.max(np.abs(true))
        rel = np.max(np.abs(rebuilt - true)) / denom
        worst = max(worst, rel)
    record("GATE A reconstruction (H0+E*V_E+B*V_B == H_function, 1e-9 rel)",
           worst < 1e-9,
           f"max relative deviation over 3 ref points = {worst:.2e} (tol 1e-9)")
    return worst


def gate_B_linearity():
    worst = 0.0
    for B in np.linspace(14000.0, 16000.0, 11):
        vb_local = (M.H_function(0.0, B) - H0) / B
        worst = max(worst, np.max(np.abs(vb_local - V_B)))
    record("GATE B V_B linearity (per-Gauss slope B-independent over 14-16 kG)",
           worst < 1e-9,
           f"max|V_B(B) - V_B| over 14000-16000 G = {worst:.2e} (machine "
           f"precision; Zeeman exactly linear)")
    return worst


def gate_C_dp_consistency():
    dev_dz = float(np.max(np.abs(d_p[1] - d['z'])))
    dev_dx = float(np.max(np.abs(d['x'] - (d_p[0] - d_p[2]) / np.sqrt(2))))
    dev_VE = float(np.max(np.abs(V_E - (-muE * d_p[1]))))
    ok = dev_dz < 1e-12 and dev_dx < 1e-12 and dev_VE < 1e-9
    record("GATE C d_p / StarkX consistency", ok,
           f"d_p[0]==StarkZ dev={dev_dz:.1e}; (d_-1-d_+1)/sqrt2==StarkX "
           f"dev={dev_dx:.1e}; V_E==-muE*d_p[0] dev={dev_VE:.1e}")


def gate_D_phase_pin():
    val_p, va2, vb2 = _dressed(offset=15.0)
    val_m, _, _ = _dressed(offset=-15.0)
    mag = 0.5 * (abs(val_p) + abs(val_m))
    coeff = (val_p / 1j).real                # <+|C|->/i ; B1 demands > 0
    rel = abs(mag - INV_SQRT3) / INV_SQRT3
    W_hz = KPRIME * WA_HZ * mag
    ok = coeff > 0 and rel < 0.03
    record("GATE D C phase pin (Karthein Eq. B1, <+|C|->/i > 0)", ok,
           f"<+|C|->/i = {coeff:+.5f} (B1: >0 -> {coeff > 0}); dressed "
           f"|<C>|={mag:.5f} vs 1/sqrt3={INV_SQRT3:.5f} ({100 * rel:.2f}%, "
           f"tol 3%); W/2pi = kprime*W_A*|<C>| = {W_hz:.4f} Hz")


def _dressed(offset):
    va2, vb2 = evec_offcrossing(M, CANON['Bc'], VEC[CANON['idx']], CANON['a'],
                                CANON['b'], offset=offset)
    if dom_N(M, va2) != 0:
        va2, vb2 = vb2, va2
    return va2 @ C @ vb2, va2, vb2


def gate_E_nuclear_ops():
    ZIZ = np.array(ham.build_operator(M.q_numbers, M.q_numbers,
                                      me.ZeemanIZ_bBJ), dtype=complex)
    d_iz = float(np.max(np.abs(I_z - ZIZ)))
    d_comm = float(np.max(np.abs(I_x @ I_y - I_y @ I_x - 1j * I_z)))
    d_i2 = float(np.max(np.abs(I_x @ I_x + I_y @ I_y + I_z @ I_z
                               - 0.75 * np.eye(SIZE))))
    ok = d_iz < 1e-12 and d_comm < 1e-12 and d_i2 < 1e-12
    record("GATE E nuclear ops (I_z == ZeemanIZ; [Ix,Iy]=iIz; I^2 = 3/4)", ok,
           f"I_z vs ZeemanIZ_bBJ = {d_iz:.1e}; [Ix,Iy]-iIz = {d_comm:.1e}; "
           f"I^2 - 3/4 = {d_i2:.1e}")


if __name__ == "__main__":
    print(f"\n  model: size={SIZE}, muE={muE:.6f} MHz/(V/cm), "
          f"g_S*mu_B={gS_muB:.6f} MHz/G")
    print(f"  canonical (0,0) crossing: Bc={CANON['Bc']:.2f} G "
          f"(pinned export Bc_G={BC_G}); pair basis idx A={A_basis} (N="
          f"{qn['N'][A_basis]:.0f}), B={B_basis} (N={qn['N'][B_basis]:.0f})\n")

    gate_A_reconstruction()
    gate_B_linearity()
    gate_C_dp_consistency()
    gate_D_phase_pin()
    gate_E_nuclear_ops()

    n_pass = sum(1 for _, p, _ in _results if p)
    n_fail = sum(1 for _, p, _ in _results if not p)
    print("\n" + "=" * 78)
    print(f"{'OK' if n_fail == 0 else 'FAIL'}: {n_pass} passed, {n_fail} failed")

    if n_fail:
        print("Gates failed -- NOT writing npz.")
        sys.exit(1)

    os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
    np.savez_compressed(
        OUT_PATH,
        H0=np.asarray(H0, dtype=complex),
        V_B=np.asarray(V_B, dtype=complex),
        V_E=np.asarray(V_E, dtype=complex),
        V_Ex=np.asarray(V_Ex, dtype=complex),
        V_Bx=np.asarray(V_Bx, dtype=complex),
        C=np.asarray(C, dtype=complex),
        I_x=np.asarray(I_x, dtype=complex),
        I_y=np.asarray(I_y, dtype=complex),
        I_z=np.asarray(I_z, dtype=complex),
        d_p=np.asarray(d_p, dtype=complex),
        parity=parity,
        qn=qn,
        pair_idx=pair_idx,
        Bc_G=np.float64(BC_G),
        ref_evals=ref_evals,
        units=units,
    )
    size_mb = os.path.getsize(OUT_PATH) / 1e6
    print(f"\nWrote {OUT_PATH}  ({size_mb:.3f} MB)")
    sys.exit(0)
