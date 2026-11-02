"""Rigorous zero-field validation of TransitionDipole_even_aBJ_noM.

Independent (non-circular) reimplementation of the thesis transition-dipole
chain in Hund's case (a) with single-nucleus hyperfine, built from scratch on
sympy Wigner symbols following:

  - Wigner-Eckart on F                : thesis Eq. A.35
  - spectator theorem (dipole on J,
    nuclear spin I a spectator)       : thesis Eq. A.39
  - case-(a) rotational reduced-D     : thesis Eq. A.54 / A.47
  - line strength / orientation avg   : thesis Eq. 3.1-3.4
        S            = sum_{Mg,Me,p} |<g,i|T1_p(d)|e,j>|^2
        gamma~_{g,e} = S / (2 F_e + 1)
        r~           = S / sum_lower S        (sum_lower r~ = 1)

The reference math does NOT import matrix_elements; the codebase functions
(TransitionDipole_even_aBJ_noM, TDM_p_even_aBJ) are the objects under test.
Basis quantum numbers are pulled from the codebase only to enumerate states.

Pass criteria (absolute, convention-independent):
  A  independent M-resolved ME  == codebase TDM_p_even_aBJ          (< 1e-10)
  B  S_ref / |noM_code|^2 == constant (2 F_e + 1) over all pairs    (< 1e-10)
  C  r~ from _noM == r~ from S_ref, and sum_lower r~ == 1           (< 1e-10)
"""
import warnings; warnings.filterwarnings("ignore")
from functools import lru_cache
import numpy as np
from config_path import add_to_sys_path; add_to_sys_path()
from sympy.physics.wigner import wigner_3j as _w3, wigner_6j as _w6
from Energy_Levels import MoleculeLevels
import matrix_elements as me   # functions UNDER TEST only


@lru_cache(maxsize=None)
def w3(a, b, c, d, e, f):
    return float(_w3(a, b, c, d, e, f))


@lru_cache(maxsize=None)
def w6(a, b, c, d, e, f):
    return float(_w6(a, b, c, d, e, f))


def tdm_ref(p, K0, Sig0, P0, J0, F0, M0, K1, Sig1, P1, J1, F1, M1, I=0.5):
    """Independent <g; J0 F0 M0 | T^1_p(d) | e; J1 F1 M1> (electronic RME = 1).

    g = double-0 (ground / lower), e = double-1 (excited / upper), matching
    the codebase argument convention.  Sigma and electronic factor conserved;
    dipole acts on the rotational part, nuclear spin I is a spectator.
    """
    if Sig0 != Sig1:
        return 0.0
    # Wigner-Eckart on total F (thesis Eq. A.35)
    we = (-1.0) ** (F0 - M0) * w3(F0, 1, F1, -M0, p, M1)
    if we == 0.0:
        return 0.0
    # spectator theorem: dipole on J, I spectator (thesis Eq. A.39)
    sp = ((-1.0) ** (F1 + J0 + I + 1)
          * np.sqrt((2 * F0 + 1) * (2 * F1 + 1))
          * w6(J1, F1, I, F0, J0, 1))
    if sp == 0.0:
        return 0.0
    # case-(a) rotational reduced-D, summed over molecule-frame q
    # (thesis Eq. A.54 / A.47); electronic <Lambda||T^1_q(d)||Lambda'> = 1
    rot = ((-1.0) ** (J0 - P0) * np.sqrt((2 * J0 + 1) * (2 * J1 + 1))
           * sum(w3(J0, 1, J1, -P0, q, P1) for q in (-1, 0, 1)))
    return we * sp * rot


def get_aBJ_qnumbers(state):
    """aBJ quantum-number dict the codebase TDM path actually feeds."""
    if "a" in state.hunds_case:
        return state.q_numbers
    return state.alt_q_numbers["aBJ"]


def build(es, Nl, Pv):
    return MoleculeLevels.initialize_state(
        molecule_name="RaF", elec_state=es, vib_state=0, N_list=Nl,
        fermion_or_boson="boson", M_sublevels="none", I_nuclei=[0, 1 / 2],
        isotope=226, round=8, P_values=Pv)


def main():
    I = 0.5  # 19F nuclear spin (226RaF: 226Ra I=0, only F nucleus)
    X = build("X", [0, 1], [1 / 2])
    A = build("A", [1], [1 / 2, 3 / 2])
    X.eigensystem(0, 1e-6)
    A.eigensystem(0, 1e-6)
    qg = get_aBJ_qnumbers(X)
    qe = get_aBJ_qnumbers(A)
    keys = ("K", "Sigma", "P", "J", "F")
    print("ground aBJ keys:", list(qg.keys()))
    print("excited aBJ keys:", list(qe.keys()))
    ng = len(qg[list(qg.keys())[0]])
    ne = len(qe[list(qe.keys())[0]])
    print(f"basis sizes  ng={ng}  ne={ne}")

    def row(qd, i):
        return {k: float(qd[k][i]) for k in keys}

    maxA = 0.0
    S_ref = np.zeros((ng, ne))
    noM2 = np.zeros((ng, ne))
    ratios = []
    for j in range(ne):
        e = row(qe, j)
        for i in range(ng):
            g = row(qg, i)
            # ---- criterion A: elementwise vs codebase TDM_p_even_aBJ ----
            Fg, Fe = g["F"], e["F"]
            for p in (-1, 0, 1):
                for Mg in np.arange(-Fg, Fg + 1):
                    Me = Mg - p
                    if abs(Me) > Fe:
                        continue
                    ref = tdm_ref(p, g["K"], g["Sigma"], g["P"], g["J"], Fg,
                                  Mg, e["K"], e["Sigma"], e["P"], e["J"], Fe,
                                  Me, I=I)
                    cod = me.TDM_p_even_aBJ(
                        p, [-1, 0, 1], g["K"], g["Sigma"], g["P"], g["J"], Fg,
                        Mg, e["K"], e["Sigma"], e["P"], e["J"], Fe, Me,
                        S=0.5, I=I)
                    maxA = max(maxA, abs(abs(ref) - abs(float(cod))))
            # ---- S_ref = sum_{Mg,Me,p} |ME_ref|^2 ----
            s = 0.0
            for p in (-1, 0, 1):
                for Mg in np.arange(-Fg, Fg + 1):
                    Me = Mg - p
                    if abs(Me) > Fe:
                        continue
                    s += tdm_ref(p, g["K"], g["Sigma"], g["P"], g["J"], Fg,
                                 Mg, e["K"], e["Sigma"], e["P"], e["J"], Fe,
                                 Me, I=I) ** 2
            S_ref[i, j] = s
            # ---- codebase _noM (M args unused inside _noM) ----
            v = me.TransitionDipole_even_aBJ_noM(
                g["K"], g["Sigma"], g["P"], g["J"], Fg, 0.0,
                e["K"], e["Sigma"], e["P"], e["J"], Fe, 0.0, S=0.5, I=I)
            noM2[i, j] = float(v) ** 2
            if noM2[i, j] > 1e-14 and S_ref[i, j] > 1e-14:
                ratios.append((S_ref[i, j] / noM2[i, j], 2 * Fe + 1))

    print(f"\n[A] max | |ME_ref| - |TDM_p_even_aBJ| | = {maxA:.3e}")

    # ---- criterion B: S_ref / |noM|^2 == 2 F_e + 1 ----
    ratios = np.array(ratios)
    devB = np.abs(ratios[:, 0] - ratios[:, 1]).max() if len(ratios) else np.nan
    print(f"[B] pairs compared = {len(ratios)} ; "
          f"max | S_ref/|noM|^2 - (2Fe+1) | = {devB:.3e} ; "
          f"ratio range [{ratios[:,0].min():.4f}, {ratios[:,0].max():.4f}]")

    # ---- criterion C: per-excited normalized branching + closure ----
    def pernorm_cols(M):
        cs = M.sum(axis=0, keepdims=True)
        return np.divide(M, cs, out=np.zeros_like(M), where=cs > 1e-30), cs

    r_ref, cs_ref = pernorm_cols(S_ref)
    r_noM, cs_noM = pernorm_cols(noM2)
    decaying = (cs_ref.ravel() > 1e-12)
    maxC = np.abs(r_ref[:, decaying] - r_noM[:, decaying]).max()
    closure = r_ref[:, decaying].sum(axis=0)
    print(f"[C] excited states with decay = {int(decaying.sum())}/{ne} ; "
          f"max |r~_ref - r~_noM| = {maxC:.3e} ; "
          f"closure sum_lower r~ in [{closure.min():.10f}, {closure.max():.10f}]")

    tol = 1e-9
    okA = maxA < tol
    okB = (not np.isnan(devB)) and devB < tol
    okC = maxC < tol and np.abs(closure - 1).max() < 1e-6
    print(f"\nA {'PASS' if okA else 'FAIL'} | "
          f"B {'PASS' if okB else 'FAIL'} | "
          f"C {'PASS' if okC else 'FAIL'}")
    print("VERDICT:", "noM VALIDATED for true zero field"
          if (okA and okB and okC) else "DISCREPANCY - see above")


if __name__ == "__main__":
    main()
