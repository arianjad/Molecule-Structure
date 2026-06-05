"""Compute 225RaF A(0)-X(0) hyperfine level energies + branching for the level diagram.

Reproduces RaF225 X-A.ipynb (cells 26-31) with DATABASE-DEFAULT A-state params
(canonical; the cell-3 dYb-flip / h1/2Yb=0 was an experiment per Arian, 2026-06-05).
Also runs the experimental param set so we can quantify the difference.

Run from this directory in the Structure env:
    conda run -n Structure python compute_225_xa_branching.py

Outputs: prints a verification report AND writes data_225_level_diagram.py
(the frozen numbers the renderer imports).

Provenance: branching_ratios(X0,A0,0,1e-8); BR[g,e]=(G_evecs_aBJ @ TDM @ E_evecs.T)^2,
shape [ground x excited], M='none' (exact closed-form, thesis Eq 3.3; columns close to 1
over v''=0 manifolds). See reports/coherent_vs_incoherent_branching.md.
"""
import os
import numpy as np
from config_path import add_to_sys_path
add_to_sys_path()
from Energy_Levels import MoleculeLevels, branching_ratios

np.set_printoptions(precision=8, suppress=True)
HERE = os.path.dirname(os.path.abspath(__file__))


def make_X():
    return MoleculeLevels.initialize_state(
        molecule_name='RaF', elec_state='X', vib_state=0,
        N_list=np.arange(0, 4, 1), fermion_or_boson='fermion',
        M_sublevels='none', M_list=[1/2], I_nuclei=[1/2, 1/2],
        isotope=225, round=8, params=None, P_values=[1/2])


def make_A():
    return MoleculeLevels.initialize_state(
        molecule_name='RaF', elec_state='A', vib_state=0,
        N_list=np.arange(1, 4, 1), fermion_or_boson='fermion',
        M_sublevels='none', M_list=[1/2], I_nuclei=[1/2, 1/2],
        isotope=225, round=8, params=None, P_values=[1/2, 3/2])


def qn(state, i, qs):
    """Dominant-basis quantum numbers of eigenstate i (same rule select_q uses)."""
    dom = int(np.argmax(state.evecs0[i] ** 2))
    out = {}
    for q in qs:
        v = state.q_numbers[q][dom]
        out[q] = float(v)
    return out


def par(state, i):
    return '+' if state.parities[i] > 0 else '-'


def ground_table(X):
    qs = ['N', 'G', 'F1', 'F']
    rows = []
    for i in range(X.size):
        d = qn(X, i, qs)
        rows.append(dict(idx=i, parity=par(X, i), E=float(X.evals0[i]), **d))
    return rows


def excited_table(A):
    qs = ['J', 'P', 'F1', 'F']
    rows = []
    for i in range(A.size):
        d = qn(A, i, qs)
        rows.append(dict(idx=i, parity=par(A, i), E=float(A.evals0[i]), **d))
    return rows


def find_cycling_upper(A, etab):
    """A J=1/2, P=1/2, parity +, F1=0 (the single F=1/2 state)."""
    cands = [r for r in etab
             if abs(r['J'] - 0.5) < 1e-6 and abs(r['P'] - 0.5) < 1e-6
             and r['parity'] == '+' and abs(r['F1'] - 0.0) < 1e-6]
    assert len(cands) == 1, f"expected 1 cycling upper state, got {len(cands)}: {cands}"
    return cands[0]['idx']


def grouped_branching(BR, gtab, e_idx):
    """Sum BR column e_idx over F within each ground (N,G,F1). Returns dict + raw column."""
    col = BR[:, e_idx]
    groups = {}
    for r in gtab:
        key = (int(r['N']), int(r['G']), int(r['F1']))
        groups.setdefault(key, 0.0)
        groups[key] += float(col[r['idx']])
    return groups, col


def run(param_set):
    """param_set: 'default' or 'experimental'."""
    X = make_X()
    A = make_A()
    if param_set == 'experimental':
        A.update_params({'dYb': -A.parameters['dYb'], 'h1/2Yb': 0}, recompute=True)
    elif param_set == 'physical':
        # flip d to negative (225Ra mu<0), KEEP real Skripnikov h1/2Yb=-1426
        A.update_params({'dYb': -A.parameters['dYb']}, recompute=True)
    BR = branching_ratios(X, A, 0, 1e-8)   # also sets evals0/evecs0/parities
    gtab = ground_table(X)
    etab = excited_table(A)
    e_up = find_cycling_upper(A, etab)
    groups, col = grouped_branching(BR, gtab, e_up)
    return dict(X=X, A=A, BR=BR, gtab=gtab, etab=etab, e_up=e_up,
                groups=groups, col=col, param_set=param_set)


def a_j12_levels(etab):
    """A J=1/2 parity+ F1 levels, energy relative to the lowest of them (MHz)."""
    sel = [r for r in etab if abs(r['J'] - 0.5) < 1e-6 and abs(r['P'] - 0.5) < 1e-6
           and r['parity'] == '+']
    e0 = min(r['E'] for r in sel)
    return [dict(F1=int(r['F1']), F=r['F'], E_rel=r['E'] - e0, idx=r['idx']) for r in sel]


def x_n1_levels(gtab):
    """X N=1 levels (both G), energy relative to the lowest N=1 level (MHz)."""
    sel = [r for r in gtab if int(r['N']) == 1]
    e0 = min(r['E'] for r in sel)
    return [dict(G=int(r['G']), F1=int(r['F1']), F=r['F'], E_rel=r['E'] - e0, idx=r['idx'])
            for r in sel]


def report(res):
    ps = res['param_set']
    print(f"\n{'='*70}\n  PARAM SET: {ps}\n{'='*70}")
    print(f"A dYb = {res['A'].parameters['dYb']}, h1/2Yb = {res['A'].parameters.get('h1/2Yb')}")
    print(f"X size={res['X'].size}, A size={res['A'].size}, BR shape={res['BR'].shape}")

    print("\n-- A J=1/2 parity+ levels (rel to lowest, MHz; /1000 = GHz) --")
    for r in sorted(a_j12_levels(res['etab']), key=lambda d: d['E_rel']):
        print(f"  F1={r['F1']} F={r['F']:.1f}   {r['E_rel']:12.3f} MHz  = {r['E_rel']/1e3:8.4f} GHz")

    print("\n-- X N=1 levels (rel to lowest N=1, MHz) --")
    for r in sorted(x_n1_levels(res['gtab']), key=lambda d: d['E_rel']):
        print(f"  G={r['G']} F1={r['F1']} F={r['F']:.1f}   {r['E_rel']:10.3f} MHz")

    e_up = res['e_up']
    print(f"\n-- cycling upper state idx={e_up} (A J=1/2 P=1/2 + F1=0) --")
    print(f"   BR column sum over ALL ground = {res['col'].sum():.8f}  (should be ~1, M='none' closure)")

    print("\n-- grouped branching from cycling upper state, by ground (N,G,F1) --")
    for key in sorted(res['groups'], key=lambda k: -res['groups'][k]):
        N, G, F1 = key
        f = res['groups'][key]
        tag = ''
        if N == 1 and G == 1 and F1 == 1:
            tag = '  <- main return (closed loop)'
        elif N == 1 and G == 1 and F1 in (0, 2):
            tag = '  <- expect ~0 (selection rule: F1=0 ->/-> F1=' + str(F1) + ')'
        elif N == 1 and G == 0:
            tag = '  <- G=0 leak'
        elif N == 3:
            tag = '  <- N=3 leak'
        print(f"  (N={N}, G={G}, F1={F1}):  {f:.3e}{tag}")
    return res


def main():
    res_def = report(run('default'))
    res_exp = report(run('experimental'))
    res_phys = report(run('physical'))

    print(f"\n{'='*70}\n  DEFAULT vs EXPERIMENTAL comparison\n{'='*70}")
    print("A J=1/2+ F1 splittings (GHz):")
    for label, res in [('default', res_def), ('experimental', res_exp), ('physical', res_phys)]:
        lv = {r['F1']: r['E_rel']/1e3 for r in a_j12_levels(res['etab'])
              if abs(r['F'] - 0.5) < 1e-6 or True}
        # report per (F1,F)
        s = ', '.join(f"F1={r['F1']}/F={r['F']:.0f}:{r['E_rel']/1e3:.4f}"
                      for r in sorted(a_j12_levels(res['etab']), key=lambda d: d['E_rel']))
        print(f"  {label:12s}: {s}")
    print("\nLeak fractions (grouped):")
    for key in [(1, 1, 1), (1, 1, 2), (1, 1, 0), (1, 0, 1), (3, 1, 2)]:
        d = res_def['groups'].get(key, 0.0)
        e = res_exp['groups'].get(key, 0.0)
        p = res_phys['groups'].get(key, 0.0)
        print(f"  {str(key):14s} default={d:.4e}  exp={e:.4e}  phys={p:.4e}")

    # ---- write the frozen data module (PHYSICAL = canonical per Arian 2026-06-05) ----
    # negative d (225Ra mu<0) + real Skripnikov h1/2Yb=-1426. Branching is identical
    # to default/experimental (insensitive to A hyperfine); only A energies differ.
    res = res_phys
    a_levels = sorted(a_j12_levels(res['etab']), key=lambda d: d['E_rel'])
    x_levels = sorted(x_n1_levels(res['gtab']), key=lambda d: d['E_rel'])
    groups = res['groups']

    def pyfloat(x):
        return repr(round(float(x), 8))

    lines = []
    lines.append('"""Frozen 225RaF A(0)-X(0) hyperfine numbers for the level diagram.')
    lines.append('')
    lines.append('AUTO-GENERATED by compute_225_xa_branching.py (PHYSICAL/canonical params).')
    lines.append('Do not hand-edit; re-run the compute script.')
    lines.append('')
    lines.append('PROVENANCE: RaF225 X-A.ipynb cells 26-31, A-state with dYb -> -dYb = -0.076*c')
    lines.append('(225Ra mu_I<0, opposite 171Yb; the +0.076*c default is a borrowed "Wilkins"')
    lines.append('171Yb placeholder) and the REAL Skripnikov h1/2Yb=-1426 (NOT the scratch h1/2=0).')
    lines.append('Canonical per Arian 2026-06-05. branching_ratios(X0,A0,0,1e-8), M=none.')
    lines.append('NB: branching is INSENSITIVE to the A hyperfine params (default==exp==phys to')
    lines.append('4 sig figs); only the A-state energies depend on this choice.')
    lines.append('All energies MHz, relative to the lowest level of their manifold group.')
    lines.append('Branching = conditional on decay to v\'\'=0 (hyperfine fraction only).')
    lines.append('"""')
    lines.append('')
    lines.append(f'CYCLING_UPPER_BR_SUM = {pyfloat(res["col"].sum())}  # column closure check')
    lines.append('')
    lines.append('# A 2Pi_1/2 (v=0) J=1/2 parity+ levels: list of (F1, F, E_rel_MHz)')
    lines.append('A_LEVELS = [')
    for r in a_levels:
        lines.append(f'    ({r["F1"]}, {pyfloat(r["F"])}, {pyfloat(r["E_rel"])}),')
    lines.append(']')
    lines.append('')
    lines.append('# X 2Sigma (v=0) N=1 levels: list of (G, F1, F, E_rel_MHz)')
    lines.append('X_N1_LEVELS = [')
    for r in x_levels:
        lines.append(f'    ({r["G"]}, {r["F1"]}, {pyfloat(r["F"])}, {pyfloat(r["E_rel"])}),')
    lines.append(']')
    lines.append('')
    lines.append('# Branching from cycling upper state (A J=1/2+ F1=0), grouped by ground')
    lines.append('# (N, G, F1), summed over F. Conditional on v\'\'=0.')
    lines.append('BRANCHING_BY_F1 = {')
    for key in sorted(groups, key=lambda k: -groups[k]):
        lines.append(f'    {key}: {repr(float(groups[key]))},')  # full precision (keep tiny leaks)
    lines.append('}')
    lines.append('')

    out_path = os.path.join(HERE, 'data_225_level_diagram.py')
    with open(out_path, 'w') as f:
        f.write('\n'.join(lines) + '\n')
    print(f"\nWROTE {out_path}")


if __name__ == '__main__':
    main()
