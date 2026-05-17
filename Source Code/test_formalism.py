"""Plain-assert tests for formalism.convert_formalism (no pytest dep)."""
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
from formalism import convert_formalism

_passed = 0
def check(name, fn):
    global _passed
    fn()
    _passed += 1
    print(f"  pass: {name}")

def expect_valueerror(fn):
    try:
        fn()
    except ValueError:
        return
    raise AssertionError("expected ValueError, none raised")

def test_r2_to_n2():
    # 'R2' set converts R²→N² (inverse B-row +2Λ²D); tag flips, Λ kept.
    src = {'formalism': 'R2', 'Lambda': 1, 'Be': 99.0, 'D': 0.5}
    out = convert_formalism(src)
    assert out['Be'] == 99.0 + 2 * 1 * 0.5, out
    assert out['D'] == 0.5, out
    assert out['formalism'] == 'N2' and out['Lambda'] == 1, out

def test_absent_formalism_is_identity():
    # No declared convention ⇒ returned unchanged (legacy R² back-compat).
    out = convert_formalism({'Be': 7.0})
    assert out == {'Be': 7.0}, out
    assert 'formalism' not in out and 'Lambda' not in out

def test_unknown_formalism_raises():
    expect_valueerror(lambda: convert_formalism({'formalism': 'X2'}))

def test_n2_without_lambda_raises():
    expect_valueerror(lambda: convert_formalism(
        {'formalism': 'N2', 'Be': 1.0}))

def test_caller_dict_not_mutated():
    src = {'formalism': 'R2', 'Lambda': 1, 'Be': 1.0}
    convert_formalism(src)
    assert src == {'formalism': 'R2', 'Lambda': 1, 'Be': 1.0}, src

import warnings as _w

def test_sigma_is_identity():
    src = {'formalism': 'N2', 'Lambda': 0, 'Be': 100.0, 'p+2q': -3.0,
           'p2q_D': 0.05}
    out = convert_formalism(src)
    assert out == {'formalism': 'R2', 'Lambda': 0, 'Be': 100.0,
                   'p+2q': -3.0, 'p2q_D': 0.05}, out   # Σ: values identical, tag flipped

def test_be_d_row():
    out = convert_formalism(
        {'formalism': 'N2', 'Lambda': 1, 'Be': 100.0, 'D': 0.5})
    assert out['Be'] == 100.0 - 2 * 1 * 0.5, out
    assert out['D'] == 0.5, out                     # D unchanged

def test_generic_x_row_lambda2():
    out = convert_formalism(
        {'formalism': 'N2', 'Lambda': 2, 'p+2q': -10.0, 'p2q_D': 0.25})
    assert out['p+2q'] == -10.0 + 4 * 0.25, out     # Λ²=4
    assert out['p2q_D'] == 0.25, out                # X_D unchanged

def test_all_partners_and_passthrough():
    src = {'formalism': 'N2', 'Lambda': 1,
           'Gamma_SR': 7.0, 'Gamma_D': 0.1,
           'q_lD': -4.0, 'q_lD_D': 0.2,
           'ASO': 5e6, 'muE': 1.23}             # ASO/muE: no X_D → identity
    out = convert_formalism(src)
    assert out['Gamma_SR'] == 7.0 + 1 * 0.1, out
    assert out['q_lD'] == -4.0 + 1 * 0.2, out
    assert out['ASO'] == 5e6 and out['muE'] == 1.23, out

def test_sextic_key_raises():
    expect_valueerror(lambda: convert_formalism(
        {'formalism': 'N2', 'Lambda': 1, 'Be': 1.0, 'p2q_H': 1e-9}))

def test_orphan_xd_warns():
    with _w.catch_warnings(record=True) as rec:
        _w.simplefilter('always')
        out = convert_formalism(
            {'formalism': 'N2', 'Lambda': 1, 'p2q_D': 0.05})
    assert any('p2q_D' in str(r.message) for r in rec), rec
    assert out['p2q_D'] == 0.05, out

def test_g_row_lambda1():
    c = 29979.2458
    out = convert_formalism(
        {'formalism': 'N2', 'Lambda': 1,
         'Be': 5743.96, 'D': 1.4e-7 * c, 'Origin': 13284.427})
    exp = 13284.427 + 1 * 5743.96 / c - 1 * (1.4e-7 * c) / c
    assert abs(out['Origin'] - exp) < 1e-12, (out['Origin'], exp)

def test_g_row_lambda2_powers():
    c = 1000.0
    out = convert_formalism(
        {'formalism': 'N2', 'Lambda': 2, 'Be': 10.0, 'D': 2.0,
         'Origin': 100.0}, c_cm=c)
    exp = 100.0 + 4 * 10.0 / c - 16 * 2.0 / c          # Λ²=4, Λ⁴=16
    assert abs(out['Origin'] - exp) < 1e-12, (out['Origin'], exp)

def test_g_row_sigma_unchanged():
    out = convert_formalism(
        {'formalism': 'N2', 'Lambda': 0, 'Be': 50.0, 'Origin': 9.0})
    assert out['Origin'] == 9.0, out

def test_g_row_no_be_no_shift():
    out = convert_formalism(
        {'formalism': 'N2', 'Lambda': 1, 'Origin': 9.0})
    assert out['Origin'] == 9.0, out

def test_g_row_custom_c():
    out = convert_formalism(
        {'formalism': 'N2', 'Lambda': 1, 'Be': 300.0, 'Origin': 0.0},
        c_cm=100.0)
    assert abs(out['Origin'] - 3.0) < 1e-12, out      # 300/100, no D

def test_backcompat_all_entries_unchanged():
    import importlib
    mp = importlib.import_module('molecule_parameters')
    n = 0
    for mol, spins in mp.molecules.items():
        for spin, states in spins.items():
            for st, d in states.items():
                if not isinstance(d, dict):
                    continue
                if (mol, spin, st) in {('RaF','boson','A0'),
                                       ('RaF','fermion','A0')}:
                    continue   # migrated to N2 — covered by oracle test
                before = dict(d)
                after = convert_formalism(dict(d), c_cm=mp.c)
                # all current entries are default-R2 ⇒ identical
                assert after == before, (mol, spin, st,
                    set(before) ^ set(after))
                n += 1
    assert n > 0
    print(f"    swept {n} entries")

def test_raf_a0_regression_oracle():
    import importlib
    mp = importlib.import_module('molecule_parameters')
    importlib.reload(mp)
    c = mp.c
    # boson
    b = mp.get_molecule_params('RaF', 'A', '0', 'boson')
    assert abs(b['Be'] - (0.191015*c - 2*1.4e-7*c)) < 1e-9, b['Be']
    assert abs(b['p+2q'] - (-0.41071*c + 1.9e-7*c)) < 1e-9, b['p+2q']
    old_buggy_b = 13284.427 + 5755.56/c
    corrected_b = 13284.427 + (0.191015*c)/c - 1.4e-7
    assert abs(b['Origin'] - corrected_b) < 1e-9, b['Origin']
    # shift vs. buggy = Λ²ΔBe/c − Λ⁴D/c; D-term = (1.4e-7*c)/c = 1.4e-7
    assert abs((b['Origin'] - old_buggy_b) - ((0.191015*c-5755.56)/c - 1.4e-7)) < 1e-9
    assert b['formalism'] == 'R2' and b['Lambda'] == 1, b   # tag follows: N2→R2, Λ kept
    # fermion
    f = mp.get_molecule_params('RaF', 'A', '0', 'fermion')
    assert abs(f['Be'] - (5729.03 - 2*1.4e-7*c)) < 1e-9, f['Be']
    assert abs(f['p+2q'] - (-0.4109*c + 1.9e-7*c)) < 1e-9, f['p+2q']
    old_buggy_f = 13284.427 + 5755.56/c
    corrected_f = 13284.427 + 5729.03/c - 1.4e-7
    assert abs(f['Origin'] - corrected_f) < 1e-9, f['Origin']
    assert abs((f['Origin'] - old_buggy_f) - ((5729.03-5755.56)/c - 1.4e-7)) < 1e-9

def test_r2_without_lambda_raises():
    expect_valueerror(lambda: convert_formalism({'formalism': 'R2', 'Be': 1.0}))

def test_n2_to_r2_flips_tag_keeps_lambda():
    out = convert_formalism(
        {'formalism': 'N2', 'Lambda': 1, 'Be': 100.0, 'D': 0.5})
    assert out['formalism'] == 'R2' and out['Lambda'] == 1, out
    assert out['Be'] == 100.0 - 2 * 1 * 0.5, out          # N²→R² numerics unchanged

def test_roundtrip_involution():
    # N²→R²→N² is the identity (Table 7.2 quartic is exactly invertible).
    c = 29979.2458
    src = {'formalism': 'N2', 'Lambda': 1, 'Be': 5743.96, 'D': 1.4e-7 * c,
           'p+2q': -12312.77, 'p2q_D': 0.0057, 'Origin': 13284.427, 'ASO': 4.0e7}
    r2 = convert_formalism(src)
    assert r2['formalism'] == 'R2', r2
    back = convert_formalism(r2)
    assert back['formalism'] == 'N2' and back['Lambda'] == 1, back
    for k in ('Be', 'D', 'p+2q', 'p2q_D', 'Origin', 'ASO'):
        assert abs(back[k] - src[k]) < 1e-9, (k, back[k], src[k])

if __name__ == '__main__':
    check('r2_to_n2', test_r2_to_n2)
    check('absent_formalism_is_identity', test_absent_formalism_is_identity)
    check('r2_without_lambda_raises', test_r2_without_lambda_raises)
    check('n2_to_r2_flips_tag', test_n2_to_r2_flips_tag_keeps_lambda)
    check('roundtrip_involution', test_roundtrip_involution)
    check('unknown_formalism_raises', test_unknown_formalism_raises)
    check('n2_without_lambda_raises', test_n2_without_lambda_raises)
    check('caller_dict_not_mutated', test_caller_dict_not_mutated)
    check('sigma_is_identity', test_sigma_is_identity)
    check('be_d_row', test_be_d_row)
    check('generic_x_row_lambda2', test_generic_x_row_lambda2)
    check('all_partners_and_passthrough', test_all_partners_and_passthrough)
    check('sextic_key_raises', test_sextic_key_raises)
    check('orphan_xd_warns', test_orphan_xd_warns)
    check('g_row_lambda1', test_g_row_lambda1)
    check('g_row_lambda2_powers', test_g_row_lambda2_powers)
    check('g_row_sigma_unchanged', test_g_row_sigma_unchanged)
    check('g_row_no_be_no_shift', test_g_row_no_be_no_shift)
    check('g_row_custom_c', test_g_row_custom_c)
    check('backcompat_all_entries_unchanged', test_backcompat_all_entries_unchanged)
    check('raf_a0_regression_oracle', test_raf_a0_regression_oracle)
    print(f"OK: {_passed} passed")
