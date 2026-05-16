"""Plain-assert tests for formalism.convert_params_to_engine_R2 (no pytest dep)."""
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
from formalism import convert_params_to_engine_R2

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

def test_r2_passthrough_and_strip():
    src = {'formalism': 'R2', 'Be': 100.0, 'p+2q': -3.0}
    out = convert_params_to_engine_R2(src)
    assert out == {'Be': 100.0, 'p+2q': -3.0}, out
    assert 'formalism' not in out and 'Lambda' not in out

def test_absent_formalism_is_r2():
    out = convert_params_to_engine_R2({'Be': 7.0})
    assert out == {'Be': 7.0}, out

def test_unknown_formalism_raises():
    expect_valueerror(lambda: convert_params_to_engine_R2({'formalism': 'X2'}))

def test_n2_without_lambda_raises():
    expect_valueerror(lambda: convert_params_to_engine_R2(
        {'formalism': 'N2', 'Be': 1.0}))

def test_caller_dict_not_mutated():
    src = {'formalism': 'R2', 'Be': 1.0}
    convert_params_to_engine_R2(src)
    assert src == {'formalism': 'R2', 'Be': 1.0}, src

import warnings as _w

def test_sigma_is_identity():
    src = {'formalism': 'N2', 'Lambda': 0, 'Be': 100.0, 'p+2q': -3.0,
           'p2q_D': 0.05}
    out = convert_params_to_engine_R2(src)
    assert out == {'Be': 100.0, 'p+2q': -3.0, 'p2q_D': 0.05}, out

def test_be_d_row():
    out = convert_params_to_engine_R2(
        {'formalism': 'N2', 'Lambda': 1, 'Be': 100.0, 'D': 0.5})
    assert out['Be'] == 100.0 - 2 * 1 * 0.5, out
    assert out['D'] == 0.5, out                     # D unchanged

def test_generic_x_row_lambda2():
    out = convert_params_to_engine_R2(
        {'formalism': 'N2', 'Lambda': 2, 'p+2q': -10.0, 'p2q_D': 0.25})
    assert out['p+2q'] == -10.0 + 4 * 0.25, out     # Λ²=4
    assert out['p2q_D'] == 0.25, out                # X_D unchanged

def test_all_partners_and_passthrough():
    src = {'formalism': 'N2', 'Lambda': 1,
           'Gamma_SR': 7.0, 'Gamma_D': 0.1,
           'q_lD': -4.0, 'q_lD_D': 0.2,
           'ASO': 5e6, 'muE': 1.23}             # ASO/muE: no X_D → identity
    out = convert_params_to_engine_R2(src)
    assert out['Gamma_SR'] == 7.0 + 1 * 0.1, out
    assert out['q_lD'] == -4.0 + 1 * 0.2, out
    assert out['ASO'] == 5e6 and out['muE'] == 1.23, out

def test_sextic_key_raises():
    expect_valueerror(lambda: convert_params_to_engine_R2(
        {'formalism': 'N2', 'Lambda': 1, 'Be': 1.0, 'p2q_H': 1e-9}))

def test_orphan_xd_warns():
    with _w.catch_warnings(record=True) as rec:
        _w.simplefilter('always')
        out = convert_params_to_engine_R2(
            {'formalism': 'N2', 'Lambda': 1, 'p2q_D': 0.05})
    assert any('p2q_D' in str(r.message) for r in rec), rec
    assert out['p2q_D'] == 0.05, out

def test_g_row_lambda1():
    c = 29979.2458
    out = convert_params_to_engine_R2(
        {'formalism': 'N2', 'Lambda': 1,
         'Be': 5743.96, 'D': 1.4e-7 * c, 'Origin': 13284.427})
    exp = 13284.427 + 1 * 5743.96 / c - 1 * (1.4e-7 * c) / c
    assert abs(out['Origin'] - exp) < 1e-12, (out['Origin'], exp)

def test_g_row_lambda2_powers():
    c = 1000.0
    out = convert_params_to_engine_R2(
        {'formalism': 'N2', 'Lambda': 2, 'Be': 10.0, 'D': 2.0,
         'Origin': 100.0}, c_cm=c)
    exp = 100.0 + 4 * 10.0 / c - 16 * 2.0 / c          # Λ²=4, Λ⁴=16
    assert abs(out['Origin'] - exp) < 1e-12, (out['Origin'], exp)

def test_g_row_sigma_unchanged():
    out = convert_params_to_engine_R2(
        {'formalism': 'N2', 'Lambda': 0, 'Be': 50.0, 'Origin': 9.0})
    assert out['Origin'] == 9.0, out

def test_g_row_no_be_no_shift():
    out = convert_params_to_engine_R2(
        {'formalism': 'N2', 'Lambda': 1, 'Origin': 9.0})
    assert out['Origin'] == 9.0, out

def test_g_row_custom_c():
    out = convert_params_to_engine_R2(
        {'formalism': 'N2', 'Lambda': 1, 'Be': 300.0, 'Origin': 0.0},
        c_cm=100.0)
    assert abs(out['Origin'] - 3.0) < 1e-12, out      # 300/100, no D

if __name__ == '__main__':
    check('r2_passthrough_and_strip', test_r2_passthrough_and_strip)
    check('absent_formalism_is_r2', test_absent_formalism_is_r2)
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
    print(f"OK: {_passed} passed")
