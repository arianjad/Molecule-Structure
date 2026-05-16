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

if __name__ == '__main__':
    check('r2_passthrough_and_strip', test_r2_passthrough_and_strip)
    check('absent_formalism_is_r2', test_absent_formalism_is_r2)
    check('unknown_formalism_raises', test_unknown_formalism_raises)
    check('n2_without_lambda_raises', test_n2_without_lambda_raises)
    check('caller_dict_not_mutated', test_caller_dict_not_mutated)
    print(f"OK: {_passed} passed")
