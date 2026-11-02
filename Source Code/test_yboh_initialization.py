"""Regression checks for odd-Yb isotope selection through the public API."""
import unittest

import numpy as np

from Energy_Levels import MoleculeLevels
from molecule_parameters import get_molecule_params


class YbOHInitializationTests(unittest.TestCase):
    def make_state(self, isotope, vib_state, params=None):
        return MoleculeLevels.initialize_state(
            'YbOH', 'X', vib_state, N_list=[1, 1],
            fermion_or_boson='fermion', isotope=isotope, params=params,
            I_nuclei=[2.5 if str(isotope) == '173' else 0.5, 0.5],
            M_sublevels='custom', M_list=[0.5], P_values=[0.5],
        )

    def test_default_and_explicit_isotopes_select_matching_parameters(self):
        for isotope in (None, 171, '171', 173, '173'):
            label = '171' if isotope is None else str(isotope)
            for vib_state in ('000', '010'):
                with self.subTest(isotope=isotope, vib_state=vib_state):
                    state = self.make_state(isotope, vib_state)
                    expected = get_molecule_params('YbOH', label + 'X', vib_state, 'fermion')
                    self.assertEqual(state.iso_state, label + 'X' + vib_state)
                    self.assertEqual(state.parameters, expected)
                    H = state.H_function(0, 1)
                    self.assertTrue(np.isfinite(H).all())
                    np.testing.assert_allclose(H, H.T.conj(), atol=1e-10, rtol=0)

    def test_explicit_parameter_path_agrees_and_preserves_input(self):
        for isotope in (171, 173):
            with self.subTest(isotope=isotope):
                params = get_molecule_params('YbOH', str(isotope) + 'X', '010', 'fermion')
                original = dict(params)
                explicit = self.make_state(isotope, '010', params=params)
                automatic = self.make_state(isotope, '010')
                self.assertEqual(params, original)
                np.testing.assert_array_equal(
                    automatic.H_function(0, 1), explicit.H_function(0, 1))


if __name__ == '__main__':
    unittest.main()
