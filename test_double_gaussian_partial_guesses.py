import unittest

import numpy as np

from gaussian_fit import fit_double_gaussian


class TestDoubleGaussianPartialGuesses(unittest.TestCase):
    def test_fit_double_gaussian_accepts_individual_initial_guesses(self):
        x = np.linspace(1402.2, 1403.3, 260)
        y = (
            5.0
            + 18.0 * np.exp(-0.5 * ((x - 1402.68) / 0.055) ** 2)
            + 8.0 * np.exp(-0.5 * ((x - 1402.91) / 0.09) ** 2)
        )

        result = fit_double_gaussian(
            x,
            y,
            center1_guess=1402.68,
            sigma1_guess=0.05,
            center2_guess=1402.91,
            sigma2_guess=0.08,
            positive=True,
        )

        self.assertTrue(result["success"])
        self.assertAlmostEqual(result["center1"], 1402.68, delta=0.01)
        self.assertAlmostEqual(result["center2"], 1402.91, delta=0.01)
        self.assertAlmostEqual(result["sigma1"], 0.055, delta=0.01)
        self.assertAlmostEqual(result["sigma2"], 0.09, delta=0.01)
        self.assertAlmostEqual(result["background"], 5.0, delta=0.1)

    def test_fit_double_gaussian_accepts_linear_background_guess(self):
        x = np.linspace(1393.0, 1395.0, 260)
        x_ref = np.mean(x)
        y = (
            3.0
            + 0.6 * (x - x_ref)
            + 12.0 * np.exp(-0.5 * ((x - 1393.7) / 0.12) ** 2)
            + 5.0 * np.exp(-0.5 * ((x - 1394.25) / 0.18) ** 2)
        )

        result = fit_double_gaussian(
            x,
            y,
            center1_guess=1393.7,
            center2_guess=1394.25,
            slope_guess=0.4,
            linear_background=True,
            positive=True,
        )

        self.assertTrue(result["success"])
        self.assertAlmostEqual(result["center1"], 1393.7, delta=0.02)
        self.assertAlmostEqual(result["center2"], 1394.25, delta=0.02)
        self.assertAlmostEqual(result["slope"], 0.6, delta=0.05)
        self.assertEqual(len(result["params"]), 8)


if __name__ == "__main__":
    unittest.main()
