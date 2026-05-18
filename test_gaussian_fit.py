import unittest

import numpy as np

from gaussian_fit import fit_double_gaussian, fit_gaussian


class TestGaussianFit(unittest.TestCase):
    def test_fit_gaussian_recovers_basic_parameters(self):
        rng = np.random.default_rng(123)
        x = np.linspace(-5.0, 5.0, 201)
        expected = {
            "amplitude": 4.0,
            "center": 0.8,
            "sigma": 1.2,
            "background": 0.5,
        }
        y_clean = (
            expected["background"]
            + expected["amplitude"]
            * np.exp(-0.5 * ((x - expected["center"]) / expected["sigma"]) ** 2)
        )
        y = y_clean + rng.normal(0.0, 0.03, size=x.size)

        result = fit_gaussian(x, y)

        self.assertTrue(result["success"])
        self.assertAlmostEqual(result["amplitude"], expected["amplitude"], delta=0.15)
        self.assertAlmostEqual(result["center"], expected["center"], delta=0.05)
        self.assertAlmostEqual(result["sigma"], expected["sigma"], delta=0.08)
        self.assertAlmostEqual(result["background"], expected["background"], delta=0.08)
        self.assertAlmostEqual(
            result["area"],
            expected["amplitude"] * expected["sigma"] * np.sqrt(2.0 * np.pi),
            delta=0.35,
        )
        self.assertEqual(result["yfit"].shape, x.shape)

    def test_fit_gaussian_with_linear_background_recovers_slope(self):
        x = np.linspace(1393.0, 1395.0, 180)
        amplitude = 9.0
        center = 1394.1
        sigma = 0.18
        background = 2.0
        slope = 0.7
        y = (
            background
            + slope * (x - np.mean(x))
            + amplitude * np.exp(-0.5 * ((x - center) / sigma) ** 2)
        )

        result = fit_gaussian(x, y, linear_background=True)

        self.assertTrue(result["success"])
        self.assertAlmostEqual(result["center"], center, delta=0.01)
        self.assertAlmostEqual(result["sigma"], sigma, delta=0.01)
        self.assertAlmostEqual(result["slope"], slope, delta=0.05)
        self.assertEqual(len(result["params"]), 5)

    def test_fit_double_gaussian_recovers_two_components(self):
        rng = np.random.default_rng(456)
        x = np.linspace(-4.0, 4.0, 240)
        amp1, center1, sigma1 = 7.0, -0.7, 0.45
        amp2, center2, sigma2 = 3.5, 0.9, 0.65
        background = 1.2
        y_clean = (
            background
            + amp1 * np.exp(-0.5 * ((x - center1) / sigma1) ** 2)
            + amp2 * np.exp(-0.5 * ((x - center2) / sigma2) ** 2)
        )
        y = y_clean + rng.normal(0.0, 0.02, size=x.size)

        result = fit_double_gaussian(
            x,
            y,
            p0=[amp1, center1, sigma1, amp2, center2, sigma2, background],
        )

        self.assertTrue(result["success"])
        self.assertAlmostEqual(result["center1"], center1, delta=0.05)
        self.assertAlmostEqual(result["center2"], center2, delta=0.05)
        self.assertAlmostEqual(result["sigma1"], sigma1, delta=0.05)
        self.assertAlmostEqual(result["sigma2"], sigma2, delta=0.05)
        self.assertAlmostEqual(result["background"], background, delta=0.08)
        self.assertAlmostEqual(
            result["total_area"],
            (amp1 * sigma1 + amp2 * sigma2) * np.sqrt(2.0 * np.pi),
            delta=0.35,
        )
        self.assertEqual(result["component1"].shape, x.shape)
        self.assertEqual(result["component2"].shape, x.shape)
        self.assertEqual(result["yfit"].shape, x.shape)


if __name__ == "__main__":
    unittest.main()
