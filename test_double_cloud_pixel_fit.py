"""双云加速的数值正确性、并行坐标和保存格式回归检查。"""

import csv
import inspect
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

import double_cloud_pixel_fit as dc


class DoubleCloudFitTests(unittest.TestCase):
    def setUp(self):
        self.lam0 = 6562.87
        self.wave = np.linspace(self.lam0 - 1.8, self.lam0 + 1.8, 78)
        self.vel = (self.wave - self.lam0) / self.lam0 * 299792.458
        self.bg = 1050 - 750 * np.exp(-((self.wave - self.lam0) / 0.55)**2)
        self.p = np.array([np.log(0.8), np.log(0.6), -22, 26, np.log(20), np.log(23), 70, 130])
        self.obs = dc.double_cloud_vel(self.p, self.vel, self.bg)

    def test_vectorized_matches_independent_scalar_objectives(self):
        rng = np.random.default_rng(309)
        population = self.p[:, None] + rng.normal(size=(8, 17)) * np.array([1, 1, 20, 20, 0.6, 0.6, 50, 50])[:, None]
        args = (self.vel, self.obs, self.bg, 22.42)
        expected = np.array([dc._objective(p, *args) for p in population.T])
        np.testing.assert_allclose(dc._objective(population, *args), expected, rtol=2e-14, atol=2e-12)
        self.assertEqual(dc._objective(population[:, :0], *args).shape, (0,))

    def test_analytic_jacobian_matches_central_differences(self):
        for shifts in ([0] * 8, [-5, 3, 4, -8, 0.5, -0.4, 100, -80], [3, -5, -4, 8, -0.3, 0.4, -50, 200]):
            p = self.p + shifts
            numerical = np.empty((self.wave.size, 8))
            args = (self.vel, self.obs, self.bg, 22.42)
            for k in range(8):
                h = 1e-5 * max(abs(p[k]), 1)
                delta = np.zeros(8)
                delta[k] = h
                numerical[:, k] = (dc._residual(p + delta, *args) - dc._residual(p - delta, *args)) / (2 * h)
            np.testing.assert_allclose(dc._jacobian(p, *args), numerical, rtol=2e-5, atol=2e-7)

    def test_noiseless_spectrum_recovery_with_analytic_refinement(self):
        start = self.p + np.array([0.1, -0.1, 2, -2, 0.05, -0.05, 5, -5])
        args = (self.vel, self.obs, self.bg, 22.42)
        fit = least_squares(dc._residual, start, jac=dc._jacobian, args=args, x_scale="jac", max_nfev=2000, ftol=1e-12, xtol=1e-12, gtol=1e-12)
        self.assertTrue(fit.success)
        np.testing.assert_allclose(dc.double_cloud_vel(fit.x, self.vel, self.bg), self.obs, atol=1e-6)

    def test_serial_parallel_alignment_failure_and_saved_candidates(self):
        cube = np.stack([self.obs, self.obs * 0.97, self.obs * 1.02, np.full_like(self.obs, np.nan)], axis=1).reshape(78, 2, 2)
        with tempfile.TemporaryDirectory() as temp:
            settings = dict(output_dir=temp, y_start=789, x_start=1515, lam0=self.lam0, n_global_runs=2, de_maxiter=8, de_popsize=5, max_nfev=100, save_plots=False, verbose=False)
            serial = dc.fit_double_cloud_cube(cube, self.wave, self.bg, n_jobs=1, **settings)
            parallel = dc.fit_double_cloud_cube(cube, self.wave, self.bg, n_jobs=2, **settings)
            np.testing.assert_allclose(parallel["model"], serial["model"], rtol=1e-7, atol=1e-5, equal_nan=True)
            np.testing.assert_allclose(parallel["raw_sse"], serial["raw_sse"], rtol=1e-6, atol=1e-7, equal_nan=True)
            self.assertTrue(np.isnan(parallel["params"][1, 1]).all())
            self.assertFalse(parallel["success"][1, 1])
            out = parallel["output_dir"]
            self.assertEqual(len(list((out / "pixel_fits").glob("*.npz"))), 4)
            self.assertFalse(list(out.rglob("*.png")))
            with open(out / "fit_parameters.csv", encoding="utf-8-sig", newline="") as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual({(int(row["y"]), int(row["x"])) for row in rows}, {(789, 1515), (789, 1516), (790, 1515), (790, 1516)})
            with np.load(out / "pixel_fits" / "fit_y0789_x1516.npz") as saved:
                self.assertEqual(saved["candidate_x"].shape, (2, 8))
                np.testing.assert_array_equal(saved["intensity"], cube[:, 0, 1])
                self.assertEqual(saved["candidate_seed"].min(), 20260708)
                self.assertAlmostEqual(float(saved["raw_sse"]), float(np.sum((saved["model"] - saved["intensity"])**2)), places=8)
            metadata = json.loads((out / "settings.json").read_text(encoding="utf-8"))
            self.assertEqual(metadata["effective_workers"], 2)
            self.assertTrue(metadata["vectorized"])

    def test_scalar_search_and_seed_validation(self):
        result = dc.fit_double_cloud_spectrum(self.obs, self.wave, self.bg, lam0=self.lam0, n_global_runs=1, de_maxiter=5, de_popsize=5, max_nfev=40, vectorized=False)
        self.assertFalse(result["vectorized"])
        self.assertEqual(len(result["candidates"]), 1)
        self.assertLessEqual(result["scaled_sse"], result["candidates"][0]["de_fun"] + 1e-10)
        with self.assertRaises(ValueError):
            dc.fit_double_cloud_spectrum(self.obs, self.wave, self.bg, base_seed=2**32 - 1, n_global_runs=2)
        self.assertEqual(inspect.signature(dc.fit_double_cloud_cube).parameters["n_jobs"].default, 16)


if __name__ == "__main__":
    unittest.main()
