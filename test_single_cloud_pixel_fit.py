"""单云科学计算验证：公式、解析导数、源函数消元、参数恢复和并行保存。"""
import os
os.environ.setdefault("MPLBACKEND", "Agg")
from pathlib import Path
import pickle
import tempfile
import unittest

import numpy as np
from scipy.optimize import minimize_scalar, curve_fit
from chase_background_fit import gaussian_background, estimate_background_center as estimate_gaussian_center
from check_single_cloud_searches import compare_search_budgets
from single_cloud_pixel_fit import (
    C_KM, single_cloud_vel, _residual, _jacobian, _profile_source,
    fit_single_cloud_spectrum, fit_single_cloud_cube, estimate_background_center,
)


class SingleCloudTests(unittest.TestCase):
    def setUp(self):
        self.lam = np.linspace(6560.8, 6564.8, 81)
        self.vel = (self.lam-6562.8)/6562.8*C_KM
        self.bg = 1000-700*np.exp(-(self.vel/30)**2)
        self.p = np.array([np.log(.8), -23, np.log(18), 180.])
        self.obs = single_cloud_vel(self.p, self.vel, self.bg)

    def test_paper_formula_and_jacobian(self):
        tau = .8*np.exp(-((self.lam-6562.8-(-23)*6562.8/C_KM)/(18*6562.8/C_KM))**2)
        expected = self.bg*np.exp(-tau)+180*(1-np.exp(-tau))
        np.testing.assert_allclose(self.obs, expected, atol=2e-8)
        args = self.vel, self.obs, self.bg, 150.
        numeric = np.empty((81,4))
        for k, step in enumerate([1e-5,1e-4,1e-5,1e-3]):
            delta = np.zeros(4)
            delta[k] = step
            numeric[:,k] = (_residual(self.p+delta,*args)-_residual(self.p-delta,*args))/(2*step)
        np.testing.assert_allclose(_jacobian(self.p,*args), numeric, rtol=2e-5, atol=1e-8)

    def test_source_projection_and_bounds(self):
        for offset in [0., -1000., 2000.]:
            obs = self.obs+offset
            score, S = _profile_source(self.p[:3],self.vel,obs,self.bg,100.,500.)
            def objective(s):
                return np.sum(_residual(np.r_[self.p[:3],s],self.vel,obs,self.bg,100.)**2)
            numeric = minimize_scalar(objective, bounds=(0,500), method="bounded")
            self.assertLessEqual(score, numeric.fun+1e-8)
            self.assertTrue(0 <= S <= 500)
        population = np.column_stack((self.p[:3], self.p[:3]+[.1,10,.2]))
        scores, sources = _profile_source(population,self.vel,self.obs,self.bg,100,500)
        for k in range(2):
            np.testing.assert_allclose([scores[k],sources[k]],
                _profile_source(population[:,k],self.vel,self.obs,self.bg,100,500), atol=1e-12)

    def test_known_parameter_recovery_and_invalid_samples(self):
        fit = fit_single_cloud_spectrum(self.obs,self.lam,self.bg)
        np.testing.assert_allclose(fit['params'], [.8,-23,18,180], atol=.005)
        self.assertLess(fit['rmse'],1e-5)
        self.assertAlmostEqual(fit['rmse']**2*fit['n_valid'],fit['raw_sse'])
        obs = np.ma.array(self.obs,mask=False)
        obs.mask[0] = True
        obs[1] = np.nan
        fit = fit_single_cloud_spectrum(obs,self.lam,self.bg)
        self.assertEqual(fit['n_valid'],79)
        with self.assertRaises(ValueError):
            fit_single_cloud_spectrum(np.full(81,np.nan),self.lam,self.bg)

    def test_background_center_and_reference_change(self):
        expected = 6562.873
        bg = 300+100*(self.lam-expected)**2
        est = estimate_background_center(self.lam,bg)
        self.assertAlmostEqual(est['lam0'],expected,places=7)
        k=np.argmin(bg)
        np.testing.assert_array_equal(est['core_indices'],np.arange(k-2,k+3))
        self.assertEqual(len(est['core_wavelength']),5)
        reversed_est = estimate_background_center(self.lam[::-1],bg[::-1])
        self.assertAlmostEqual(reversed_est['lam0'],expected,places=7)
        new_grid = (self.lam-expected)/expected*C_KM
        obs = single_cloud_vel(self.p,new_grid,bg)
        fit = fit_single_cloud_spectrum(obs,self.lam,bg)
        self.assertAlmostEqual(fit['lam0'],expected,places=7)
        np.testing.assert_allclose(fit['params'],[.8,-23,18,180],atol=.005)
        # 同一物理光谱改用新参考波长：速度和速度单位线宽须一致换算。
        rebased = self.p.copy()
        rebased[1] = self.p[1]*6562.8/expected+C_KM*(6562.8-expected)/expected
        rebased[2] += np.log(6562.8/expected)
        np.testing.assert_allclose(single_cloud_vel(rebased,new_grid,self.bg),self.obs,atol=1e-9)
        for bad in [np.ones(81),np.arange(81),np.full(81,np.nan)]:
            with self.assertRaises(ValueError):
                estimate_background_center(self.lam,bad)

    def test_background_matches_workshop_model_and_window(self):
        bg = self.bg + .3*np.sin(np.arange(81))
        sl = slice(15,65)
        fit = estimate_gaussian_center(self.lam,bg,fit_slice=sl)
        p,cov = curve_fit(gaussian_background,self.lam[sl],bg[sl],
                          p0=[-700,6562.82,.5,1000],method='lm')
        self.assertEqual(fit['n_fit_points'],50)
        np.testing.assert_allclose(fit['fit_wavelength'],self.lam[sl])
        self.assertAlmostEqual(fit['lam0'],p[1],places=6)
        np.testing.assert_allclose(fit['fit_model'],gaussian_background(self.lam[sl],*p),atol=.001)

    def test_parallel_equivalence_and_files(self):
        cube = np.stack((self.obs,self.obs+.1*np.sin(np.arange(81)),self.obs,np.full(81,np.nan)),axis=1).reshape(81,2,2)
        with tempfile.TemporaryDirectory(prefix='single_cloud_validation_') as temp:
            serial = fit_single_cloud_cube(cube,self.lam,self.bg,output_dir=temp,n_jobs=1,
                                           n_global_runs=2,save_plots=False,verbose=False)
            parallel = fit_single_cloud_cube(cube,self.lam,self.bg,output_dir=temp,n_jobs=2,
                                             n_global_runs=2,save_plots=True,plot_dpi=80,verbose=False)
            np.testing.assert_allclose(serial['params'],parallel['params'],atol=1e-7,equal_nan=True)
            np.testing.assert_array_equal(serial['success'],parallel['success'])
            self.assertFalse(parallel['success'][1,1])
            self.assertTrue(np.isnan(parallel['params'][1,1]).all())
            out = parallel['output_dir']
            self.assertEqual(len(list((out/'pixel_fits').glob('*.png'))),4)
            with np.load(out/'pixel_fits/fit_y0789_x1515.npz',allow_pickle=False) as f:
                self.assertEqual(f['candidate_x'].shape,(2,4))
                self.assertEqual(f['raw_sse'],f['candidate_raw_sse'].min())
            with np.load(out/'fit_maps.npz',allow_pickle=False) as f:
                np.testing.assert_allclose(f['v'],parallel['params'][...,1],equal_nan=True)
                np.testing.assert_array_equal(f['x_pixels'],[1515,1516])
            with (out/'result.pkl').open('rb') as f:
                loaded = pickle.load(f)
                np.testing.assert_allclose(loaded['params'],parallel['params'],equal_nan=True)
            self.assertTrue((out/'ha_single_cloud_maps.pdf').is_file())
            self.assertTrue((out/'background_center.png').is_file())
            self.assertEqual(parallel['lam0_method'],'background_core_quadratic')
            report = compare_search_budgets(out,budgets=(1,2))
            self.assertEqual(report['comparisons'][0]['n_pixels'],3)
            with np.load(out/'search_stability.npz') as f:
                np.testing.assert_allclose(f['best_raw_sse'][-1],parallel['raw_sse'],equal_nan=True)
                good=np.isfinite(f['best_raw_sse']).all(axis=0)
                self.assertTrue(np.all(f['best_raw_sse'][0][good]>=f['best_raw_sse'][1][good]))


if __name__ == '__main__':
    unittest.main(verbosity=2)
