"""CHASE H-alpha 单云模型：向量化搜索、解析源函数、多进程逐像素拟合。

I = I0*exp(-tau) + S*(1-exp(-tau));
tau = tau0*exp(-((v_grid-v)/W)**2).
参数顺序 [tau0, v, W, S]，v/W 单位 km/s，S 与输入强度同单位。
不做临边昏暗修正。依赖 numpy、scipy>=1.9、matplotlib、joblib>=1.4。
"""

from datetime import datetime
from pathlib import Path
from time import perf_counter
import csv
import json
import os
import pickle

import numpy as np
from scipy.optimize import differential_evolution, least_squares

C_KM = 299792.458
PARAM_NAMES = ("tau0", "v", "W", "S")
DEFAULT_OUT = Path(r"C:\Learning\PHD1st\magnetic_reconnecion\data\CHASE_Ha_image\single_cloud")


def _array(a):
    return np.asarray(np.ma.asarray(a, dtype=float).filled(np.nan))


from background_core_fit import estimate_background_center, save_background_center as _save_background_center


def single_cloud_vel(params, vel, I0):
    """优化坐标为 [log(tau0), v, log(W), S]，不是物理参数数组。"""
    ltau, v, lW, S = params
    tau = np.exp(ltau - ((vel - v) / np.exp(lW))**2)
    return I0 * np.exp(-tau) + S * (-np.expm1(-tau))


def _residual(p, vel, obs, bg, scale):
    return (single_cloud_vel(p, vel, bg) - obs) / scale


def _jacobian(p, vel, obs, bg, scale):
    """残差对四个优化参数的解析导数，避免数值差分重复求模型。"""
    ltau, v, lW, S = p
    W = np.exp(lW)
    z = (vel - v) / W
    tau = np.exp(ltau - z**2)
    common = (S - bg) * np.exp(-tau) * tau
    return np.column_stack((common, common * 2 * z / W,
                            common * 2 * z**2, -np.expm1(-tau))) / scale


def _profile_source(p, vel, obs, bg, scale, s_max):
    """固定前三个参数，解析消去有界 S；支持 (3,) 或 (3,population)。

    a=1-exp(-tau), b=I0*exp(-tau):
    S_hat=clip(sum(a*(I-b))/sum(a*a), 0, s_max).
    这是同一个有界、等权 raw-SSE 目标，不是近似源函数或固定 S。
    """
    p = np.asarray(p)
    scalar = p.ndim == 1
    q = p[:, None] if scalar else p
    tau = np.exp(q[0] - ((vel[:, None] - q[1]) / np.exp(q[2]))**2)
    a = -np.expm1(-tau)
    b = bg[:, None] * np.exp(-tau)
    denominator = np.sum(a*a, axis=0)
    S = np.divide(np.sum(a*(obs[:, None] - b), axis=0), denominator,
                  out=np.zeros_like(denominator), where=denominator > 0)
    S = np.clip(S, 0, s_max)
    score = np.sum(((b + a*S - obs[:, None]) / scale)**2, axis=0)
    return (float(score[0]), float(S[0])) if scalar else (score, S)


def _profile_objective(p, *args):
    return _profile_source(p, *args)[0]


def _check_settings(lam0, n_global_runs, base_seed, de_maxiter, de_popsize, max_nfev, w_bounds):
    if not np.isfinite(lam0) or lam0 <= 0:
        raise ValueError("lam0 必须是正的有限波长 (Angstrom)")
    for name, value in (("n_global_runs", n_global_runs), ("de_maxiter", de_maxiter),
                        ("de_popsize", de_popsize), ("max_nfev", max_nfev)):
        if not isinstance(value, (int, np.integer)) or value < 1:
            raise ValueError(f"{name} 必须为正整数")
    if not isinstance(base_seed, (int, np.integer)) or not 0 <= base_seed <= 2**32-n_global_runs:
        raise ValueError("base_seed 及后续种子必须在 [0,2**32) 内")
    if len(w_bounds) != 2 or not np.all(np.isfinite(w_bounds)) or not 0 < w_bounds[0] < w_bounds[1]:
        raise ValueError("w_bounds 必须满足 0 < W_MIN < W_MAX，单位 km/s")


def fit_single_cloud_spectrum(
    intensity, wavelength, I0, *, lam0=None, n_global_runs=30,
    base_seed=20260706, de_maxiter=500, de_popsize=15,
    max_nfev=2000, w_bounds=(0.05, 200.0), background_core_points=5,
):
    """单像素单云拟合；返回最佳候选及全部搜索结果。

    输入三个等长一维数组；波长单位 Angstrom。共同有效点数须 >4。
    lam0=None 时仅用背景谷底附近 background_core_points 个点局部二次拟合线心。
    也可以显式传入已标定的背景线心。批处理中统一估计一次，全部像素共用。
    W 是 exp[-((v-v0)/W)^2] 的宽度，不是 Gaussian sigma。
    success 只表示局部优化收敛，不表示唯一解或速度误差。
    """
    start = perf_counter()
    obs_all, lam, bg_all = map(_array, (intensity, wavelength, I0))
    if obs_all.ndim != 1 or obs_all.shape != lam.shape or bg_all.shape != lam.shape:
        raise ValueError("intensity、wavelength、I0 必须是一维等长数组")
    if lam0 is None:
        lam0 = estimate_background_center(lam,bg_all,n_core_points=background_core_points)["lam0"]
    _check_settings(lam0, n_global_runs, base_seed, de_maxiter, de_popsize, max_nfev, w_bounds)
    valid = np.isfinite(obs_all) & np.isfinite(lam) & np.isfinite(bg_all)
    if valid.sum() <= 4:
        raise ValueError("有效波长点必须多于 4 个拟合参数")
    if len(np.unique(lam[valid])) != valid.sum():
        raise ValueError("有效波长点不能重复")
    vel_all = (lam - lam0) / lam0 * C_KM
    vel, obs, bg = vel_all[valid], obs_all[valid], bg_all[valid]
    Imean = float(obs.mean())
    if Imean <= 0:
        raise ValueError("光谱平均强度须为正，以定义 S 的上界")
    scale = max(float(obs.std()), Imean * 1e-3, 1e-12)
    lb = np.array([np.log(1e-3), -200., np.log(w_bounds[0]), 0.])
    ub = np.array([np.log(100.), 200., np.log(w_bounds[1]), 10*Imean])
    args = (vel, obs, bg, scale)
    candidates = []
    for k in range(n_global_runs):
        seed = int(base_seed+k)
        de = differential_evolution(
            _profile_objective, list(zip(lb[:3], ub[:3])), args=(*args, ub[3]),
            maxiter=de_maxiter, popsize=de_popsize, tol=1e-8, atol=1e-12,
            seed=seed, init="sobol", polish=False, vectorized=True,
            updating="deferred", workers=1,
        )
        _, S = _profile_source(de.x, *args, ub[3])
        p0 = np.r_[de.x, S]
        lsq = least_squares(
            _residual, p0, jac=_jacobian, args=args, bounds=(lb, ub),
            loss="linear", x_scale=[1, 50, 1, max(Imean, 1)],
            ftol=1e-10, xtol=1e-10, gtol=1e-10, max_nfev=max_nfev,
        )
        # 保留精炼前的解，防止极端数值情形下局部精炼反而增大 SSE。
        before = float(np.sum(_residual(p0, *args)**2))
        after = float(np.sum(_residual(lsq.x, *args)**2))
        refined = np.isfinite(after) and after <= before
        p = lsq.x if refined else p0
        sse = float(np.sum((single_cloud_vel(p, vel, bg) - obs)**2))
        candidates.append(dict(
            x=p, raw_sse=sse, seed=seed, success=bool(lsq.success and refined),
            message=str(lsq.message) if refined else "Retained DE solution; local refinement did not improve SSE",
            de_success=bool(de.success), de_nit=int(de.nit),
            de_scaled_sse=float(de.fun), lsq_nfev=int(lsq.nfev),
        ))
    candidates.sort(key=lambda c: c["raw_sse"] if np.isfinite(c["raw_sse"]) else np.inf)
    best = candidates[0]
    if not np.isfinite(best["raw_sse"]):
        raise ValueError("未获得有限的拟合候选")
    p = best["x"].copy()
    physical = p.copy()
    physical[[0, 2]] = np.exp(physical[[0, 2]])
    rmse = float(np.sqrt(best["raw_sse"] / valid.sum()))
    near = (p-lb <= 1e-4*(ub-lb)) | (ub-p <= 1e-4*(ub-lb))
    return dict(
        params=physical, log_params=p, model=single_cloud_vel(p, vel_all, bg_all),
        raw_sse=best["raw_sse"], rmse=rmse, nrmse=rmse/Imean,
        scaled_sse=best["raw_sse"]/scale**2, scale=scale, n_valid=int(valid.sum()),
        valid=valid, success=best["success"], message=best["message"], seed=best["seed"],
        near_bound=near, lower_bounds=lb, upper_bounds=ub,
        W_angstrom=physical[2]*lam0/C_KM, candidates=candidates, lam0=lam0,
        elapsed_seconds=perf_counter()-start,
    )


def _plot_pixel(path, lam, obs, bg, result, x, y, lam0, dpi):
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    fig = Figure(figsize=(6.5, 4.5), layout="constrained")
    FigureCanvasAgg(fig)
    ax, res = fig.subplots(2, 1, sharex=True, gridspec_kw={"height_ratios": [3, 1]})
    order = np.argsort(lam)
    ax.plot(lam[order], obs[order], "o", ms=3, color="#0072B2", label="Observed")
    ax.plot(lam[order], bg[order], "--", color="#009E73", label="Background I0")
    if "params" in result:
        model = result["model"]
        ax.plot(lam[order], model[order], color="#D55E00", label="Single-cloud fit")
        residual = np.where(result["valid"], obs-model, np.nan)
        res.plot(lam[order], residual[order], ".-", ms=3, color="#0072B2")
        ax.set_title(f"x={x}, y={y} | v={result['params'][1]:.2f} km/s | RMSE={result['rmse']:.3g}\n"
                     f"Converged={result['success']} | Near bound={result['near_bound'].any()}", fontsize=10)
    else:
        ax.set_title(f"x={x}, y={y}: fit failed; see NPZ / CSV", fontsize=10)
    ax.axvline(lam0, ls=":", color="0.4")
    ax.set_ylabel("Intensity (DN)")
    ax.legend(fontsize=8)
    res.axhline(0, color="0.4", lw=0.8)
    res.set_ylabel("Obs - fit\n(DN)")
    res.set_xlabel(r"Wavelength ($\mathrm{\AA}$)")
    fig.savefig(path, dpi=dpi)
    fig.clear()


def _pixel_job(j, i, obs, lam, bg, settings, out, y_start, x_start, save_plots, dpi):
    # 此函数定义在可导入模块中，适用于 Windows / Jupyter 的 loky 进程。
    x, y = x_start+i, y_start+j
    try:
        result = fit_single_cloud_spectrum(obs, lam, bg, **settings)
    except (ValueError, RuntimeError, FloatingPointError) as exc:
        result = dict(success=False, message=f"{type(exc).__name__}: {exc}",
                      n_valid=int(np.sum(np.isfinite(obs) & np.isfinite(lam) & np.isfinite(bg))))
    saved = {k: v for k, v in result.items() if k != "candidates"}
    if "candidates" in result:
        for key in result["candidates"][0]:
            saved["candidate_"+key] = np.array([c[key] for c in result["candidates"]])
    stem = Path(out) / "pixel_fits" / f"fit_y{y:04d}_x{x:04d}"
    temp = stem.with_suffix(".tmp.npz")
    np.savez_compressed(temp, **saved, wavelength=lam, intensity=obs, I0=bg,
                        y=y, x=x, param_names=PARAM_NAMES)
    os.replace(temp, stem.with_suffix(".npz"))
    if save_plots:
        _plot_pixel(stem.with_suffix(".png"), lam, obs, bg, result, x, y, settings["lam0"], dpi)
    result.pop("candidates", None)
    return j, i, result


def fit_single_cloud_cube(
    data_cube, wavelength, I0, *, output_dir=DEFAULT_OUT, y_start=789, x_start=1515,
    lam0=None, n_jobs=8, n_global_runs=30, base_seed=20260706,
    de_maxiter=500, de_popsize=15, max_nfev=2000, w_bounds=(0.05, 200.0),
    save_plots=True, plot_dpi=180, verbose=True, background_core_points=5,
):
    """逐像素独立拟合 (n_lambda,ny,nx) 数据；完成后自动保存 result.pkl。

    n_jobs 为像素进程数，1=串行，-1=全部可用逻辑核；默认8。
    全局搜索内部向量化、单进程；限制各像素进程的 BLAS 线程为1。
    n_global_runs=30 为默认搜索预算；每像素保存所有候选，供检验搜索稳定性。
    每个像素独立随机种子只取决于位置，并行完成顺序不改变种子。
    所有原始候选均保存在逐像素 NPZ。失败像素保留 NaN 并继续。
    CSV 按完成顺序写入，通过 x/y 定位；数组始终按 [y,x] 排列。
    plot_dpi=180 用于批量检查，可设300；save_plots=False 可跳过逐像素图。
    lam0=None（默认）从共同I0估计一次背景线心；保存background_center.png/npz。
    background_core_points=5：背景最低点及左右各2点；云模型仍使用全部输入波长点。
    """
    from joblib import Parallel, delayed, parallel_config, effective_n_jobs
    cube, lam, bg = map(_array, (data_cube, wavelength, I0))
    if cube.ndim != 3 or min(cube.shape) < 1:
        raise ValueError("data_cube 必须是非空 (n_lambda,ny,nx) 数组")
    if lam.shape != (cube.shape[0],) or bg.shape != lam.shape:
        raise ValueError("wavelength、I0 须与数据波长轴等长")
    calibration = None
    if lam0 is None:
        calibration = estimate_background_center(lam,bg,n_core_points=background_core_points)
        lam0 = calibration["lam0"]
    _check_settings(lam0, n_global_runs, base_seed, de_maxiter, de_popsize, max_nfev, w_bounds)
    if not isinstance(n_jobs, (int, np.integer)) or n_jobs == 0:
        raise ValueError("n_jobs 必须为非零整数")
    if not np.isfinite(plot_dpi) or plot_dpi <= 0:
        raise ValueError("plot_dpi 必须为正数")
    if (np.isfinite(lam) & np.isfinite(bg)).sum() <= 4:
        raise ValueError("共同波长和 I0 有效点必须多于4个")
    finite_lam = lam[np.isfinite(lam) & np.isfinite(bg)]
    if np.unique(finite_lam).size != finite_lam.size:
        raise ValueError("有效波长不能重复")
    ny, nx = cube.shape[1:]
    if base_seed + ny*nx*n_global_runs > 2**32:
        raise ValueError("随机种子超出范围")
    workers = min(effective_n_jobs(n_jobs), ny*nx)
    out = Path(output_dir) / datetime.now().strftime("run_%Y%m%d_%H%M%S_%f")
    (out / "pixel_fits").mkdir(parents=True, exist_ok=False)
    settings = dict(lam0=lam0, n_global_runs=n_global_runs, base_seed=base_seed,
                    de_maxiter=de_maxiter, de_popsize=de_popsize, max_nfev=max_nfev,
                    w_bounds=w_bounds)
    metadata = dict(settings, shape=list(cube.shape), y_start=y_start, x_start=x_start,
                    n_jobs=workers, save_plots=save_plots, plot_dpi=plot_dpi,
                    param_names=PARAM_NAMES, log_param_names=["ltau0", "v", "lW", "S"],
                    model="I=I0*exp(-tau)+S*(1-exp(-tau)); no limb-darkening correction",
                    width="W in km/s; tau=tau0*exp(-((v_grid-v)/W)**2)",
                    zero_point="background_core_quadratic" if calibration is not None else "explicit_lam0",
                    background_core_points=background_core_points if calibration is not None else None,
                    optimizer="Bounded analytical S; vectorized DE in 3 dimensions; analytic-Jacobian 4D least squares")
    (out / "settings.json").write_text(json.dumps(metadata, indent=2, default=lambda a: a.item()), encoding="utf-8")
    np.savez_compressed(out / "input_spectra.npz", data_cube=cube, wavelength=lam, I0=bg)
    if calibration is not None:
        _save_background_center(out,lam,bg,calibration)
    result = {name: np.full((ny, nx), np.nan) for name in
              ("raw_sse", "rmse", "nrmse", "W_angstrom", "elapsed_seconds")}
    result.update(params=np.full((ny, nx, 4), np.nan), log_params=np.full((ny, nx, 4), np.nan),
                  model=np.full(cube.shape, np.nan), success=np.zeros((ny, nx), dtype=bool),
                  near_bound=np.zeros((ny, nx, 4), dtype=bool), n_valid=np.zeros((ny, nx), dtype=int))
    fields = ["y", "x", *PARAM_NAMES, "W_angstrom", "raw_sse", "rmse", "nrmse",
              "n_valid", "success", "near_bound", "seed", "elapsed_seconds", "message"]
    start = perf_counter()
    if verbose:
        print(f"Single cloud: {ny*nx} pixels, {workers} workers, {n_global_runs} searches/pixel\n{out}", flush=True)
        print(f"Shared wavelength zero point: {lam0:.8f} Angstrom ({metadata['zero_point']})", flush=True)
    jobs = (delayed(_pixel_job)(j, i, cube[:, j, i], lam, bg,
                dict(settings, base_seed=base_seed+(j*nx+i)*n_global_runs),
                out, y_start, x_start, save_plots, plot_dpi) for j, i in np.ndindex(ny, nx))
    with (out / "fit_parameters.csv").open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        with parallel_config(backend="loky", inner_max_num_threads=1):
            with Parallel(n_jobs=workers, return_as="generator_unordered", batch_size=1) as parallel:
                for done, (j, i, fit) in enumerate(parallel(jobs), 1):
                    for name, values in result.items():
                        if name in fit:
                            if name == "model":
                                values[:, j, i] = fit[name]
                            else:
                                values[j, i] = fit[name]
                    row = {name: fit.get(name, np.nan) for name in fields}
                    row.update(y=y_start+j, x=x_start+i,
                               **dict(zip(PARAM_NAMES, result["params"][j, i])),
                               near_bound=";".join(np.array(PARAM_NAMES)[result["near_bound"][j, i]]))
                    writer.writerow(row)
                    stream.flush()
                    if verbose:
                        print(f"[{done}/{ny*nx}] y={y_start+j}, x={x_start+i}, "
                              f"v={result['params'][j,i,1]:.3f}, RMSE={result['rmse'][j,i]:.4g}, "
                              f"success={result['success'][j,i]}, elapsed={perf_counter()-start:.1f}s", flush=True)
    result.update(v=result["params"][..., 1], tau0=result["params"][..., 0],
                  W=result["params"][..., 2], S=result["params"][..., 3],
                  wavelength=lam, I0=bg, lam0=lam0, lam0_method=metadata["zero_point"], param_names=np.array(PARAM_NAMES),
                  x_pixels=np.arange(x_start, x_start+nx), y_pixels=np.arange(y_start, y_start+ny),
                  fit_and_pixel_output_seconds=perf_counter()-start)
    np.savez_compressed(out / "fit_maps.npz", **result)
    result["output_dir"] = out
    with (out / "result.pkl").open("wb") as stream:
        pickle.dump(result, stream, protocol=pickle.HIGHEST_PROTOCOL)
    from plot_single_cloud_maps import plot_single_cloud_maps
    fig, _, _ = plot_single_cloud_maps(out, save=True)
    import matplotlib.pyplot as plt
    plt.close(fig)
    if verbose:
        print(f"Saved {out}; converged={result['success'].sum()}/{ny*nx}; total={perf_counter()-start:.1f}s", flush=True)
    return result
