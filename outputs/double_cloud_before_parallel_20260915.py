"""CHASE H-alpha 双云模型逐像素拟合；输入顺序为 (wavelength, y, x)。

用法见 double_cloud_pixel_fit_example.py。只依赖 numpy/scipy/matplotlib。
保持原模型的层序：背景 -> L（下层）-> U（上层）-> 观测者。
"""

from pathlib import Path
from datetime import datetime
import csv
import json

import numpy as np
from scipy.optimize import differential_evolution, least_squares
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg


PARAM_NAMES = ("tau0L", "tau0U", "vL", "vU", "WL", "WU", "SL", "SU")
LOG_PARAM_NAMES = ("ltauL", "ltauU", "vL", "vU", "lWL", "lWU", "SL", "SU")
DEFAULT_OUT = Path(r"C:\Learning\PHD1st\magnetic_reconnecion\data\CHASE_Ha_image\double_cloud")


def _float_array(value):
    return np.asarray(np.ma.asarray(value, dtype=float).filled(np.nan))


def double_cloud_vel(params, vel, I0):
    """原始八参数模型；W 是 exp[-((v-v0)/W)^2] 中的宽度，单位 km/s。"""
    ltauL, ltauU, vL, vU, lWL, lWU, SL, SU = params
    tauL = np.exp(ltauL) * np.exp(-((vel - vL) / np.exp(lWL)) ** 2)
    tauU = np.exp(ltauU) * np.exp(-((vel - vU) / np.exp(lWU)) ** 2)
    return (I0 * np.exp(-(tauL + tauU))
            + SL * (-np.expm1(-tauL)) * np.exp(-tauU)
            + SU * (-np.expm1(-tauU)))


def _residual(params, vel, intensity, background, scale):
    return (double_cloud_vel(params, vel, background) - intensity) / scale


def _objective(params, *args):
    residual = _residual(params, *args)
    return float(residual @ residual)


def fit_double_cloud_spectrum(
    intensity, wavelength, I0, *, lam0=6562.8, n_global_runs=30,
    base_seed=20260706, de_maxiter=1500, de_popsize=25,
    max_nfev=10000, w_bounds=(0.05, 200.0),
):
    """独立拟合一条光谱，返回最小原始 SSE 候选及所有候选。

    wavelength/lam0: Angstrom；I0 与 intensity 必须同单位、同波长采样。
    scale 仅为原代码的残差归一化量，不是测量噪声，不能当作 chi-square。
    有限次搜索不保证数学上的全局最优，也不保证参数唯一。
    """
    intensity, wavelength, I0 = map(_float_array, (intensity, wavelength, I0))
    if intensity.ndim != 1 or intensity.shape != wavelength.shape or I0.shape != intensity.shape:
        raise ValueError("intensity、wavelength、I0 必须是相同长度的一维数组")
    if not np.isfinite(lam0) or lam0 <= 0:
        raise ValueError("lam0 必须为正的有限数")
    for name, value in (("n_global_runs", n_global_runs), ("de_maxiter", de_maxiter),
                        ("de_popsize", de_popsize), ("max_nfev", max_nfev)):
        if not isinstance(value, (int, np.integer)) or value < 1:
            raise ValueError(f"{name} 必须为正整数")
    if len(w_bounds) != 2 or not np.all(np.isfinite(w_bounds)) or not 0 < w_bounds[0] < w_bounds[1]:
        raise ValueError("w_bounds 必须满足 0 < W_MIN < W_MAX")
    vel = (wavelength - lam0) / lam0 * 299792.458
    valid = np.isfinite(vel) & np.isfinite(intensity) & np.isfinite(I0)
    if valid.sum() <= 8:
        raise ValueError("有效波长点必须多于 8 个拟合参数")
    v, obs, bg = vel[valid], intensity[valid], I0[valid]
    if np.unique(v).size != v.size:
        raise ValueError("有效波长点不能重复")
    mean_intensity = float(obs.mean())
    if mean_intensity <= 0:
        raise ValueError("有效光谱平均强度必须为正，以定义源函数上界")
    scale = float(np.std(obs - np.median(obs)))
    if not np.isfinite(scale) or scale <= 0:
        scale = max(float(np.mean(np.abs(obs))) * 1e-3, 1e-12)
    lb = np.array([np.log(1e-3), np.log(1e-3), -200, -200,
                   np.log(w_bounds[0]), np.log(w_bounds[0]), 0, 0])
    ub = np.array([np.log(1e2), np.log(1e2), 200, 200,
                   np.log(w_bounds[1]), np.log(w_bounds[1]),
                   mean_intensity * 10, mean_intensity * 10])
    args = (v, obs, bg, scale)
    candidates = []
    for run in range(n_global_runs):
        seed = int(base_seed + run)
        de = differential_evolution(
            _objective, list(zip(lb, ub)), args=args, strategy="best1bin",
            maxiter=de_maxiter, popsize=de_popsize, tol=1e-8, atol=0,
            mutation=(0.5, 1.0), recombination=0.7, seed=seed,
            init="sobol", polish=False, updating="immediate", workers=1,
        )
        lsq = least_squares(
            _residual, de.x, args=args, bounds=(lb, ub), loss="linear",
            x_scale=[1, 1, 50, 50, 1, 1, max(mean_intensity, 1), max(mean_intensity, 1)],
            ftol=1e-12, xtol=1e-12, gtol=1e-12, max_nfev=max_nfev,
        )
        raw_residual = double_cloud_vel(lsq.x, v, bg) - obs
        candidates.append(dict(
            seed=seed, x=lsq.x, raw_sse=float(raw_residual @ raw_residual),
            success=bool(lsq.success), message=str(lsq.message),
            de_success=bool(de.success), de_message=str(de.message),
            de_fun=float(de.fun), lsq_nfev=int(lsq.nfev),
        ))
    candidates.sort(key=lambda item: item["raw_sse"] if np.isfinite(item["raw_sse"]) else np.inf)
    best = candidates[0]
    if not np.isfinite(best["raw_sse"]) or not np.all(np.isfinite(best["x"])):
        raise ValueError("未获得有限的拟合候选")
    params = best["x"].copy()
    params[[0, 1, 4, 5]] = np.exp(params[[0, 1, 4, 5]])
    # 相对整个搜索范围的 0.01% 内，标记为接近边界；不自动剔除。
    near_bound = ((best["x"] - lb) <= 1e-4 * (ub - lb)) | ((ub - best["x"]) <= 1e-4 * (ub - lb))
    rmse = np.sqrt(best["raw_sse"] / valid.sum())
    return dict(
        params=params, log_params=best["x"], model=double_cloud_vel(best["x"], vel, I0),
        valid=valid, n_valid=int(valid.sum()), raw_sse=best["raw_sse"],
        scaled_sse=best["raw_sse"] / scale**2, scale=scale,
        rmse=float(rmse), nrmse=float(rmse / mean_intensity),
        success=best["success"], message=best["message"], seed=best["seed"],
        near_bound=near_bound, lower_bounds=lb, upper_bounds=ub, candidates=candidates,
    )


def _plot_spectrum(path, wavelength, intensity, I0, result, x, y, lam0, error):
    fig = Figure(figsize=(7, 4.8), layout="constrained")
    FigureCanvasAgg(fig)
    ax, residual_ax = fig.subplots(2, 1, sharex=True, gridspec_kw={"height_ratios": [3, 1]})
    order = np.argsort(wavelength)
    ax.plot(wavelength[order], intensity[order], "o", ms=3, color="#0072B2", label="Observed")
    ax.plot(wavelength[order], I0[order], "--", color="#009E73", label="Background I0")
    if result is not None:
        model = result["model"]
        ax.plot(wavelength[order], model[order], color="#D55E00", label="Double-cloud fit")
        residual = np.where(result["valid"], intensity - model, np.nan)
        residual_ax.plot(wavelength[order], residual[order], ".-", ms=3, color="#0072B2")
        p = result["params"]
        ax.set_title(f"Pixel x={x}, y={y} | vL={p[2]:.2f}, vU={p[3]:.2f} km/s\n"
                     f"RMSE={result['rmse']:.3g} | converged={result['success']} | "
                     f"near bound={result['near_bound'].any()}", fontsize=10)
    else:
        ax.set_title(f"Pixel x={x}, y={y}: fit failed", fontsize=10)
        residual_ax.text(0.01, 0.5, "No valid fit; see CSV / NPZ for the reason.",
                         transform=residual_ax.transAxes, fontsize=8)
    ax.axvline(lam0, ls=":", color="0.3")
    ax.set_ylabel("Intensity (DN)")
    ax.legend(fontsize=8)
    residual_ax.axhline(0, color="0.4", lw=0.8)
    residual_ax.set_ylabel("Obs - fit\n(DN)")
    residual_ax.set_xlabel(r"Wavelength ($\mathrm{\AA}$)")
    fig.savefig(path, dpi=300)
    fig.clear()


def _plot_maps(out, params, success, x_start, y_start):
    ny, nx = success.shape
    velocities = np.where(success[..., None], params[..., 2:4], np.nan)
    finite = np.abs(velocities[np.isfinite(velocities)])
    vmax = max(float(finite.max()), 1.0) if finite.size else 1.0
    fig = Figure(figsize=(7, 4.8), layout="constrained")
    FigureCanvasAgg(fig)
    axes = fig.subplots(1, 2, sharex=True, sharey=True)
    from matplotlib import colormaps
    cmap = colormaps["RdBu_r"].copy()
    cmap.set_bad("0.85")
    for k, ax in enumerate(axes):
        im = ax.imshow(velocities[..., k], origin="lower", interpolation="nearest",
                       extent=(x_start - 0.5, x_start + nx - 0.5,
                               y_start - 0.5, y_start + ny - 0.5),
                       cmap=cmap, vmin=-vmax, vmax=vmax)
        ax.set_title(("Lower cloud: vL", "Upper cloud: vU")[k])
        ax.set_xlabel("X (pixel)")
    axes[0].set_ylabel("Y (pixel)")
    fig.colorbar(im, ax=list(axes), shrink=0.8, label="Doppler velocity (km/s)")
    fig.suptitle("Positive = redshift; gray = failed / not converged", fontsize=10)
    fig.savefig(out / "doppler_maps.png", dpi=300)
    fig.savefig(out / "doppler_maps.pdf")
    fig.clear()


def fit_double_cloud_cube(
    data_cube, wavelength, I0, *, output_dir=DEFAULT_OUT,
    y_start=789, x_start=1515, lam0=6562.8, n_global_runs=30,
    base_seed=20260706, de_maxiter=1500, de_popsize=25,
    max_nfev=10000, w_bounds=(0.05, 200.0), verbose=True,
):
    """拟合矩形内全部像素，不应用 contour mask，也不做空间平均。

    data_cube: (n_lambda, ny, nx)，例如 rsm[1].data[44:96,789:806,1515:1526]。
    wavelength: (n_lambda,)，例如 lam00[44:96]，单位 Angstrom。
    I0: (n_lambda,)，同一波长切片的参考区域平均光谱。
    y_start/x_start: 切片在原图中的起始像素，用于命名和图轴。

    每次调用在 output_dir 下创建独立时间戳目录，避免覆盖历史结果。
    每完成一个像素即写 NPZ、CSV 和 PNG；异常像素记录失败后继续。
    返回 params[ny,nx,8]、vL/vU[ny,nx]、success 和 output_dir 等。
    success 仅表示最优候选的局部优化收敛，不表示物理可信度。
    速度正值为红移、负值为蓝移，相对于 lam0，不额外校准零点。
    """
    cube, wavelength, I0 = map(_float_array, (data_cube, wavelength, I0))
    if cube.ndim != 3 or min(cube.shape) < 1:
        raise ValueError("data_cube 必须是非空的 (n_lambda, ny, nx) 三维数组")
    if wavelength.shape != (cube.shape[0],) or I0.shape != wavelength.shape:
        raise ValueError("wavelength 和 I0 必须与 data_cube 的波长维长度一致")
    settings = dict(lam0=lam0, n_global_runs=n_global_runs, base_seed=base_seed,
                    de_maxiter=de_maxiter, de_popsize=de_popsize,
                    max_nfev=max_nfev, w_bounds=w_bounds)
    # 在开始批处理前检查优化设置，避免将配置错误误记为逐像素失败。
    for name in ("n_global_runs", "de_maxiter", "de_popsize", "max_nfev"):
        if not isinstance(settings[name], (int, np.integer)) or settings[name] < 1:
            raise ValueError(f"{name} 必须为正整数")
    if not np.isfinite(lam0) or lam0 <= 0 or len(w_bounds) != 2 or not np.all(np.isfinite(w_bounds)) or not 0 < w_bounds[0] < w_bounds[1]:
        raise ValueError("检查 lam0 和 w_bounds")
    if not isinstance(base_seed, (int, np.integer)) or not 0 <= base_seed < 2**32:
        raise ValueError("base_seed 必须是 [0, 2**32) 内的整数")
    ny, nx = cube.shape[1:]
    if base_seed + ny * nx * n_global_runs > 2**32:
        raise ValueError("随机种子范围超出上限")
    out = Path(output_dir) / datetime.now().strftime("run_%Y%m%d_%H%M%S_%f")
    (out / "pixel_fits").mkdir(parents=True, exist_ok=False)
    metadata = dict(settings, shape=list(cube.shape), y_start=y_start, x_start=x_start,
                    param_names=PARAM_NAMES, log_param_names=LOG_PARAM_NAMES,
                    reference="Shared input I0; no spatial mask", velocity_unit="km/s",
                    width_definition="tau=tau0*exp(-((v-v0)/W)**2)",
                    map_axes="Original array pixels, not solar coordinates",
                    quality="success is local convergence only; scale is not measurement noise")
    (out / "settings.json").write_text(json.dumps(metadata, indent=2, default=lambda v: v.item()), encoding="utf-8")
    np.savez_compressed(out / "input_spectra.npz", data_cube=cube, wavelength=wavelength, I0=I0)
    params = np.full((ny, nx, 8), np.nan)
    log_params = np.full_like(params, np.nan)
    model = np.full_like(cube, np.nan)
    sse, rmse, nrmse = (np.full((ny, nx), np.nan) for _ in range(3))
    success = np.zeros((ny, nx), dtype=bool)
    near_bound = np.zeros((ny, nx, 8), dtype=bool)
    fields = ["y", "x", *PARAM_NAMES, "raw_sse", "rmse", "nrmse", "n_valid",
              "success", "near_bound", "seed", "message"]
    with (out / "fit_parameters.csv").open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for j, i in np.ndindex(ny, nx):
            y, x = y_start + j, x_start + i
            stem = out / "pixel_fits" / f"fit_y{y:04d}_x{x:04d}"
            result, error = None, ""
            try:
                pixel_settings = dict(settings, base_seed=base_seed + (j * nx + i) * n_global_runs)
                result = fit_double_cloud_spectrum(cube[:, j, i], wavelength, I0, **pixel_settings)
            except (ValueError, RuntimeError, FloatingPointError) as exc:
                error = f"{type(exc).__name__}: {exc}"
            if result is not None:
                params[j, i], log_params[j, i] = result["params"], result["log_params"]
                model[:, j, i] = result["model"]
                sse[j, i], rmse[j, i], nrmse[j, i] = result["raw_sse"], result["rmse"], result["nrmse"]
                success[j, i], near_bound[j, i] = result["success"], result["near_bound"]
                row = dict(y=y, x=x, **dict(zip(PARAM_NAMES, result["params"])),
                           **{k: result[k] for k in ("raw_sse", "rmse", "nrmse", "n_valid", "success", "seed", "message")},
                           near_bound=";".join(np.array(PARAM_NAMES)[result["near_bound"]]))
                saved = {k: v for k, v in result.items() if k != "candidates"}
                for key in result["candidates"][0]:
                    saved["candidate_" + key] = np.array([c[key] for c in result["candidates"]])
            else:
                row = dict.fromkeys(fields, np.nan)
                row.update(y=y, x=x, success=False, near_bound="", message=error,
                           n_valid=int(np.sum(np.isfinite(cube[:, j, i]) & np.isfinite(wavelength) & np.isfinite(I0))))
                saved = dict(success=False, message=error)
            np.savez_compressed(stem.with_suffix(".npz"), **saved, wavelength=wavelength,
                                intensity=cube[:, j, i], I0=I0, y=y, x=x, param_names=PARAM_NAMES)
            writer.writerow(row)
            stream.flush()
            _plot_spectrum(stem.with_suffix(".png"), wavelength, cube[:, j, i], I0, result, x, y, lam0, error)
            if verbose:
                print(f"[{j * nx + i + 1}/{ny * nx}] y={y}, x={x}, "
                      f"success={success[j, i]}, RMSE={rmse[j, i]:.5g}", flush=True)
    maps = dict(params=params, log_params=log_params, vL=params[..., 2], vU=params[..., 3],
                raw_sse=sse, rmse=rmse, nrmse=nrmse, success=success, near_bound=near_bound,
                model=model, wavelength=wavelength, I0=I0, param_names=np.array(PARAM_NAMES),
                y_pixels=np.arange(y_start, y_start + ny), x_pixels=np.arange(x_start, x_start + nx))
    np.savez_compressed(out / "fit_maps.npz", **maps)
    _plot_maps(out, params, success, x_start, y_start)
    if verbose:
        print(f"Saved: {out}\nConverged: {success.sum()}/{ny * nx}", flush=True)
    return dict(maps, output_dir=out)
