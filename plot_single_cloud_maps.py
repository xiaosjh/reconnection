"""从单云拟合结果目录读取 H-alpha 图像及单云速度图，无需重新拟合。"""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator


def plot_single_cloud_maps(run_dir, *, ha_slice=slice(0, 10), velocity_limit=None,
                           contour_level=800, exclude_near_bound=False, save=True):
    """ha_slice 相对保存的光谱起点；原输入44:96时，0:10对应FITS的44:54。

    默认只屏蔽未收敛值，near_bound 可额外筛选；不做空间平滑或重采样。
    返回 fig, axes, data；save=True 保存 PNG/PDF 到同一 run 目录。
    """
    run_dir = Path(run_dir)
    with np.load(run_dir / "fit_maps.npz", allow_pickle=False) as f:
        velocity = f["v"].copy()
        valid = f["success"] & np.isfinite(velocity)
        if exclude_near_bound:
            valid &= ~f["near_bound"].any(axis=-1)
        x, y = f["x_pixels"], f["y_pixels"]
        lam0 = float(f["lam0"])
        reference = str(f["lam0_method"]) if "lam0_method" in f else "explicit_lam0"
    with np.load(run_dir / "input_spectra.npz", allow_pickle=False) as f:
        band = f["data_cube"][ha_slice]
        wavelengths = f["wavelength"][ha_slice]
    if band.ndim != 3 or not band.shape[0]:
        raise ValueError("ha_slice 必须是非空波长切片")
    count = np.isfinite(band).sum(axis=0)
    total = np.where(np.isfinite(band), band, 0).sum(axis=0)
    ha = np.divide(total, count, out=np.full(total.shape, np.nan), where=count > 0)
    velocity[~valid] = np.nan
    values = velocity[valid]
    if velocity_limit is None:
        velocity_limit = max(float(np.max(np.abs(values))), 1.) if values.size else 1.
    if not np.isfinite(velocity_limit) or velocity_limit <= 0:
        raise ValueError("velocity_limit 必须为正数")
    fig, axes = plt.subplots(1, 2, figsize=(7, 5.3), sharex=True, sharey=True, layout="constrained")
    extent = (x[0]-.5, x[-1]+.5, y[0]-.5, y[-1]+.5)
    ha_cmap, vel_cmap = plt.get_cmap("afmhot").copy(), plt.get_cmap("RdBu_r").copy()
    ha_cmap.set_bad("0.85")
    vel_cmap.set_bad("0.85")
    kwargs = dict(origin="lower", extent=extent, interpolation="nearest", aspect="equal")
    im0 = axes[0].imshow(ha, cmap=ha_cmap, **kwargs)
    im1 = axes[1].imshow(velocity, cmap=vel_cmap, vmin=-velocity_limit, vmax=velocity_limit, **kwargs)
    axes[0].set_title(r"(a) H$\alpha$ wing mean"+f"\n{wavelengths.min():.3f}–{wavelengths.max():.3f}"+r" $\mathrm{\AA}$", fontsize=10)
    axes[1].set_title("(b) Single-cloud velocity", fontsize=10)
    reference_label = "Background center" if reference.startswith("background_") else "Reference wavelength"
    fig.suptitle(f"{reference_label}: {lam0:.6f} Angstrom", fontsize=10)
    finite = ha[np.isfinite(ha)]
    for ax in axes:
        if contour_level is not None and min(ha.shape) >= 2 and finite.size and finite.min() < contour_level < finite.max():
            ax.contour(x, y, ha, levels=[contour_level], colors="white", linewidths=.8)
        ax.set_xlabel("X (pixel)")
        ax.xaxis.set_major_locator(MaxNLocator(nbins=4, integer=True))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=6, integer=True))
        ax.set_xlim(extent[:2])
        ax.set_ylim(extent[2:])
    axes[0].set_ylabel("Y (pixel)")
    fig.colorbar(im0, ax=axes[0], orientation="horizontal", label="Intensity (DN)", shrink=.9)
    below, above = np.any(values < -velocity_limit), np.any(values > velocity_limit)
    extend = "both" if below and above else "min" if below else "max" if above else "neither"
    fig.colorbar(im1, ax=axes[1], orientation="horizontal", label="Velocity (km/s); + = redshift", shrink=.9, extend=extend)
    if save:
        fig.savefig(run_dir / "ha_single_cloud_maps.png", dpi=300)
        fig.savefig(run_dir / "ha_single_cloud_maps.pdf")
    return fig, axes, dict(ha_image=ha, v=velocity, valid=valid, x_pixels=x, y_pixels=y)
