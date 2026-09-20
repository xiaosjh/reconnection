"""读取一次逐像素拟合的输出目录，画 H-alpha、vL、vU 三联图。"""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.ticker import MaxNLocator


def plot_double_cloud_triptych(
    run_dir, *, ha_slice=slice(0, 10), ha_image=None,
    velocity_limit=None, contour_level=800, exclude_near_bound=False,
    max_nrmse=None, save=True,
):
    """返回 (fig, axes, data)，不重新拟合。

    run_dir: 包含 fit_maps.npz 和 input_spectra.npz 的 run_时间戳目录。
    ha_slice: 相对于保存的 data_cube 的波长索引。对于输入 44:96，
        slice(0,10) 对应原 FITS 的 44:54，复用原 notebook 的翼部平均。
        这不是对整个拟合波段求平均，也不是单色线心图。
    ha_image: 可选的 (ny,nx) 原始强度图；提供时替代 ha_slice 的平均图。
        必须和拟合区域同范围、同排列，单位 DN。
    velocity_limit: 两张速度图共同的对称色标范围 [-limit,+limit]，km/s。
        None 表示使用有效速度的最大绝对值；固定范围时色条标注超界值。
    contour_level: 在三个面板叠加同一强度等值线；None 表示不画。
        等值线只供定位，不限制速度显示区域。
    exclude_near_bound/max_nrmse: 可选筛选，默认只屏蔽未收敛和非有限值。
        收敛并不等于参数唯一或具有可靠的物理解释。

    图轴为原始数组像素；不会自动做旋转、配准或太阳角秒坐标转换。
    """
    run_dir = Path(run_dir)
    if not (run_dir / "fit_maps.npz").is_file():
        raise FileNotFoundError(
            f"尚未找到 {run_dir / 'fit_maps.npz'}。该汇总文件在全部像素拟合完成后生成；"
            "请等待结束，并确认选择的是具体的 run_时间戳目录。"
        )
    with np.load(run_dir / "fit_maps.npz", allow_pickle=False) as f:
        vL, vU = f["vL"].copy(), f["vU"].copy()
        valid = f["success"].copy()
        x, y = f["x_pixels"].copy(), f["y_pixels"].copy()
        if exclude_near_bound:
            valid &= ~f["near_bound"].any(axis=-1)
        if max_nrmse is not None:
            if not np.isfinite(max_nrmse) or max_nrmse <= 0:
                raise ValueError("max_nrmse 必须是正的有限数")
            valid &= np.isfinite(f["nrmse"]) & (f["nrmse"] <= max_nrmse)
    if ha_image is None:
        with np.load(run_dir / "input_spectra.npz", allow_pickle=False) as f:
            band = f["data_cube"][ha_slice]
            wavelengths = f["wavelength"][ha_slice]
        if band.ndim != 3 or band.shape[0] == 0:
            raise ValueError("ha_slice 应为非空波长切片，例如 slice(0,10)")
        count = np.isfinite(band).sum(axis=0)
        total = np.where(np.isfinite(band), band, 0).sum(axis=0)
        ha_image = np.divide(total, count, out=np.full(total.shape, np.nan), where=count > 0)
        ha_title = (r"(a) H$\alpha$ wing mean" + "\n"
                    + f"{wavelengths.min():.3f}–{wavelengths.max():.3f}" + r" $\mathrm{\AA}$")
    else:
        ha_image = np.asarray(np.ma.asarray(ha_image, dtype=float).filled(np.nan))
        ha_title = r"(a) H$\alpha$ intensity"
    if ha_image.shape != vL.shape or vL.shape != (len(y), len(x)):
        raise ValueError("H-alpha 图像和速度图必须具有相同的 (ny,nx) 形状")
    valid &= np.isfinite(vL) & np.isfinite(vU)
    vL, vU = np.where(valid, vL, np.nan), np.where(valid, vU, np.nan)
    values = np.concatenate((vL[valid], vU[valid]))
    if velocity_limit is None:
        velocity_limit = max(float(np.max(np.abs(values))), 1) if values.size else 1
    if not np.isfinite(velocity_limit) or velocity_limit <= 0:
        raise ValueError("velocity_limit 必须是正的有限数")
    norm = Normalize(vmin=-velocity_limit, vmax=velocity_limit)
    extent = (x[0] - 0.5, x[-1] + 0.5, y[0] - 0.5, y[-1] + 0.5)
    fig = plt.figure(figsize=(10, 5.4), layout="constrained")
    grid = fig.add_gridspec(2, 3, height_ratios=[1, 0.045])
    axes = np.array([fig.add_subplot(grid[0, 0])])
    axes = np.append(axes, [fig.add_subplot(grid[0, k], sharex=axes[0], sharey=axes[0])
                            for k in (1, 2)])
    ha_cmap = plt.get_cmap("afmhot").copy()
    vel_cmap = plt.get_cmap("RdBu_r").copy()
    ha_cmap.set_bad("0.85")
    vel_cmap.set_bad("0.85")
    image_kwargs = dict(origin="lower", extent=extent, interpolation="nearest", aspect="equal")
    im_ha = axes[0].imshow(ha_image, cmap=ha_cmap, **image_kwargs)
    im_vel = axes[1].imshow(vL, cmap=vel_cmap, norm=norm, **image_kwargs)
    axes[2].imshow(vU, cmap=vel_cmap, norm=norm, **image_kwargs)
    titles = (ha_title, r"(b) Lower cloud: $v_L$", r"(c) Upper cloud: $v_U$")
    finite_ha = ha_image[np.isfinite(ha_image)]
    draw_contour = (contour_level is not None and min(ha_image.shape) >= 2
                    and finite_ha.size and finite_ha.min() < contour_level < finite_ha.max())
    for k, ax in enumerate(axes):
        ax.set_title(titles[k], fontsize=11)
        ax.set_xlabel("X (pixel)")
        ax.xaxis.set_major_locator(MaxNLocator(nbins=4, integer=True))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=6, integer=True))
        if k:
            ax.tick_params(labelleft=False)
        if draw_contour:
            ax.contour(x, y, ha_image, levels=[contour_level], colors="green", linewidths=0.8)
        ax.set_xlim(extent[:2])
        ax.set_ylim(extent[2:])
    axes[0].set_ylabel("Y (pixel)")
    fig.colorbar(im_ha, cax=fig.add_subplot(grid[1, 0]), orientation="horizontal", label="Intensity (DN)")
    below, above = np.any(values < -velocity_limit), np.any(values > velocity_limit)
    extend = "both" if below and above else "min" if below else "max" if above else "neither"
    fig.colorbar(im_vel, cax=fig.add_subplot(grid[1, 1:]), orientation="horizontal",
                 label="Doppler velocity (km/s); positive = redshift", extend=extend)
    if save:
        fig.savefig(run_dir / "ha_double_cloud_triptych.png", dpi=300)
        fig.savefig(run_dir / "ha_double_cloud_triptych.pdf")
    return fig, axes, dict(ha_image=ha_image, vL=vL, vU=vU, valid=valid, x_pixels=x, y_pixels=y)
