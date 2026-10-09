import pickle
from copy import copy
from pathlib import Path

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from matplotlib.colors import Normalize
from matplotlib.patches import ConnectionPatch, Rectangle
import sunpy.map
from sunpy.coordinates import frames

mpl.rcParams.update({'text.usetex': False, 'font.size': 11, 'axes.labelsize': 11, 'xtick.labelsize': 10, 'ytick.labelsize': 10})

# 1. 路径和显示参数；这里只读取结果、画图，不重新拟合。
run_dir = Path(r'C:\Learning\PHD1st\magnetic_reconnecion\data\CHASE_Ha_image\hybrid_cloud\run_20260915_204941_009138')
chase_file = Path(r'C:\Learning\PHD1st\magnetic_reconnecion\data\CHASE_Ha\RSM20240618T211711_0000_HA.fits')
aia_file = Path(r'C:\Learning\PHD1st\magnetic_reconnecion\data\AIA2\131\aia.lev1_euv_12s.2024-06-18T211708Z.131.image_lev1.fits')
save_dir = Path(r'C:\Learning\PHD1st\magnetic_reconnecion\paper\fig6')
save_figure = True
velocity_limit = 70.0  # 对称色标 [-70, 70] km/s；改成 50 即为 [-50, 50]。
contour_level = 800.0  # 这份已保存结果实际采用的阈值。

# 2. 读取“单云 + 一阶矩”的速度，不用双云结果替换其中的像素。
with (run_dir / 'result_hybrid.pkl').open('rb') as f:
    result_hybrid = pickle.load(f)

velocity = np.asarray(result_hybrid['velocity_stage1'], dtype=float)
ha_roi = np.asarray(result_hybrid['ha_image'], dtype=float)
x_pixels = np.asarray(result_hybrid['x_pixels'], dtype=int)
y_pixels = np.asarray(result_hybrid['y_pixels'], dtype=int)
ny, nx = velocity.shape
assert velocity.shape == ha_roi.shape == (y_pixels.size, x_pixels.size)
assert np.all(np.diff(x_pixels) == 1) and np.all(np.diff(y_pixels) == 1)
x0, x1 = int(x_pixels[0]), int(x_pixels[-1]) + 1
y0, y1 = int(y_pixels[0]), int(y_pixels[-1]) + 1
print(f'读取速度图：{velocity.shape}，y={y0}:{y1}，x={x0}:{x1}')

# 3. CHASE 翼部强度图仍用 44:54；这是成像波段，不改变拟合结果。
with fits.open(chase_file) as hdul:
    chase_header = hdul[1].header.copy()
    hawing = np.mean(hdul[1].section[44:54, :, :], axis=0, dtype=np.float64)

lam00 = chase_header['CRVAL3'] + np.arange(chase_header['NAXIS3']) * chase_header['CDELT3']
chase_time = chase_header['DATE_OBS']
coord_HIS = SkyCoord(0 * u.arcsec, 0 * u.arcsec, obstime=chase_time, observer='earth', frame=frames.Helioprojective)

# 保留原图的 CHASE reference_pixel 和 plate scale，使 a、b、c 沿用同一套配准。
headerwing = sunpy.map.make_fitswcs_header(hawing, coord_HIS, reference_pixel=[chase_header['CRPIX1'], chase_header['CRPIX2']] * u.pixel, scale=[0.5218 * 2, 0.5218 * 2] * u.arcsec / u.pixel, telescope='CHASE', instrument='RSM')
hawing_map = sunpy.map.Map(hawing, headerwing)
bottom_left = SkyCoord(Tx=200 * u.arcsec, Ty=-520 * u.arcsec, frame=hawing_map.coordinate_frame)
top_right = SkyCoord(Tx=420 * u.arcsec, Ty=-300 * u.arcsec, frame=hawing_map.coordinate_frame)
sub_hawing_map = hawing_map.submap(bottom_left, top_right=top_right)

# 用同一张 CHASE 图的精确像素裁剪取得 c 图 WCS；submap 的右上像素包含在内。
roi_map = hawing_map.submap([x0, y0] * u.pixel, top_right=[x1 - 1, y1 - 1] * u.pixel)
assert roi_map.data.shape == velocity.shape
assert np.allclose(roi_map.data, ha_roi, rtol=1e-5, atol=1e-4, equal_nan=True), 'FITS 的翼部强度与保存结果不匹配，请核对文件和区域。'

# 4. AIA 131 图像。
aia = sunpy.map.Map(aia_file)
bottom_left = SkyCoord(Tx=200 * u.arcsec, Ty=-520 * u.arcsec, frame=aia.coordinate_frame)
top_right = SkyCoord(Tx=420 * u.arcsec, Ty=-300 * u.arcsec, frame=aia.coordinate_frame)
sub_aia = aia.submap(bottom_left, top_right=top_right)
aia_norm = copy(sub_aia.plot_settings['norm'])
aia_norm.vmin, aia_norm.vmax = 0.01, 2000

# 5. 一行三幅图；c 保持真实像素比例，因此比 a、b 窄。
fig = plt.figure(figsize=(16.5, 5.4), facecolor='white')
gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, nx / ny], left=0.05, right=0.92, bottom=0.16, top=0.96, wspace=0.28)
ax1 = fig.add_subplot(gs[0, 0], projection=sub_aia)
ax2 = fig.add_subplot(gs[0, 1], projection=sub_hawing_map)
ax3 = fig.add_subplot(gs[0, 2], projection=roi_map)
sub_aia.plot(axes=ax1, norm=aia_norm)
sub_hawing_map.plot(axes=ax2, cmap='afmhot', vmin=0, vmax=2 * np.nanmean(sub_hawing_map.data))

velocity_cmap = mpl.colormaps['bwr'].resampled(257).copy()
velocity_cmap.set_bad('0.75')
velocity_norm = Normalize(vmin=-velocity_limit, vmax=velocity_limit)
im = ax3.imshow(np.ma.masked_invalid(velocity), origin='lower', cmap=velocity_cmap, norm=velocity_norm, interpolation='nearest', aspect='equal')
ax3.set_xlim(-0.5, nx - 0.5)
ax3.set_ylim(-0.5, ny - 0.5)

# b 图参考区域黑框，保留原来的位置。
background_bl = SkyCoord(Tx=380 * u.arcsec, Ty=-370 * u.arcsec, frame=hawing_map.coordinate_frame)
background_tr = SkyCoord(Tx=410 * u.arcsec, Ty=-340 * u.arcsec, frame=hawing_map.coordinate_frame)
sub_hawing_map.draw_quadrangle(background_bl, top_right=background_tr, axes=ax2, edgecolor='black', linewidth=1)

# a、b 的白框都沿 c 图的像素外边缘绘制，包含完整的 17 × 11 个像素。
for ax in (ax1, ax2):
    box = Rectangle((-0.5, -0.5), nx, ny, facecolor='none', edgecolor='white', linewidth=1.1, transform=ax.get_transform(roi_map.wcs), zorder=4)
    ax.add_patch(box)

ax2.contour(ha_roi, levels=[contour_level], colors='blue', linewidths=0.5, transform=ax2.get_transform(roi_map.wcs), zorder=5)
ax3.contour(ha_roi, levels=[contour_level], colors='0.3', linewidths=0.9, zorder=3)

# 6. 三幅图的坐标轴统一显示太阳角秒坐标。
for ax in (ax1, ax2, ax3):
    ax.set_title('')
    ax.coords[0].set_format_unit(u.arcsec)
    ax.coords[1].set_format_unit(u.arcsec)
    ax.coords[0].set_major_formatter('s')
    ax.coords[1].set_major_formatter('s')
    ax.coords[0].set_axislabel('X (arcsec)', minpad=0.7)
    ax.coords[1].set_axislabel('Y (arcsec)', minpad=0.7)
    ax.coords[0].set_ticks_position('b')
    ax.coords[1].set_ticks_position('l')
    ax.coords[0].set_ticklabel_position('b')
    ax.coords[1].set_ticklabel_position('l')
    ax.coords[0].set_axislabel_position('b')
    ax.coords[1].set_axislabel_position('l')
    ax.coords.grid(False)

for ax in (ax1, ax2):
    ax.coords[0].set_ticks(spacing=50 * u.arcsec)
    ax.coords[1].set_ticks(spacing=50 * u.arcsec)
ax3.coords[0].set_ticks(spacing=4 * u.arcsec)
ax3.coords[1].set_ticks(spacing=4 * u.arcsec)

# 7. 图题放在各图内部；derived 的拼写和一阶矩方法的英文已调整。
aia_clock = aia.date.strftime('%H:%M:%S')
chase_clock = hawing_map.date.strftime('%H:%M:%S')
ax1.text(0.02, 0.975, f'(a)  AIA 131  {aia_clock} UT', color='white', transform=ax1.transAxes, fontsize=12, va='top', zorder=6)
b_title = rf'(b)  CHASE H$\alpha$ {lam00[44]:.2f}–{lam00[54]:.2f} $\mathrm{{\AA}}$' + f'\n       {chase_clock} UT'
ax2.text(0.02, 0.975, b_title, color='white', transform=ax2.transAxes, fontsize=12, va='top', zorder=6)
c_title = '(c) Doppler velocity derived\nfrom single-cloud modeling\nand moment'
ax3.text(0.025, 0.975, c_title, color='black', transform=ax3.transAxes, fontsize=10, va='top', linespacing=1.2, bbox=dict(facecolor='white', edgecolor='none', alpha=0.8, pad=2), zorder=6)

# 8. 长方形 colorbar；extend='neither' 去掉上下两端的尖角。
fig.canvas.draw()
pos3 = ax3.get_position()
cax = fig.add_axes([pos3.x1 + 0.012, pos3.y0, 0.012, pos3.height])
cbar = fig.colorbar(im, cax=cax, extend='neither')
cbar.set_label(r'Doppler velocity (km s$^{-1}$)')

# 按你的要求：b 中白框左上、左下，分别连到 c 图左上、左下。
for roi_corner, panel_corner in [((-0.5, ny - 0.5), (0, 1)), ((-0.5, -0.5), (0, 0))]:
    connection = ConnectionPatch(xyA=roi_corner, xyB=panel_corner, coordsA=ax2.get_transform(roi_map.wcs), coordsB='axes fraction', axesA=ax2, axesB=ax3, color='black', linewidth=0.9, arrowstyle='-', clip_on=False, zorder=5)
    fig.add_artist(connection)

if save_figure:
    save_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(save_dir / 'aia_chase_single_cloud_moment.pdf', bbox_inches='tight', pad_inches=0.1)
    fig.savefig(save_dir / 'aia_chase_single_cloud_moment.png', dpi=300, bbox_inches='tight', pad_inches=0.1)
    print(f'图片已保存到：{save_dir}')

plt.show()
