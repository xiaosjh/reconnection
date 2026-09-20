"""将下面代码复制到已有 rsm、lam00 的 notebook 单元中运行。"""

from double_cloud_pixel_fit import fit_double_cloud_cube


def run_chase_region(rsm, lam00, **fit_options):
    """使用现有 FITS 对象和已有波长标定，拟合指定的全部 187 个像素。"""
    return fit_double_cloud_cube(data_cube=rsm[1].data[44:96, 789:806, 1515:1526], wavelength=lam00[44:96], I0=rsm[1].data[44:96, 797:826, 1540:1569].mean(axis=(1, 2)), y_start=789, x_start=1515, **fit_options)


# 在 notebook 中：
# from double_cloud_pixel_fit_example import run_chase_region
# result = run_chase_region(rsm, lam00)  # 默认16个像素进程，每像素30次全局搜索
# print(result['output_dir'])
# vL_map = result['vL']                 # shape = (17, 11)，单位 km/s
# vU_map = result['vU']
# params = result['params']             # shape = (17, 11, 8)
# # 参数顺序：tau0L, tau0U, vL, vU, WL, WU, SL, SU
# # 原图 (y=795, x=1521) 的参数：params[795-789, 1521-1515]
#
# 首次可用单像素检查耗时（降低搜索次数只用于试运行）：
# result_test = fit_double_cloud_cube(rsm[1].data[44:96, 795:796, 1521:1522], lam00[44:96], rsm[1].data[44:96, 797:826, 1540:1569].mean(axis=(1, 2)), y_start=795, x_start=1521, n_global_runs=1)
#
# 对已有 hybrid Notebook：保持 double_mask 与原来的外层 Parallel 循环，把 n_jobs 改为16。
# fit_double_cloud_spectrum 自动启用向量化和解析导数，不在单像素内部再开进程。
# 全矩形入口支持 save_plots=False 跳过绘图，但仍保存逐像素NPZ、CSV和汇总NPZ。
