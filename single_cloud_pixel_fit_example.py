"""导入 run_chase_region_single(rsm, lam00) 拟合同一矩形的全部187像素。"""
from single_cloud_pixel_fit import fit_single_cloud_cube


def run_chase_region_single(rsm, lam00, **fit_options):
    return fit_single_cloud_cube(
        rsm[1].data[44:96, 789:806, 1515:1526],
        lam00[44:96],
        rsm[1].data[44:96, 797:826, 1540:1569].mean(axis=(1, 2)),
        y_start=789, x_start=1515, **fit_options,
    )


# 在已有 rsm、lam00 的 notebook 中运行：
# from single_cloud_pixel_fit_example import run_chase_region_single
# result_single = run_chase_region_single(rsm, lam00, n_jobs=8, n_global_runs=30)
# # 默认lam0=None：背景谷底5点局部二次拟合线心，全部像素共用。
# print('Background center:', result_single['lam0'], 'Angstrom')
# print(result_single['output_dir'])
# velocity = result_single['v']    # (17,11), km/s
# params = result_single['params'] # (17,11,4): tau0,v,W,S
#
# 默认30次搜索；可按保存的逐候选结果比较前3/10/30次的最佳参数。
# 查看图：
# from plot_single_cloud_maps import plot_single_cloud_maps
# fig, axes, data = plot_single_cloud_maps(result_single['output_dir'], velocity_limit=50)
#
# 重启后读取自动保存的完整字典：
# import pickle
# with open(r'实际run目录/result.pkl', 'rb') as f:
#     result_single = pickle.load(f)
