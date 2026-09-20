# CHASE H-alpha 单云逐像素拟合

在已有 `rsm` 和 `lam00` 的 notebook 中：

```python
from single_cloud_pixel_fit import fit_single_cloud_cube

result_single = fit_single_cloud_cube(
    rsm[1].data[44:96, 789:806, 1515:1526],
    lam00[44:96],
    rsm[1].data[44:96, 797:826, 1540:1569].mean(axis=(1, 2)),
    n_jobs=8,
    n_global_runs=30,
    lam0=None,                  # 自动用参考背景线心作为共同速度零点
    background_core_points=5,   # 仅用背景谷底附近5个波长点定位参考线心
)
```

模型为 `I=I0*exp(-tau)+S*(1-exp(-tau))`，不做临边昏暗修正。
`tau=tau0*exp(-((v_grid-v)/W)**2)`，其中 `v_grid=(lambda-lam0)/lam0*c`。
物理参数顺序 `[tau0,v,W,S]`，速度和宽度 W 单位均为 km/s。
论文截图的 W 是波长单位，转换结果另存于 `W_angstrom`。
当前已按用户要求恢复五点局部二次拟合。`lam0=None` 时，先在共同的实测背景平均谱中
找到最低点，再选择该点及左右各2个波长点拟合抛物线，用其顶点作为参考线心。
全部像素共用这一参考值；每个像素的单云拟合仍使用44:96全部52点，I0仍为实测平均谱。
这里的五点是五个波长采样，不是五个空间像素。背景诊断图中的完整蓝色曲线是
实测数据的连线，不是拟合曲线；绿色曲线才是仅在局部窗口内计算的二次拟合。
不能把这种局部拟合的残差与整个背景波段的高斯拟合残差直接比较。

保留 `n_global_runs=30`。`background_core_points=5` 已恢复；删除先前调用中的
`background_fit_slice`。如已import模块，请先 `importlib.reload(single_cloud_pixel_fit)`。
也可显式传入标定值 `lam0`；自动拟合失败会抛错，不回退到固定波长。
`background_center.png/pdf/npz` 保存局部诊断及实际使用的五点、二次多项式系数。
实际线心为 `result_single['lam0']`，这次数据恢复为6562.872972907947 Angstrom。

恢复五点方法并保留30次搜索后的结果目录：
`C:\Learning\PHD1st\magnetic_reconnecion\data\CHASE_Ha_image\single_cloud\run_20260909_193002_558142`。
已完成187像素，拟合及逐像素输出约27.47秒，全部报告局部收敛。
五个参考点对应原始波长索引69:74，即69、70、71、72、73。
五点仅用于背景线心的局部二次拟合。蓝色完整背景曲线是实测数据的连线，不能作为整条谱线拟合的质量证据。

论文核对：[Qiu et al. 2024，第3.1节](https://arxiv.org/html/2401.16730v1#S3.SS1)。
该节把lambda0定义为背景区域平均谱线心，待拟合参数为S、tau0、v、W。
没有指定用背景区云模型拟合获得lambda0，也没有指定采用五点算法。
由公式可知tau的中心lambda_c=lambda0*(1+v/c)，必须给定参考零点才能解释速度。

代码需要 numpy、scipy>=1.9、matplotlib、joblib>=1.4。
如当前 notebook 内核提示缺少 joblib，可在单元运行 `%pip install "joblib>=1.4"`。
已在本次测试使用的 Python312 环境安装 joblib。

## 运行与输出

- `n_jobs=8`：8个像素进程；1为串行，-1使用全部可用CPU。
- `n_global_runs=30`：默认全局搜索次数；不保证有限搜索的全局最优。
- `de_maxiter=500`：每次差分进化最大代数。
- `save_plots=True`：每像素保存拟合/残差PNG；False跳过逐像素图。
- `plot_dpi=180`：默认诊断图分辨率；可设300。
- 每像素随机种子固定，与进程完成顺序无关。

输出默认位于 `C:\Learning\PHD1st\magnetic_reconnecion\data\CHASE_Ha_image\single_cloud\run_时间戳`。

`pixel_fits/` 保存各像素PNG和NPZ（含每次搜索候选）；CSV按完成顺序写入，使用x/y定位。
`input_spectra.npz` 保存输入立方体、波长、参考谱。
`fit_maps.npz` 保存所有参数/速度/质量数组；完成后保存完整 `result.pkl`。
`ha_single_cloud_maps.png/pdf` 为H-alpha翼部强度和单云速度图，使用原始像素坐标。
所有像素独立拟合，不做空间平均、平滑或只保留contour内的像素。

```python
velocity = result_single['v']       # (17,11), km/s
params = result_single['params']    # (17,11,4): tau0,v,W,S
success = result_single['success']

from plot_single_cloud_maps import plot_single_cloud_maps
fig, axes, data = plot_single_cloud_maps(
    result_single['output_dir'], velocity_limit=50, contour_level=800
)

from check_single_cloud_searches import compare_search_budgets
report = compare_search_budgets(result_single['output_dir'], budgets=(3,10,30))
# 读取保存的候选，生成 search_stability.csv/npz/json，不重新拟合。
```

`success`仅表示局部数值收敛，不等于模型充分解释数据。
RMSE=sqrt(raw_sse/n_valid)，单位与输入强度一致。并非速度误差。
`near_bound` 为参数接近边界的标志，并不自动剔除。
图中默认H-alpha翼部使用保存数据的前10个波长点，即FITS的44:54；它不包含索引54。

## 初版固定零点的实测记录

输入文件 `C:\Learning\PHD1st\magnetic_reconnecion\data\CHASE_Ha\RSM20240618T211711_0000_HA.fits`，
波长按旧 notebook 的 `CRVAL3+arange(NAXIS3)*CDELT3` 构造，使用44:96。
结果目录：`C:\Learning\PHD1st\magnetic_reconnecion\data\CHASE_Ha_image\single_cloud\run_20260909_184451_674713`。
187像素，8进程，3次搜索，拟合和逐像素输出约6.47秒，连同汇总图约7.3秒（单次实测）。
187个局部优化收敛，1个像素有参数接近边界；RMSE中位数13.54 DN，最大44.80 DN。
最大RMSE像素(y=792,x=1522)仍有系统性残差，不能将收敛等同于拟合质量良好。
抽查4像素改用10次搜索，速度变化均小于0.001 km/s；不代表全区域的参数不确定度。

提速包括模型由8参数降为4参数、默认搜索预算减少，以及以下实现优化：
固定tau0/v/W时解析求有界S，全球搜索仅3维；种群批量向量化；局部精炼使用解析导数；
像素间loky多进程，各进程的BLAS线程限制为1。
相同目标、种子与迭代设置的单像素全局搜索对照中，向量化耗时0.023秒，
逐候选计算0.108秒，目标值一致（单次测试，不是整个程序的加速比）。

验证：`python -m unittest test_single_cloud_pixel_fit -v`。
测试覆盖论文公式、解析Jacobian、S的有界解析解、已知参数恢复、无效数据、串/并行一致性与结果保存。

实现参考：
[SciPy vectorized differential_evolution](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.differential_evolution.html)，
[joblib Parallel](https://joblib.readthedocs.io/en/stable/generated/joblib.Parallel.html)，
[joblib parallel_config](https://joblib.readthedocs.io/en/stable/generated/joblib.parallel_config.html)。

## CHASE高斯背景参考 + 30次全局搜索的实测记录

方法来源：[CHASE workshop Python 第4.1节](https://ssdc.nju.edu.cn/web-service/software/Demo/CHASE_workshop_Python.html)。
背景区域仍为y=797:826,x=1540:1569，背景拟合波长索引50:90，使用40个采样点。
参考线心6562.84973863939 Angstrom，背景高斯拟合RMSE=22.5853 DN。
背景线心形状存在系统性残差，该高斯中心是一种模型定义的参考值，不等于绝对无误差线心。
背景窗口由50:90扩为48:92或缩为52:88，线心分别为6562.850511和6562.848786 Angstrom。

完整187像素、8进程、30次搜索，拟合及逐像素输出耗时26.60秒。
目录：`C:\Learning\PHD1st\magnetic_reconnecion\data\CHASE_Ha_image\single_cloud\run_20260909_192223_353462`。
全部像素报告局部收敛，已保存187张拟合图、完整result.pkl、背景诊断和二维图。

按同一组嵌套随机种子比较每像素前3/10/30次搜索：
- 前3次相对30次：1个像素(y=804,x=1525)漏掉更好解，速度变化6.354995 km/s，SSE降低12.9809%。
- 前10次相对30次：全区域最大速度差0.00047735 km/s，最大相对SSE改善3.50e-11。
- 这些是数值搜索稳定性检查，不是观测噪声下的置信区间，也不能证明全局最优。
