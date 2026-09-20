# 双云拟合的并行与加速

`double_cloud_pixel_fit.py` 保持八参数模型、层序、原始 SSE 目标和参数范围。默认每像素仍做30次独立全局搜索，每次 DE 默认最多1500代。更新不修改已经保存的拟合结果。

## 已有的混合方法 Notebook

保留现有 `double_mask`、波长范围、背景线心和双云结果保存循环。该循环已经使用 joblib 在像素之间并行，设置 `n_jobs = 16` 即可。不要用全矩形入口替代这个筛选循环。

当前正在运行的拟合先结束，再执行下面的代码，关闭之前缓存旧模块的 loky 工作进程并重新载入模块。无需清空 Notebook 变量。

```python
import importlib
import double_cloud_pixel_fit as dc
from joblib.externals.loky import get_reusable_executor

get_reusable_executor().shutdown(wait=True)
importlib.reload(dc)
n_jobs = 16
```

然后仅重新运行原来的双云拟合和后续保存、绘图单元格。`Parallel` 的上下文仍保持 `parallel_config(backend="loky", inner_max_num_threads=1)`。`fit_double_cloud_spectrum` 默认 `vectorized=True`，内部固定 `workers=1`，避免每个像素再启动16个进程。

## 整个矩形区域的拟合

`fit_double_cloud_cube(..., n_jobs=16)` 默认在像素之间并行；可以设置 `n_jobs=1` 串行，或 `n_jobs=-1` 使用可用逻辑处理器。进程数不超过像素数。`save_plots=False` 跳过逐像素图和总图，但仍保存数据。

## 加速实现与结果差异

- DE 一次用 NumPy 计算一批候选光谱，降低 Python 逐候选调用开销。
- 局部最小二乘使用八个参数的解析导数。
- 如果局部精炼反而增加 SSE，则保留精炼前解并标记局部优化未成功。
- 全矩形结果按任务完成顺序写CSV，参数图及逐像素文件仍按原图坐标对应。失败像素记录错误，其他像素继续。
- 向量化 DE 使用 `deferred` 更新，旧实现为 `immediate`。同样的随机种子可能得到不同搜索结果，不能宣称参数逐位不变或保证全局最优。`vectorized=False` 可以恢复逐候选、立即更新的 DE 路径，局部精炼仍使用解析导数。
- 单像素函数无需 joblib；全矩形并行入口需要支持 `generator_unordered` 的 joblib（本机1.6）。

验证命令：`python -m unittest test_double_cloud_pixel_fit -v`。测试覆盖解析导数、批量目标函数、无噪声谱线恢复、串行与并行坐标/结果一致性、失败像素及候选文件保存。速度对照结果保存在 `outputs/double_cloud_speed_comparison_20260915.json`，其结果仅代表注明的测试预算。
