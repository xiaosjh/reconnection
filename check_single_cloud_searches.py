"""按每像素随机种子顺序，比较前3/10/30次全局搜索的最佳结果。"""
from pathlib import Path
import csv
import json
import numpy as np


def compare_search_budgets(run_dir, budgets=(3,10,30)):
    """只读已有逐候选文件，不重新拟合；保存CSV/NPZ/JSON稳定性报告。

    所有预算使用同一组嵌套随机种子，因此SSE不能随预算增加而变差。
    搜索间的差异是数值稳定性诊断，不是观测误差或参数置信区间。
    """
    out = Path(run_dir)
    budgets = np.array(budgets,dtype=int)
    if np.any(budgets < 1) or np.any(np.diff(budgets) <= 0):
        raise ValueError("budgets必须为递增正整数")
    with np.load(out/'fit_maps.npz',allow_pickle=False) as f:
        x,y = f['x_pixels'],f['y_pixels']
    params = np.full((len(budgets),len(y),len(x),4),np.nan)
    sse = np.full(params.shape[:-1],np.nan)
    for j,yv in enumerate(y):
        for i,xv in enumerate(x):
            with np.load(out/'pixel_fits'/f'fit_y{yv:04d}_x{xv:04d}.npz',allow_pickle=False) as f:
                if 'candidate_x' not in f:
                    continue
                order = np.argsort(f['candidate_seed'])
                if order.size < budgets[-1]:
                    raise ValueError(f'像素y={yv},x={xv}只有{order.size}个候选，少于预算{budgets[-1]}')
                for k,n in enumerate(budgets):
                    subset = order[:n]
                    best = subset[np.argmin(f['candidate_raw_sse'][subset])]
                    p = f['candidate_x'][best].copy()
                    p[[0,2]] = np.exp(p[[0,2]])
                    params[k,j,i] = p
                    sse[k,j,i] = f['candidate_raw_sse'][best]
    valid = np.isfinite(sse).all(axis=0)
    reports = []
    for k,n in enumerate(budgets[:-1]):
        diff = np.abs(params[k]-params[-1])
        relative_gain = np.divide(sse[k]-sse[-1],sse[k],out=np.zeros_like(sse[k]),where=sse[k]>0)
        reports.append(dict(
            first_runs=int(n),reference_runs=int(budgets[-1]),n_pixels=int(valid.sum()),
            max_abs_delta_v_km_s=float(np.max(diff[...,1][valid])) if valid.any() else None,
            pixels_delta_v_gt_0p1=int(np.sum((diff[...,1]>.1)&valid)),
            max_abs_delta_params_tau_v_W_S=np.max(diff[valid],axis=0).tolist() if valid.any() else None,
            max_relative_sse_improvement=float(np.max(relative_gain[valid])) if valid.any() else None,
            pixels_relative_sse_improvement_gt_1e_6=int(np.sum((relative_gain>1e-6)&valid))))
    report = dict(run_dir=str(out),budgets=budgets.tolist(),comparisons=reports,
                  interpretation='Nested numerical search budgets; not parameter uncertainty or proof of global optimum')
    (out/'search_stability.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    np.savez_compressed(out/'search_stability.npz',budgets=budgets,best_params=params,best_raw_sse=sse,x_pixels=x,y_pixels=y)
    with (out/'search_stability.csv').open('w',newline='',encoding='utf-8-sig') as f:
        writer=csv.writer(f)
        writer.writerow(['y','x','first_runs','tau0','v','W','S','raw_sse'])
        for j,yv in enumerate(y):
            for i,xv in enumerate(x):
                for k,n in enumerate(budgets):
                    writer.writerow([yv,xv,n,*params[k,j,i],sse[k,j,i]])
    return report
