"""按 CHASE workshop 4.1 的单高斯+常数模型拟合背景参考线心。

https://ssdc.nju.edu.cn/web-service/software/Demo/CHASE_workshop_Python.html
示例使用原始波长索引50:90。若传入44:96切片，对应fit_slice=slice(6,46)。
"""
import warnings
import numpy as np
from scipy.optimize import curve_fit, OptimizeWarning

SOURCE_URL = "https://ssdc.nju.edu.cn/web-service/software/Demo/CHASE_workshop_Python.html"


def gaussian_background(wavelength, amplitude, center, width, constant):
    """A*exp[-((lambda-center)/width)^2]+C；吸收线A<0，width单位Angstrom。"""
    return amplitude*np.exp(-((np.asarray(wavelength)-center)/width)**2)+constant


def estimate_background_center(wavelength, I0, *, fit_slice=None):
    """返回高斯拟合线心、四个参数、协方差、拟合曲线与残差。

    与官方示例相同模型和curve_fit(method='lm')；为数值稳定先平移波长。
    用三个宽度初值检验局部优化，并选SSE最小的有效吸收线解。
    None使用全部输入波长点；slice只作用于背景线心拟合，云模型I0不变。
    covariance/center_formal_std是该模型下的形式误差，不含波长标定、
    背景不对称或拟合窗口选择造成的系统误差。
    """
    lam = np.asarray(np.ma.asarray(wavelength,dtype=float).filled(np.nan))
    bg = np.asarray(np.ma.asarray(I0,dtype=float).filled(np.nan))
    if lam.ndim != 1 or bg.shape != lam.shape:
        raise ValueError("wavelength 和 I0 必须为一维等长数组")
    if fit_slice is not None and not isinstance(fit_slice,slice):
        raise ValueError("fit_slice 必须为 slice 或 None")
    indices = np.arange(lam.size)[slice(None) if fit_slice is None else fit_slice]
    good = np.isfinite(lam[indices]) & np.isfinite(bg[indices])
    indices = indices[good]
    indices = indices[np.argsort(lam[indices])]
    x,y = lam[indices],bg[indices]
    if x.size <= 4 or np.any(np.diff(x) <= 0):
        raise ValueError("背景高斯拟合需要多于4个不重复的有效波长点")
    k = int(np.argmin(y))
    if k in (0,len(y)-1):
        raise ValueError("背景最低点位于拟合波段边缘，请调整background_fit_slice")
    if np.ptp(y) <= 1e-10*max(float(np.max(np.abs(y))),1.):
        raise ValueError("背景近似平坦，无法确定吸收线心")
    anchor = float(x[k])
    relative_x = x-anchor
    amplitude = float(y.min()-y.max())
    constant = float(y.max())
    candidates = []
    for width0 in (0.5,0.3,0.8):
        try:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always",OptimizeWarning)
                p,cov = curve_fit(gaussian_background,relative_x,y,
                                  p0=[amplitude,0.,width0,constant],method="lm",maxfev=10000)
            p[1] += anchor
            # 模型只依赖width平方，将负宽度等价映射为正，同时变换协方差。
            signs = np.array([1.,1.,1. if p[2] >= 0 else -1.,1.])
            cov = cov*signs[:,None]*signs[None,:]
            p[2] = abs(p[2])
            if not np.all(np.isfinite(p)) or p[0] >= 0 or p[2] <= 0 or not x[0] < p[1] < x[-1]:
                continue
            fitted = gaussian_background(x,*p)
            candidates.append(dict(params=p,covariance=cov,raw_sse=float(np.sum((fitted-y)**2)),
                                   initial_width=width0,warning="; ".join(str(w.message) for w in caught)))
        except (RuntimeError,ValueError,FloatingPointError):
            continue
    if not candidates:
        raise ValueError("背景高斯拟合未得到有效吸收线解；请检查拟合窗口和I0")
    best = min(candidates,key=lambda c:c["raw_sse"])
    p,cov = best["params"],best["covariance"]
    formal_std = np.sqrt(cov[1,1]) if np.isfinite(cov[1,1]) and cov[1,1]>=0 else np.nan
    return dict(lam0=float(p[1]),method="background_gaussian_constant",source_url=SOURCE_URL,
                params=p,param_names=np.array(["amplitude","center","width_angstrom","constant"]),
                covariance=cov,center_formal_std=float(formal_std),
                fit_indices=indices,fit_wavelength=x,fit_intensity=y,
                fit_model=gaussian_background(x,*p),fit_residual=y-gaussian_background(x,*p),
                model=gaussian_background(lam,*p),n_fit_points=int(x.size),
                raw_sse=best["raw_sse"],rmse=float(np.sqrt(best["raw_sse"]/x.size)),
                fit_range_angstrom=np.array([x[0],x[-1]]),
                candidate_params=np.array([c["params"] for c in candidates]),
                candidate_raw_sse=np.array([c["raw_sse"] for c in candidates]),
                warning=best["warning"])


def save_background_center(out, lam, bg, calibration):
    """保存完整诊断：背景谱/高斯拟合、放大线心、残差、参数与协方差。"""
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    np.savez_compressed(out/"background_center.npz",**calibration,wavelength=lam,I0=bg)
    fig = Figure(figsize=(9,5.6),layout="constrained")
    FigureCanvasAgg(fig)
    grid = fig.add_gridspec(2,2,height_ratios=[3,1])
    ax = fig.add_subplot(grid[0,0])
    zoom = fig.add_subplot(grid[0,1])
    res = fig.add_subplot(grid[1,:])
    good = np.isfinite(lam)&np.isfinite(bg)
    x,y = lam[good],bg[good]
    order = np.argsort(x)
    x,y = x[order],y[order]
    fx = calibration["fit_wavelength"]
    fine = np.linspace(fx[0],fx[-1],500)
    for a in (ax,zoom):
        a.plot(x,y,".",color="#0072B2",ms=4,label="Measured background I0")
        a.plot(fine,gaussian_background(fine,*calibration["params"]),color="#D55E00",label="Gaussian + constant")
        a.axvline(calibration["lam0"],color="0.3",ls="--",label="Fitted center")
        a.set_xlabel(r"Wavelength ($\mathrm{\AA}$)")
        a.ticklabel_format(axis="x",useOffset=False)
    ax.axvspan(fx[0],fx[-1],color="0.5",alpha=.06)
    ax.set_ylabel("Intensity (DN)")
    ax.legend(fontsize=7)
    zoom.set_xlim(calibration["lam0"]-.25,calibration["lam0"]+.25)
    near = np.abs(x-calibration["lam0"]) <= .3
    local = np.r_[y[near],gaussian_background(x[near],*calibration["params"])]
    if local.size:
        pad = max(float(np.ptp(local))*.15,1.)
        zoom.set_ylim(local.min()-pad,local.max()+pad)
    zoom.set_title("Line core",fontsize=10)
    res.plot(fx,calibration["fit_residual"],".-",color="#0072B2",ms=3)
    res.axhline(0,color="0.4",lw=.8)
    res.set_xlabel(r"Wavelength ($\mathrm{\AA}$)")
    res.set_ylabel("Obs - fit (DN)")
    fig.suptitle(f"Background Gaussian center = {calibration['lam0']:.6f} Angstrom\n"
                 f"N={calibration['n_fit_points']}; RMSE={calibration['rmse']:.3f} DN",fontsize=11)
    fig.savefig(out/"background_center.png",dpi=250)
    fig.savefig(out/"background_center.pdf")
    fig.clear()
