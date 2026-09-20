"""用背景吸收谷附近5点的局部抛物线顶点确定共同参考波长。"""
import numpy as np


def estimate_background_center(wavelength, I0, *, n_core_points=5):
    """返回线心及选中采样点；只拟合局部谷底，不拟合整条背景谱。

    默认选择最低点及左右各2点，局部拟合 I=a*(lambda-anchor)^2+
    b*(lambda-anchor)+c，lambda_ref=anchor-b/(2*a)。
    全部像素共用同一个参考值；云模型继续使用完整的实测I0。
    最低点靠边、平坦、重复波长或顶点异常时抛错，不回退到固定波长。
    """
    lam,bg = [np.asarray(np.ma.asarray(a,dtype=float).filled(np.nan)) for a in (wavelength,I0)]
    if lam.ndim != 1 or bg.shape != lam.shape:
        raise ValueError("wavelength 和 I0 必须为一维等长数组")
    if not isinstance(n_core_points,(int,np.integer)) or n_core_points < 3 or n_core_points%2 != 1:
        raise ValueError("n_core_points 必须为不小于3的奇数")
    indices = np.flatnonzero(np.isfinite(lam)&np.isfinite(bg))
    indices = indices[np.argsort(lam[indices])]
    x,y = lam[indices],bg[indices]
    if len(x) < n_core_points or np.any(np.diff(x)<=0):
        raise ValueError("背景线心拟合需要足够的不重复有效波长点")
    k,half = int(np.argmin(y)),n_core_points//2
    if k < half or k+half >= len(x):
        raise ValueError("背景最低点过于靠近波段边缘，请扩大波段或传入已标定lam0")
    core_indices = indices[k-half:k+half+1]
    cx,cy = lam[core_indices],bg[core_indices]
    if np.ptp(cy) <= 1e-10*max(float(np.max(np.abs(cy))),1.):
        raise ValueError("背景线心附近强度近似平坦，无法确定线心")
    anchor = float(x[k])
    coeff = np.polyfit(cx-anchor,cy,2)
    if coeff[0] <= 0:
        raise ValueError("背景线心局部二次拟合曲率非正")
    center = float(anchor-coeff[1]/(2*coeff[0]))
    if not np.isfinite(center) or not x[k-1] <= center <= x[k+1] or center <= 0:
        raise ValueError("拟合顶点偏离背景最低点，请检查I0或手动提供lam0")
    return dict(lam0=center,method="background_core_quadratic",n_core_points=int(n_core_points),
                core_indices=core_indices,core_wavelength=cx,core_intensity=cy,
                discrete_min_wavelength=anchor,polynomial_coefficients=coeff,
                core_rmse=float(np.sqrt(np.mean((np.polyval(coeff,cx-anchor)-cy)**2))))


def save_background_center(out,lam,bg,calibration):
    """标明实测背景与局部拟合，避免把数据连线误认成整体拟合。"""
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    np.savez_compressed(out/'background_center.npz',**calibration,wavelength=lam,I0=bg)
    fig=Figure(figsize=(9,4),layout='constrained')
    FigureCanvasAgg(fig)
    ax,zoom=fig.subplots(1,2)
    valid=np.isfinite(lam)&np.isfinite(bg)
    x,y=lam[valid],bg[valid]
    order=np.argsort(x)
    x,y=x[order],y[order]
    core=calibration['core_wavelength']
    for a in (ax,zoom):
        a.plot(x,y,'.-',ms=4,color='#0072B2',label='Measured I0 (samples joined)')
        a.axvline(calibration['lam0'],color='#D55E00',ls='--',label='Reference center')
        a.set_xlabel(r'Wavelength ($\mathrm{\AA}$)')
        a.ticklabel_format(axis='x',useOffset=False)
    ax.axvspan(core.min(),core.max(),color='#009E73',alpha=.12,label='Local fit window')
    fine=np.linspace(core.min(),core.max(),200)
    zoom.plot(fine,np.polyval(calibration['polynomial_coefficients'],fine-calibration['discrete_min_wavelength']),
              color='#009E73',lw=2,label='Local quadratic fit')
    zoom.plot(core,calibration['core_intensity'],'o',mfc='none',mec='#D55E00',ms=7,label='Samples used in fit')
    zoom.set_xlim(core.min()-.03,core.max()+.03)
    lo,hi=np.min(calibration['core_intensity']),np.max(calibration['core_intensity'])
    zoom.set_ylim(lo-.2*(hi-lo),hi+.2*(hi-lo))
    ax.set_title('Full measured background; no full-profile fit',fontsize=10)
    zoom.set_title(f"Only {calibration['n_core_points']} core samples fitted",fontsize=10)
    ax.set_ylabel('Intensity (DN)')
    ax.legend(fontsize=7)
    zoom.legend(fontsize=7)
    fig.suptitle(f"Background reference = {calibration['lam0']:.6f} Angstrom",fontsize=12)
    fig.savefig(out/'background_center.png',dpi=250)
    fig.savefig(out/'background_center.pdf')
    fig.clear()
