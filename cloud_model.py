import sunpy
import sunpy.map
import numpy as np
from math import *
import astropy.units as u
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sunpy.coordinates import frames
import matplotlib.gridspec as gridspec
from scipy.interpolate import interp1d
from astropy.coordinates import SkyCoord
from matplotlib.patches import ConnectionPatch
import glob
from scipy.io import readsav
from scipy.ndimage import zoom
from scipy.special import wofz
from functools import partial
from astropy import constants as const


def voigt(lambda1, I_voigt, gamma_voigt, sigma_voigt):
    """
    Voigt 分量函数
    """
    z = (lambda1 + 1j * gamma_voigt) / (sigma_voigt * np.sqrt(2))
    voigt_component = np.real(wofz(z)) / (sigma_voigt * np.sqrt(2 * np.pi)) * I_voigt
    return voigt_component


def voigt_max_coefficient(gamma_voigt, sigma_voigt):
    """
    计算 Voigt 的最大系数，用于归一化
    """
    z_max = (1j * gamma_voigt) / (sigma_voigt * np.sqrt(2))
    return (sigma_voigt * np.sqrt(2 * np.pi)) / (np.real(wofz(z_max)))

# 单高斯 + Voigt + 常数
def voigt_single_cloud(lambda1, I_voigt, gamma_voigt, sigma_voigt,
                       I_cloud, mu_cloud, sigma_cloud, I_constant, S_l, I_0):

    # z = (lambda1 + 1j * gamma_voigt) / (sigma_voigt * np.sqrt(2))
    # voigt_component = np.real(wofz(z)) / (sigma_voigt * np.sqrt(2 * np.pi)) * I_voigt
    voigt_component = voigt(lambda1, I_voigt, gamma_voigt, sigma_voigt)
    tau1 = I_cloud * np.exp(-0.5*((lambda1 - mu_cloud) / sigma_cloud) ** 2)
    cloud_component = (S_l - I_0) * (1 - np.exp(-tau1))
    return voigt_component + cloud_component + I_constant

#数据准备
#读取数据
dir=r'D:\Learning\PHD1st\magnetic_reconnecion\data\CHASE_Ha\RSM20240618T211711_0000_HA.fits'
rsm=fits.open(dir)
spectrum=rsm[1].data[:,795,1521]
lam = rsm[1].header['CRVAL3'] + np.arange(rsm[1].header['NAXIS3']) * rsm[1].header['CDELT3']

#计算I_0
wavelength_point_list_ha_t0=lam
F_t0_l_TR=rsm[1].data[:,795,1521]
F_t0_l_TR_interpolated = interp1d(wavelength_point_list_ha_t0, F_t0_l_TR,
                                kind='cubic', fill_value='extrapolate')(
                                wavelength_point_list_ha_t0)
# TR_width = 70
# TR_height = 70
# psr = 1.0436
# r_fe=rsm[1].header['R_SUN']
# f_t0_lcont_TR=float(interp1d(wavelength_point_list_ha_t0, F_t0_l_TR,
#                               kind='cubic', fill_value='extrapolate')(6568))
# #这个I0是一小片区域相对于太阳整体的相对强度，我有必要用这种相对强度么
# I_0 = F_t0_l_TR_interpolated * (2 * TR_width / psr + 1) * (2 * TR_height / psr + 1) / (
#                 f_t0_lcont_TR * (2 * r_fe + 1) ** 2)
I_0=rsm[1].data[:,775:815,1501:1541].mean(axis=(1,2))
#待拟合数据
ha_left_index=50
ha_right_index=100
x_fit = lam[ha_left_index:ha_right_index + 1] - 6562.8 
y_fit = spectrum[ha_left_index:ha_right_index + 1]

#画图
def cloud_component(I_cloud,lambda1,mu_cloud,sigma_cloud,S_l,I_0):
    tau1 = I_cloud * np.exp(-0.5*((lambda1 - mu_cloud) / sigma_cloud) ** 2)
    cloud_component = (S_l - I_0) * (1 - np.exp(-tau1))
    return cloud_component

def fit_single_all(paras,x_fit,y_fit,I_0):
    voigt_single_cloud1 = partial(voigt_single_cloud, I_0=I_0[ha_left_index:ha_right_index + 1])
    popt,pcov=curve_fit(
        voigt_single_cloud1,x_fit,y_fit,
        p0=paras[0::3],
        bounds=(paras[1::3],paras[2::3])
    )
    c=const.c
    v1=c.value*popt[4]/(6562.8*1000)
    print(f"Velocity component 1 = {v1:.3f} km/s")
    I_voigt_fit=popt[0]
    I_cloud1_fit=popt[3]
    mu_cloud1_fit=popt[4]
    sigma_cloud1_fit=popt[5]
    cloud_component1=cloud_component(I_cloud1_fit,x_fit,mu_cloud1_fit,
                                 sigma_cloud1_fit,S_l_fit,
                                 I_0[ha_left_index:ha_right_index + 1])
    voigt1=voigt(x_fit,I_voigt_fit,gamma_voigt_fit,sigma_voigt_fit)

    plt.plot(x_fit,voigt1,label='Voigt profile')
    plt.plot(x_fit,cloud_component1,label='Cloud component profile')
    plt.plot(x_fit,voigt1+cloud_component1+I_constant_fit,label='voigt+cloud+I_constant')
    plt.plot(x_fit,y_fit,label='Total profile')
    plt.legend()
    plt.show()
    print(popt)
    return 