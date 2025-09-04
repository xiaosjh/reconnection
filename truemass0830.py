# ==============================
# Input - begin

input_dir_ha = r"/Volumes/T7Shield/CHASEdata/CME1/20240830/data00/Ha"
input_dir_fe = r"/Volumes/T7Shield/CHASEdata/CME1/20240830/data00/Fe"
output_dir_figures = r"/Users/yvonne/Documents/sunas/CME1/output0830/Figures"

image_half_side = 1100  # The half side length of the cropped image. Unit: pixels.
psr = 1.0436  # Pixel spatial resolution. Unit: arcsec per pixel.

TR_center = [-460,
             220]  # The center of the target region. Relative to the center of the solar disk. [relative_x , relative_y]. Unit: arcsec.
TR_half_side = 70  # The half side length of the target region. Unit: arcsec.
TR_height = 70
TR_width = 70

ha_left_index = 10  # The index of the left  border of Ha line. Minimum: 30
ha_right_index = 117  # The index of the right border of Ha line. Maximum: 117

fit_starttime_index = 1  # The index of the start time for fitting. Minimum: 0
fit_endtime_index = 25  # The index of the end   time for fitting. Maximum: 22

dt_format = '%Y-%m-%d %H:%M:%S'

# Input - end
# ==============================

# ==============================
# import - begin

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os
import time
from scipy.interpolate import interp1d
from scipy.special import wofz
from scipy.optimize import curve_fit
import gc
from datetime import datetime
from matplotlib.font_manager import FontProperties
from functools import partial
import scipy.integrate as integrate

# import - end
# ==============================

# ============================================================== #
# The following codes cannot be modified, unless there are bugs. #
# ============================================================== #

# ==============================
# Functions - begin

# 单高斯 + Voigt + 常数
def voigt_single_cloud(lambda1, I_voigt, gamma_voigt, sigma_voigt,
                       I_cloud, mu_cloud, sigma_cloud, I_constant, S_l, I_0):

    z = (lambda1 + 1j * gamma_voigt) / (sigma_voigt * np.sqrt(2))
    voigt_component = np.real(wofz(z)) / (sigma_voigt * np.sqrt(2 * np.pi)) * I_voigt
    tau1 = I_cloud * np.exp(-0.5*((lambda1 - mu_cloud) / sigma_cloud) ** 2)
    cloud_component = (S_l - I_0) * (1 - np.exp(-tau1))
    return voigt_component + cloud_component + I_constant
# 双高斯 + Voigt + 常数
def voigt_double_cloud(lambda1, I_voigt, gamma_voigt, sigma_voigt,
                       I_cloud1, mu_cloud1, sigma_cloud1,
                       I_cloud2, mu_cloud2, sigma_cloud2, I_constant):
    z = (lambda1 + 1j * gamma_voigt) / (sigma_voigt * np.sqrt(2))
    voigt_component = np.real(wofz(z)) / (sigma_voigt * np.sqrt(2 * np.pi)) * I_voigt
    tau1 = I_cloud1 * np.exp(-((lambda1 - mu_cloud1) / sigma_cloud1) ** 2)
    cloud_component1 = S_l * np.exp(-tau1)
    tau2 = I_cloud2 * np.exp(-((lambda1 - mu_cloud2) / sigma_cloud2) ** 2)
    cloud_component2 = S_l * np.exp(-tau2)
    return voigt_component + cloud_component1 + cloud_component2 + I_constant


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



line_center = 6562.8
line_width = 6.0
def calculate_ew(spectrum, wavelength, line_center, line_width):
    line_region = (wavelength >= line_center - line_width / 2) & (wavelength <= line_center + line_width / 2)
    intensity_line = spectrum[line_region]
    wavelength_line = wavelength[line_region]
    ew = np.trapz(intensity_line, wavelength_line)
    return ew

def tau_lambda(lambda_vals, lambda_0, tau_0, v_shift, W):
    return tau_0 * np.exp(-0.5 * ((lambda_vals / lambda_0 - (1 + v_shift / c)) / (W / c))**2)

def intensity_lambda(tau_vals, I0_lambda, S):
    return (S - I0_lambda) * (1 - np.exp(-tau_vals))

def compute_EW(lambda_vals, I_lambda):
    EW = integrate.simps(I_lambda, lambda_vals)
    return EW
# Functions - end
# ==============================

print('==============================')
try:
    os.mkdir(output_dir_figures + '/image_hac')  # hac: H alpha line
    print('Directory {} has been created. '.format(output_dir_figures + '/image_hac'))
except:
    print('Directory {} already exists.   '.format(output_dir_figures + '/image_hac'))

try:
    os.mkdir(output_dir_figures + '/image_few')  # few: Fe I line wing
    print('Directory {} has been created. '.format(output_dir_figures + '/image_few'))
except:
    print('Directory {} already exists.   '.format(output_dir_figures + '/image_few'))

try:
    os.mkdir(output_dir_figures + '/ds')  # ds: dynamic spectra
    print('Directory {} has been created. '.format(output_dir_figures + '/ds'))
except:
    print('Directory {} already exists.   '.format(output_dir_figures + '/ds'))

try:
    os.mkdir(output_dir_figures + '/ds_snapshots')
    print('Directory {} has been created. '.format(output_dir_figures + '/ds_snapshots'))
except:
    print('Directory {} already exists.   '.format(output_dir_figures + '/ds_snapshots'))

try:
    os.mkdir(output_dir_figures + '/ds_TReqFD')  # ds: dynamic spectra
    print('Directory {} has been created. '.format(output_dir_figures + '/ds_TReqFD'))
except:
    print('Directory {} already exists.   '.format(output_dir_figures + '/ds_TReqFD'))

try:
    os.mkdir(output_dir_figures + '/ds_snapshots_TReqFD')
    print('Directory {} has been created. '.format(output_dir_figures + '/ds_snapshots_TReqFD'))
except:
    print('Directory {} already exists.   '.format(output_dir_figures + '/ds_snapshots_TReqFD'))

try:
    os.mkdir(output_dir_figures + '/EW')
    print('Directory {} has been created. '.format(output_dir_figures + '/EW'))
except:
    print('Directory {} already exists.   '.format(output_dir_figures + '/EW'))

# try:
#    os.mkdir(output_dir_figures+'/ew')
#    print('Directory {} has been created. '.format(output_dir_figures+'/ew'))
# except:
#    print('Directory {} already exists.   '.format(output_dir_figures+'/ew'))
print('==============================')
TR_center_x = TR_center[0]
TR_center_y = TR_center[1]
input_filename_ha_list = os.listdir(input_dir_ha)
input_filename_fe_list = os.listdir(input_dir_fe)
input_fitsfilename_ha_list = []
input_fitsfilename_fe_list = []
for i in range(len(input_filename_ha_list)):
    if input_filename_ha_list[i][0:3] == 'RSM' and input_filename_ha_list[i][-7:] == 'HA.fits':
        input_fitsfilename_ha_list.append(input_filename_ha_list[i])
input_fitsfilename_ha_list = sorted(input_fitsfilename_ha_list)
for i in range(len(input_filename_fe_list)):
    if input_filename_fe_list[i][0:3] == 'RSM' and input_filename_fe_list[i][-7:] == 'FE.fits':
        input_fitsfilename_fe_list.append(input_filename_fe_list[i])
input_fitsfilename_fe_list = sorted(input_fitsfilename_fe_list)
print('==============================')
print('Total number of Ha FITS files:', len(input_fitsfilename_ha_list))
print('Total number of Fe FITS files:', len(input_fitsfilename_fe_list))
print('==============================')

start_time_str_list = []
EW_t1_list = []
EW_t1_TReqFD_list = []
mu_cloud1_list = []
sigma_cloud1_list = []
I_cloud1_list = []
I_cloud2_list = []
mu_cloud2_list = []
sigma_cloud2_list = []
I_voigt_list = []
gamma_voigt_list = []
sigma_voigt_list = []
S_l = []
MT = [0]
ew_list = []
MF = []
MT10 = []
for i in range(len(input_fitsfilename_ha_list)):
    print('------------------------------')
    print('Progress:', str(i + 1) + ' / ' + str(len(input_fitsfilename_ha_list)))
    input_fitsfilename_ha = input_fitsfilename_ha_list[i]
    print('Current Ha FITS filename:', input_fitsfilename_ha)
    input_fitsfilename_fe = input_fitsfilename_fe_list[i]
    print('Current Fe FITS filename:', input_fitsfilename_fe)

    os.chdir(input_dir_ha)
    input_fitsfile_ha = fits.open(input_fitsfilename_ha)
    input_fitsfile_ha_0 = input_fitsfile_ha[0]
    input_fitsfile_ha_1 = input_fitsfile_ha[1]

    os.chdir(input_dir_fe)
    input_fitsfile_fe = fits.open(input_fitsfilename_fe)
    input_fitsfile_fe_0 = input_fitsfile_fe[0]
    input_fitsfile_fe_1 = input_fitsfile_fe[1]

    cx_ha = float(input_fitsfile_ha_1.header['CRPIX1'])  # Unit: pixels.
    cy_ha = float(input_fitsfile_ha_1.header['CRPIX2'])  # Unit: pixels.
    r_ha = float(input_fitsfile_ha_1.header['R_SUN'])  # Unit: pixels.

    cx_fe = float(input_fitsfile_fe_1.header['CRPIX1'])  # Unit: pixels.
    cy_fe = float(input_fitsfile_fe_1.header['CRPIX2'])  # Unit: pixels.
    r_fe = float(input_fitsfile_fe_1.header['R_SUN'])  # Unit: pixels.

    start_time_year_int = int(input_fitsfilename_ha[3:7])
    start_time_month_int = int(input_fitsfilename_ha[7:9])
    start_time_day_int = int(input_fitsfilename_ha[9:11])
    start_time_hour_int = int(input_fitsfilename_ha[12:14])
    start_time_minute_int = int(input_fitsfilename_ha[14:16])
    start_time_second_int = int(input_fitsfilename_ha[16:18])
    start_time_str = str(start_time_year_int).zfill(4) + '-' + str(start_time_month_int).zfill(2) + '-' + str(
        start_time_day_int).zfill(2) + ' ' + str(start_time_hour_int).zfill(2) + ':' + str(start_time_minute_int).zfill(
        2) + ':' + str(start_time_second_int).zfill(2)
    start_time_str_list.append(start_time_str)
    print('The start time of scanning:', start_time_str)
    if i == 1:
        start_time_str_0 = start_time_str
    else:
        start_time_str_f = start_time_str  # f: final

    end_time_str = input_fitsfile_ha_1.header['END_TIME'][:10] + ' ' + input_fitsfile_ha_1.header['END_TIME'][11:]
    print('The end time of scanning:  ', end_time_str)

    wavelength_point_list_ha_t1 = np.array(
        [input_fitsfile_ha_1.header['CRVAL3'] + i * input_fitsfile_ha_1.header['CDELT3'] for i in range(118)])
    wavelength_point_list_fe_t1 = np.array(
        [input_fitsfile_fe_1.header['CRVAL3'] + i * input_fitsfile_fe_1.header['CDELT3'] for i in range(46)])
    wavelength_interval_ha_t1 = input_fitsfile_ha_1.header['CDELT3']
    if i == 0:
        wavelength_point_list_ha_t0 = wavelength_point_list_ha_t1
        wavelength_point_list_fe_t0 = wavelength_point_list_fe_t1
        wavelength_interval_ha_t0 = wavelength_interval_ha_t1

    # ------------------------------
    # Image of Ha line center (full disk) - begin
    data_hac_image = input_fitsfile_ha_1.data[68,
                     int(round(cy_ha - image_half_side, 0)):int(round(cy_ha + image_half_side, 0)) + 1,
                     int(round(cx_ha - image_half_side, 0)):int(
                         round(cx_ha + image_half_side)) + 1]  # hac: Ha line center

    os.chdir(output_dir_figures + '/image_hac')

    plt.figure(figsize=(6.8897638, 3.3464567), dpi=450)
    plt.title(start_time_str + ' UT', fontsize=6, fontname='Arial')
    plt.axis('off')
    plt.subplot(1, 2, 1)
    plt.xticks([-1000, -750, -500, -250, 0, 250, 500, 750, 1000], fontsize=5, fontname='Arial')
    plt.yticks([-1000, -750, -500, -250, 0, 250, 500, 750, 1000], fontsize=5, fontname='Arial')
    plt.minorticks_on()
    plt.xlabel('Solar X (arcsec) from center', fontsize=6, fontname='Arial')
    plt.ylabel('Solar Y (arcsec) from center', fontsize=6, fontname='Arial')
    plt.imshow(data_hac_image, origin='lower', cmap='afmhot', vmin=0, vmax=4 * data_hac_image.mean(), extent=(
    (-image_half_side - 0.5) * psr, (image_half_side + 0.5) * psr, (-image_half_side - 0.5) * psr,
    (image_half_side + 0.5) * psr))

    # Target region
    plt.plot([TR_center_x - TR_width, TR_center_x + TR_width],
             [TR_center_y + TR_height, TR_center_y + TR_height], linestyle='-', color='cyan', linewidth=1)  # top
    plt.plot([TR_center_x - TR_width, TR_center_x + TR_width],
             [TR_center_y -TR_height, TR_center_y - TR_height], linestyle='-', color='cyan',
             linewidth=1)  # bottom
    plt.plot([TR_center_x - TR_width, TR_center_x - TR_width],
             [TR_center_y - TR_height, TR_center_y + TR_height], linestyle='-', color='cyan', linewidth=1)  # left
    plt.plot([TR_center_x + TR_width, TR_center_x + TR_width],
             [TR_center_y - TR_height, TR_center_y + TR_height], linestyle='-', color='cyan',
             linewidth=1)  # right

    # Frame of full solar disk
    # plt.plot(r_ha*psr * np.cos(np.linspace(0 , 2*np.pi , 512)) , r_ha*psr * np.sin(np.linspace(0 , 2*np.pi , 512)) , linestyle='--' , color='white' , linewidth=1)  # circle
    plt.plot([-r_ha * psr, r_ha * psr], [r_ha * psr, r_ha * psr], linestyle='--', color='white', linewidth=1)  # top
    plt.plot([-r_ha * psr, r_ha * psr], [-r_ha * psr, -r_ha * psr], linestyle='--', color='white',
             linewidth=1)  # bottom
    plt.plot([-r_ha * psr, -r_ha * psr], [-r_ha * psr, r_ha * psr], linestyle='--', color='white', linewidth=1)  # left
    plt.plot([r_ha * psr, r_ha * psr], [-r_ha * psr, r_ha * psr], linestyle='--', color='white', linewidth=1)  # right

    # Image of Ha line center (full disk) - end
    # ------------------------------

    # ------------------------------
    # Image of Ha line center (TR) - begin
    data_ha_TR = input_fitsfile_ha_1.data[:, int(round(cy_ha + TR_center_y / psr - TR_height / psr)):int(
        round(cy_ha + TR_center_y / psr + TR_height / psr)) + 1,
                 int(round(cx_ha + TR_center_x / psr - TR_width / psr)):int(
                     round(cx_ha + TR_center_x / psr + TR_width / psr)) + 1]
    data_hac_TR_image = data_ha_TR[68, :, :]
    data_ha_FD = input_fitsfile_ha_1.data[:, int(round(cy_ha - r_ha)):int(round(cy_ha + r_ha)) + 1,
                 int(round(cx_ha - r_ha)):int(round(cx_ha + r_ha)) + 1]

    plt.subplot(1, 2, 2)
    plt.xticks([-100, -50, 0, 50, 100, 150, 200], fontsize=5, fontname='Arial')
    plt.yticks([-500, -450, -400, -350, -300, -250, -200], fontsize=5, fontname='Arial')
    plt.minorticks_on()
    plt.xlabel('Solar X (arcsec) from center', fontsize=6, fontname='Arial')
    plt.ylabel('Solar Y (arcsec) from center', fontsize=6, fontname='Arial')
    plt.imshow(data_hac_TR_image, origin='lower', cmap='afmhot', vmin=0, vmax=4 * data_hac_image.mean(), extent=(
    TR_center_x - TR_width, TR_center_x + TR_width, TR_center_y - TR_height, TR_center_y + TR_height))

    plt.subplots_adjust(left=0.10, right=0.95, bottom=0.10, top=0.90, hspace=0.20)

    plt.savefig(input_fitsfilename_ha[:-5] + '_image_hac_TR.png')
    plt.close()

    # Image of Ha line center (TR) - end
    # ------------------------------

    # ------------------------------
    # Image of Fe line wing (full disk) - begin
    data_few_image = input_fitsfile_fe_1.data[5,
                     int(round(cy_fe - image_half_side, 0)):int(round(cy_fe + image_half_side, 0)) + 1,
                     int(round(cx_fe - image_half_side, 0)):int(
                         round(cx_fe + image_half_side)) + 1]  # few: Fe I line wing

    os.chdir(output_dir_figures + '/image_few')

    plt.figure(figsize=(6.8897638, 3.3464567), dpi=450)
    plt.title(start_time_str + ' UT', fontsize=6, fontname='Arial')
    plt.axis('off')
    plt.subplot(1, 2, 1)
    plt.xticks([-1000, -750, -500, -250, 0, 250, 500, 750, 1000], fontsize=5, fontname='Arial')
    plt.yticks([-1000, -750, -500, -250, 0, 250, 500, 750, 1000], fontsize=5, fontname='Arial')
    plt.minorticks_on()
    plt.xlabel('Solar X (arcsec) from center', fontsize=6, fontname='Arial')
    plt.ylabel('Solar Y (arcsec) from center', fontsize=6, fontname='Arial')
    plt.imshow(data_few_image, origin='lower', cmap='afmhot', vmin=0, vmax=4 * data_few_image.mean(), extent=(
    (-image_half_side - 0.5) * psr, (image_half_side + 0.5) * psr, (-image_half_side - 0.5) * psr,
    (image_half_side + 0.5) * psr))

    plt.plot([TR_center_x - TR_width, TR_center_x + TR_width],
             [TR_center_y + TR_height, TR_center_y + TR_height], linestyle='-', color='cyan', linewidth=1)  # top
    plt.plot([TR_center_x - TR_width, TR_center_x + TR_width],
             [TR_center_y -TR_height, TR_center_y - TR_height], linestyle='-', color='cyan',
             linewidth=1)  # bottom
    plt.plot([TR_center_x - TR_width, TR_center_x - TR_width],
             [TR_center_y - TR_height, TR_center_y + TR_height], linestyle='-', color='cyan', linewidth=1)  # left
    plt.plot([TR_center_x + TR_width, TR_center_x + TR_width],
             [TR_center_y - TR_height, TR_center_y + TR_height], linestyle='-', color='cyan',
             linewidth=1)  # right

    # Frame of full solar disk
    # plt.plot(r_fe*psr * np.cos(np.linspace(0 , 2*np.pi , 512)) , r_fe*psr * np.sin(np.linspace(0 , 2*np.pi , 512)) , linestyle='--' , color='white' , linewidth=1)  # circle
    plt.plot([-r_fe * psr, r_fe * psr], [r_fe * psr, r_fe * psr], linestyle='--', color='white', linewidth=1)  # top
    plt.plot([-r_fe * psr, r_fe * psr], [-r_fe * psr, -r_fe * psr], linestyle='--', color='white',
             linewidth=1)  # bottom
    plt.plot([-r_fe * psr, -r_fe * psr], [-r_fe * psr, r_fe * psr], linestyle='--', color='white', linewidth=1)  # left
    plt.plot([r_fe * psr, r_fe * psr], [-r_fe * psr, r_fe * psr], linestyle='--', color='white', linewidth=1)  # right

    # Image of Fe line wing (full disk) - end
    # ------------------------------

    # ------------------------------
    # Image of Fe line wing (TR) - begin
    data_fe_TR = input_fitsfile_fe_1.data[:, int(round(cy_fe + TR_center_y / psr - TR_height / psr)):int(
        round(cy_fe + TR_center_y / psr + TR_height / psr)) + 1,
                 int(round(cx_fe + TR_center_x / psr - TR_width / psr)):int(
                     round(cx_fe + TR_center_x / psr + TR_width / psr)) + 1]
    data_few_TR_image = data_fe_TR[5, :, :]
    data_fe_TR_ave = data_fe_TR.mean(axis=(1, 2))
    data_few_TR_ave = float(
        interp1d(wavelength_point_list_fe_t1, data_fe_TR_ave, kind='cubic', fill_value='extrapolate')(6568))

    data_fe_FD = input_fitsfile_fe_1.data[:, int(round(cy_fe - r_fe)):int(round(cy_fe + r_fe)) + 1,
                 int(round(cx_fe - r_fe)):int(round(cx_fe + r_fe)) + 1]
    data_fe_FD_ave = data_fe_FD.mean(axis=(1, 2))
    data_few_FD_ave = float(
        interp1d(wavelength_point_list_fe_t1, data_fe_FD_ave, kind='cubic', fill_value='extrapolate')(6568))

    plt.subplot(1, 2, 2)
    plt.xticks([-100, -50, 0, 50, 100, 150, 200], fontsize=5, fontname='Arial')
    plt.yticks([-500, -450, -400, -350, -300, -250, -200], fontsize=5, fontname='Arial')
    plt.minorticks_on()
    plt.xlabel('Solar X (arcsec) from center', fontsize=6, fontname='Arial')
    plt.ylabel('Solar Y (arcsec) from center', fontsize=6, fontname='Arial')
    plt.imshow(data_few_TR_image, origin='lower', cmap='afmhot', vmin=0, vmax=4 * data_few_image.mean(), extent=(
    TR_center_x - TR_width, TR_center_x + TR_width, TR_center_y - TR_height, TR_center_y + TR_height))

    plt.subplots_adjust(left=0.10, right=0.95, bottom=0.10, top=0.90, hspace=0.20)

    plt.savefig(input_fitsfilename_fe[:-5] + '_image_few_TR.png')
    plt.close()

    # Image of Fe line wing (TR) - end
    # ------------------------------
    data_hac_TR_image1 = np.array(data_ha_TR[68, :, :], dtype=np.float64)  # 复制第 68 层数据
    for n in range(1, 21):  # 从 1 到 10 进行累加
        data_hac_TR_image1 += data_ha_TR[68 - n, :, :]
        data_hac_TR_image1 += data_ha_TR[68 + n, :, :]
    data_hac_TR_image1 /= 41

    low_flux_indices = np.where(data_hac_TR_image1 < 500)
    i_indices, j_indices = low_flux_indices
    high_flux_indices = np.where(data_hac_TR_image1 >= 500)
    k_indices, l_indices = high_flux_indices
    indices = np.where(input_fitsfile_ha_1.data > 0)

    N_y, N_x = data_hac_TR_image1.shape
    x_min = TR_center_x - TR_width
    x_max = TR_center_x + TR_width
    y_min = TR_center_y - TR_height
    y_max = TR_center_y + TR_height
    x_coords = x_min + (j_indices / N_x) * (x_max - x_min)
    y_coords = y_min + (i_indices / N_y) * (y_max - y_min)

    # ------------------------------
    # Calculation of SaaS dynamic spectra (inside loop) - begin
    if i == 0:
        f_t1_l_TR = data_ha_TR[:, k_indices, l_indices].mean(axis=1) # t1: current time
    else:
        f_t1_l_TR = data_ha_TR.mean(axis=(1, 2))

    f_t1_lcont_TR = data_few_TR_ave
    if i == 0:
        f_t0_lcont_TR = f_t1_lcont_TR  # t0: initial time
    F_t1_l_TR = f_t1_l_TR * (f_t0_lcont_TR / f_t1_lcont_TR)
    if i == 0:
        F_t0_l_TR = F_t1_l_TR

    # ..............................

    f_t1_l_FD = data_ha_FD.mean(axis=(1, 2))

    f_t1_lcont_FD = data_few_FD_ave  # FD: full disk
    if i == 0:
        f_t0_lcont_FD = f_t1_lcont_FD

    F_t1_l_FD = f_t1_l_FD * (f_t0_lcont_FD / f_t1_lcont_FD)
    if i == 0:
        F_t0_l_FD = F_t1_l_FD

    # ..............................
    F_t0_l_TR_interpolated = interp1d(wavelength_point_list_ha_t0, F_t0_l_TR, kind='cubic', fill_value='extrapolate')(
        wavelength_point_list_ha_t1)
    F_t0_l_FD_interpolated = interp1d(wavelength_point_list_ha_t0, F_t0_l_FD, kind='cubic', fill_value='extrapolate')(
        wavelength_point_list_ha_t1)

    Delta_S_t1_l = ((F_t1_l_TR - F_t0_l_TR_interpolated) * (2 * TR_width / psr + 1) *(2 * TR_height / psr + 1)) / (
                f_t0_lcont_TR * (2 * r_fe + 1) ** 2)

    k = (2 * TR_width / psr + 1) *(2 * TR_height / psr + 1) / (
                f_t0_lcont_TR * (2 * r_fe + 1) ** 2)
    #Delta_S = F_t1_l_TR - F_t0_l_TR_interpolated
    Delta_S_t1_l_TReqFD = (F_t1_l_FD - F_t0_l_FD_interpolated) / f_t0_lcont_TR
    #Delta_S_interpolated = interp1d(wavelength_point_list_ha_t1, Delta_S, kind='cubic',
                                         #fill_value='extrapolate')(wavelength_point_list_ha_t0)

    Delta_S_t1_l_interpolated = interp1d(wavelength_point_list_ha_t1, Delta_S_t1_l, kind='cubic',
                                         fill_value='extrapolate')(wavelength_point_list_ha_t0)
    Delta_S_t1_l_interpolated_TReqFD = interp1d(wavelength_point_list_ha_t1, Delta_S_t1_l_TReqFD, kind='cubic',
                                                fill_value='extrapolate')(wavelength_point_list_ha_t0)
    EW_delta_S_t1 = interp1d(wavelength_point_list_ha_t1, ((F_t1_l_TR - F_t0_l_TR_interpolated)/f_t0_lcont_TR)* (2 * TR_width / psr + 1) *(2 * TR_height / psr + 1) / (
                (2 * r_fe + 1) ** 2), kind='cubic',
                                         fill_value='extrapolate')(wavelength_point_list_ha_t0)
    ew = calculate_ew(EW_delta_S_t1, wavelength_point_list_ha_t0, line_center, line_width)
    ew_list.append(ew)
    delta_ew_list = np.array(ew_list) - ew_list[0]

    if i == 0:
        Delta_S_t_l_interpolated = Delta_S_t1_l_interpolated
        Delta_S_t_l_interpolated_TReqFD = Delta_S_t1_l_interpolated_TReqFD
        I_0 = F_t0_l_TR_interpolated * (2 * TR_width / psr + 1) * (2 * TR_height / psr + 1) / (
                f_t0_lcont_TR * (2 * r_fe + 1) ** 2)
    elif i >= 1:
        Delta_S_t_l_interpolated = np.vstack((Delta_S_t_l_interpolated, Delta_S_t1_l_interpolated))
        Delta_S_t_l_interpolated_TReqFD = np.vstack((Delta_S_t_l_interpolated_TReqFD, Delta_S_t1_l_interpolated_TReqFD))

    def cloud(x, tau_0, lam_0, delta_lam, S_l):
        """
        云模型公式：
        I(λ) = S_l * exp(-tau(λ))
        tau(λ) = tau_0 * exp(-((λ - λ0) / delta_λ)^2)
        """
        tau = tau_0 * np.exp(-0.5 * ((x - lam_0) / delta_lam) ** 2)
        intensity = (S_l - I_0[ha_left_index:ha_right_index + 1]) * (1 - np.exp(-tau))

        return intensity
    # Calculation of SaaS dynamic spectra (inside loop) - end
    # ------------------------------

    # ------------------------------
    # Fit Delta_S_t1_l_interpolated - begin
    if i >= fit_starttime_index and i <= fit_endtime_index:
        x_fit = np.array(wavelength_point_list_ha_t0)[ha_left_index:ha_right_index + 1] - 6562.8
        y_fit = Delta_S_t1_l_interpolated[ha_left_index:ha_right_index + 1]

        # Voigt 分量参数
        gamma_voigt_initial = 0.5
        gamma_voigt_min = 0.0
        gamma_voigt_max = 7

        sigma_voigt_initial = 0.5
        sigma_voigt_min = 0.0
        sigma_voigt_max = 7

        I_voigt_initial = 1.0 * np.max(y_fit) * voigt_max_coefficient(gamma_voigt_initial, sigma_voigt_initial)
        I_voigt_min = 0.0
        I_voigt_max = 5.0 * np.max(y_fit) * voigt_max_coefficient(gamma_voigt_initial, sigma_voigt_initial)
        print(I_voigt_max)
        # 第一个高斯分量参数
        I_cloud1_initial = 0.6 * np.min(y_fit)
        I_cloud1_min = 0.5 * np.min(y_fit)
        I_cloud1_max = 1.1 * np.max(y_fit)

        mu_cloud1_initial = -0.5
        mu_cloud1_min = -1.0
        if i > 7:
            mu_cloud1_max = -0.5
        else:
            mu_cloud1_max = -0.2

        sigma_cloud1_initial = 0.15
        sigma_cloud1_min = 0.05
        sigma_cloud1_max = 0.25

        # 第二个高斯分量参数
        I_cloud2_initial = 2  # 第二个分量振幅（示例值，可根据蓝翼特征调整）
        I_cloud2_min = 0.01
        I_cloud2_max = 15

        mu_cloud2_initial = -0.5  # 第二个分量的中心位置
        mu_cloud2_min = -1.5
        mu_cloud2_max = -0.05

        sigma_cloud2_initial = 0.10  # 第二个分量的宽度
        sigma_cloud2_min = 0.01
        sigma_cloud2_max = 0.35

        # 常数分量参数
        I_constant_initial = 0.0
        I_constant_min = -0.001
        I_constant_max = 0.002
        S_l_initial = 0.003
        S_l_min = 0.001
        S_l_max = 1

        voigt_single_cloud1 = partial(voigt_single_cloud, I_0=I_0[ha_left_index:ha_right_index + 1])

        popt, pcov = curve_fit(
            voigt_single_cloud1, x_fit, y_fit,
            p0=(
                I_voigt_initial, gamma_voigt_initial, sigma_voigt_initial,
                I_cloud2_initial, mu_cloud2_initial, sigma_cloud2_initial,
                I_constant_initial, S_l_initial
            ),
            bounds=(
                [I_voigt_min, gamma_voigt_min, sigma_voigt_min,
                 I_cloud2_min, mu_cloud2_min, sigma_cloud2_min,
                 I_constant_min, S_l_min],
                [I_voigt_max, gamma_voigt_max, sigma_voigt_max,
                 I_cloud2_max, mu_cloud2_max, sigma_cloud2_max,
                 I_constant_max, S_l_max]
            )
        )

        # 提取双高斯参数
        I_voigt_fit = popt[0]
        gamma_voigt_fit = popt[1]
        sigma_voigt_fit = popt[2]
        I_cloud2_fit = popt[3]
        mu_cloud2_fit = popt[4]
        sigma_cloud2_fit = popt[5]
        I_constant_fit = popt[6]
        S_l_fit = popt[7]
        print("I_cloud2_fit", I_cloud2_fit, "mu_cloud2_fit", mu_cloud2_fit, "sigma_cloud2_fit", sigma_cloud2_fit, "S",
              S_l_fit,"I_voigt_fit",I_voigt_fit,"gamma_voigt_fit",gamma_voigt_fit,"sigma_voigt_fit",sigma_voigt_fit)
        d = 2000 * len(i_indices) / (2 * r_fe + 1) ** 2
        N2 = 7.26e7 * (I_cloud2_fit * sigma_cloud2_fit) / d
        Ne = 3.2e8 * np.sqrt(N2)
        NH = 5e8 * 10 ** (0.5 * np.log10(N2))
        m_H = 1.67e-24
        M = (NH * m_H + 0.0851 * NH * 3.97 * m_H) * d * 1e5

        c = 3e5
        lambda_0 = 6562.8
        line_width = 6
        line_region = (wavelength_point_list_ha_t1 >= 6562.8 - line_width / 2) & (
                wavelength_point_list_ha_t1 <= 6562.8 + line_width / 2)
        lambda_vals = wavelength_point_list_ha_t1[line_region]  # Wavelength range in nm
        I_cont = f_t1_lcont_TR * (2 * TR_width / psr + 1) * (2 * TR_height / psr + 1) / (
                F_t0_l_FD_interpolated[line_region] * (2 * r_fe + 1) ** 2)

        tau_0 = I_cloud2_fit  # Central optical depth
        v_shift = c * mu_cloud2_fit / 6562.8  # Shift velocity in km/s
        W = c * sigma_cloud2_fit / 6562.8  # Line width in km/s
        S = S_l_fit  # Source function (normalized)
        I_lambda = I_0[line_region]  # Reference intensity
        MF0 = M * (2 * r_fe + 1) ** 2 * (756.7e5) ** 2 * 10e-3
        MF.append(MF0)
        print(f"{MF0:.6e}""kg")

        # 保存结果
        mu_cloud2_list.append(mu_cloud2_fit)
        sigma_cloud2_list.append(sigma_cloud2_fit)
        mu_cloud1_list.append(None)
        sigma_cloud1_list.append(None)
        I_cloud2_list.append(I_cloud2_fit)
        I_cloud1_list.append(None)
        I_voigt_list.append(I_voigt_fit)
        S_l.append(S_l_fit)
        gamma_voigt_list.append(gamma_voigt_fit)
        sigma_voigt_list.append(sigma_voigt_fit)

        data_fe_TR = input_fitsfile_fe_1.data[:, int(round(cy_ha + TR_center_y / psr - TR_height / psr)):int(
            round(cy_ha + TR_center_y / psr + TR_height / psr)) + 1,
                     int(round(cx_ha + TR_center_x / psr - TR_width / psr)):int(
                         round(cx_ha + TR_center_x / psr + TR_width / psr)) + 1]
        data_few_TR_image = data_fe_TR[5, :, :]
        data_fe_TR_ave = data_fe_TR.mean(axis=(1, 2))
        data_few_TR_ave = float(
            interp1d(wavelength_point_list_fe_t1, data_fe_TR_ave, kind='cubic', fill_value='extrapolate')(6568))


    # Fit Delta_S_t1_l_interpolated - end
    # ------------------------------

    # ------------------------------
    # Plot the snapshots of Delta_S_t_l - begin
    # TR

    os.chdir(output_dir_figures + '/ds_snapshots')

    x_min = wavelength_point_list_ha_t0[0] - wavelength_interval_ha_t0 / 2 - 6562.8
    x_max = wavelength_point_list_ha_t0[-1] + wavelength_interval_ha_t0 / 2 - 6562.8

    plt.figure(figsize=(3.3464567, 3.3464567 * 0.75), dpi=450)
    plt.title(start_time_str + ' UT', fontsize=6, fontname='Arial')
    plt.xlim(x_min, x_max)
    plt.ylim(-0.003, 0.011)
    plt.xticks(fontsize=5, fontname='Arial')
    plt.yticks([0.000, 0.002, 0.004, 0.006, 0.008, 0.010], fontsize=5, fontname='Arial')
    plt.minorticks_on()
    plt.xlabel('Wavelength ($\mathrm{\AA}$) from $6562.8 \ \mathrm{\AA}$', fontsize=6, fontname='Arial')
    plt.ylabel(r'$\Delta S (t,\lambda,\mathrm{TR})$', fontsize=6, fontname='Arial')

    # 绘制原始数据
    plt.plot(np.array(wavelength_point_list_ha_t0) - 6562.8, Delta_S_t1_l_interpolated, linestyle='-', linewidth=1,
             color='black', zorder=1)
    plt.axvline(0, linewidth=1, linestyle='-.', color='darkviolet', zorder=0.1)
    if i >= fit_starttime_index and i <= fit_endtime_index:
        plt.plot(
            x_fit,
            voigt_single_cloud1(x_fit, I_voigt_fit, gamma_voigt_fit, sigma_voigt_fit,
                                I_cloud2_fit, mu_cloud2_fit, sigma_cloud2_fit,
                                I_constant_fit, S_l_fit),
            linestyle='--',
            linewidth=1,
            color='yellowgreen',
            zorder=2.1,
            label='Fitting result'
        )

        plt.plot(
            x_fit,
            cloud(x_fit, I_cloud2_fit, mu_cloud2_fit, sigma_cloud2_fit, S_l_fit),
            linestyle=':',
            linewidth=1,
            color='magenta',
            zorder=2.4,
            label='Gauss 2 component'
        )

        plt.axvline(mu_cloud2_fit, linewidth=1, linestyle='-.', color='magenta', zorder=0.2)
        # Voigt + 常数分量
        plt.plot(
            x_fit,
            voigt(x_fit, I_voigt_fit, gamma_voigt_fit, sigma_voigt_fit) + I_constant_fit,
            linestyle='--',
            linewidth=1,
            color='darkorange',
            zorder=2.2,
            label='Voigt + Constant components')
        plt.legend(loc='upper left', frameon=True, prop={'family': 'Arial', 'weight': 'normal', 'size': 5})

        plt.axvline(mu_cloud2_fit, linewidth=1, linestyle='-.', color='magenta', zorder=0.2)

    twiny = plt.twiny()
    plt.xlim(x_min / 6562.8 * 3e5, x_max / 6562.8 * 3e5)
    plt.xticks(fontsize=5, fontname='Arial')
    plt.minorticks_on()
    plt.xlabel('Velocity ($\mathrm{km \ s^{-1}}$)', fontsize=6, fontname='Arial')

    # 保存图像
    plt.savefig(input_fitsfilename_ha[:-5] + '_ds_snapshot.png', bbox_inches='tight')
    plt.close()

    # ..............................
    # TReqFD

    os.chdir(output_dir_figures + '/ds_snapshots_TReqFD')

    plt.figure(figsize=(3.3464567, 3.3464567 * 0.75), dpi=450)
    plt.title(start_time_str + ' UT', fontsize=6, fontname='Arial')
    plt.xlim(x_min, x_max)
    plt.ylim(-0.010, 0.025)
    plt.xticks(fontsize=5, fontname='Arial')
    plt.yticks([-0.010, -0.005, 0.000, 0.005, 0.010, 0.015, 0.020, 0.025], fontsize=5, fontname='Arial')
    plt.minorticks_on()
    plt.xlabel('Wavelength ($\mathrm{\AA}$) from $6562.8 \ \mathrm{\AA}$', fontsize=6, fontname='Arial')
    plt.ylabel(r'$\Delta S (t,\lambda,\mathrm{FD})$', fontsize=6, fontname='Arial')

    plt.plot(np.array(wavelength_point_list_ha_t0) - 6562.8, Delta_S_t1_l_interpolated_TReqFD, linestyle='-',
             linewidth=1, color='black', zorder=1)
    plt.axvline(0, linewidth=1, linestyle='-.', color='darkviolet', zorder=0.1)

    # if i >= fit_starttime_index and i <= fit_endtime_index:
    #    plt.plot(x_fit , voigt_gauss(x_fit,I_voigt_fit,gamma_voigt_fit,sigma_voigt_fit,I_gauss_fit,mu_gauss_fit,sigma_gauss_fit,I_constant_fit)                  , linestyle='--' , linewidth=1 , color='yellowgreen' , zorder=2.1 , label='Fitting result')
    #    plt.plot(x_fit , voigt(      x_fit,I_voigt_fit,gamma_voigt_fit,sigma_voigt_fit                                                        ) + I_constant_fit , linestyle='--' , linewidth=1 , color='darkorange'  , zorder=2.2 , label='Voigt component + Constant component')
    #    plt.plot(x_fit , gauss(      x_fit,                                            I_gauss_fit,mu_gauss_fit,sigma_gauss_fit               )                  , linestyle=':'  , linewidth=1 , color='dodgerblue'  , zorder=2.3 , label='Gauss component')
    #    plt.axvline(mu_gauss_fit , line width=1 , linestyle='-.' , color='dodgerblue' , zorder=0.2)
    #    plt.legend(loc='upper left' , frameon=True , prop={'family':'Arial' , 'weight':'normal' , 'size':5})

    twiny = plt.twiny()
    plt.xlim(x_min / 6562.8 * 3e5, x_max / 6562.8 * 3e5)
    plt.xticks(fontsize=5, fontname='Arial')
    plt.minorticks_on()
    plt.xlabel('Velocity ($\mathrm{km \ s^{-1}}$)', fontsize=6, fontname='Arial')

    plt.savefig(input_fitsfilename_ha[:-5] + '_ds_snapshot_TReqFD.png', bbox_inches='tight')
    plt.close()

    # Plot the snapshots of Delta_S_t_l - end
    # ------------------------------

    # ------------------------------
    # Free memory - begin
    del input_fitsfile_ha
    del input_fitsfile_ha_0
    del input_fitsfile_ha_1
    del input_fitsfile_fe
    del input_fitsfile_fe_0
    del input_fitsfile_fe_1
    gc.collect()
    time.sleep(1)
    # Free memory - end
    # ------------------------------

    print('------------------------------')
# end for

Delta_S_t_l_interpolated = Delta_S_t_l_interpolated.T
Delta_S_t_l_interpolated_TReqFD = Delta_S_t_l_interpolated_TReqFD.T

# ==============================
# Time calculation - begin
start_time_interval_s_list = []
start_time_s_list = []

for i in range(len(start_time_str_list)):
    dt_start_time = datetime.strptime(start_time_str_list[i], dt_format)
    if i == 1:  # 找到索引1的时间点
        dt_start_time_0 = dt_start_time
    if i > 0:  # 计算与前一个时间点的时间间隔
        start_time_interval_s = float((dt_start_time - dt_start_time_last).total_seconds())
        start_time_interval_s_list.append(start_time_interval_s)
    dt_start_time_last = dt_start_time

# 第二次遍历：计算每个时间点相对于索引1的偏移量
for i in range(len(start_time_str_list)):
    dt_start_time = datetime.strptime(start_time_str_list[i], dt_format)
    # 计算与索引1的时间差
    start_time_s = float((dt_start_time - dt_start_time_0).total_seconds())
    start_time_s_list.append(start_time_s)

# 转换为NumPy数组
start_time_interval_s_array = np.array(start_time_interval_s_list)
start_time_s_array = np.array(start_time_s_list)

dt_start_0 = datetime.strptime(start_time_str_0, dt_format)
dt_start_f = datetime.strptime(start_time_str_f, dt_format)
print(start_time_interval_s_list)
a = start_time_interval_s_list[-1]
b = start_time_interval_s_list[-2]

if start_time_interval_s_list[-1] == a:
    duration_s = float((dt_start_f - dt_start_0).total_seconds()) + b
elif start_time_interval_s_list[-1] == b:
    duration_s = float((dt_start_f - dt_start_0).total_seconds()) + a

# Time calculation - end
# ==============================

# ==============================
# Plot the dynamic spectra and EW-t - begin
# ..............................
# TR - begin
# ------------------------------
# Plot the dynamic spectra - begin

os.chdir(output_dir_figures + '/ds')

x_min = 0
x_max = duration_s
y_min = wavelength_point_list_ha_t0[0] - wavelength_interval_ha_t0 / 2 - 6562.8
y_max = wavelength_point_list_ha_t0[-1] + wavelength_interval_ha_t0 / 2 - 6562.8

fig = plt.figure(figsize=(3.3464567+1, 3.3464567 * 0.75), dpi=450)

plt.xlim(x_min, x_max)
plt.ylim(y_min, y_max)
plt.yticks(fontsize=5, fontname='Arial')
plt.minorticks_on()
# plt.xlabel('Time ($\mathrm{s}$) ' + 'since {} UT'.format(start_time_str_0) , fontsize=6 , fontname='Arial')
plt.ylabel('Wavelength ($\mathrm{\AA}$) from $6562.8 \ \mathrm{\AA}$', fontsize=6, fontname='Arial')

ds_imshow = plt.imshow(Delta_S_t_l_interpolated, origin='lower', cmap='PuOr_r', aspect='auto',
                       vmin=-np.max(Delta_S_t_l_interpolated) * 1.75, vmax=np.max(Delta_S_t_l_interpolated) * 1.75,
                       extent=(x_min, x_max, y_min, y_max))
plt.axhline(0, linewidth=1, linestyle='--', color='yellowgreen')

plt.plot(start_time_s_array[fit_starttime_index: fit_endtime_index + 1] + (a+b) / 2, mu_cloud2_list, linestyle='-',
         linewidth=1, color='palevioletred', label='Cloud fit Center')
plt.scatter(start_time_s_array[fit_starttime_index: fit_endtime_index + 1] + (a+b) / 2, mu_cloud2_list, s=3,
            color='palevioletred')
plt.errorbar(start_time_s_array[fit_starttime_index: fit_endtime_index + 1] + (a+b) / 2, mu_cloud2_list,
             yerr=sigma_cloud2_list, color='palevioletred', fmt='none', linewidth=0.5, capsize=0.5, capthick=0.5)

# 添加图例
plt.legend(loc='upper right', fontsize=5, frameon=True)
twinx = plt.twinx()
plt.ylim(y_min / 6562.8 * 3e5, y_max / 6562.8 * 3e5)
plt.yticks(fontsize=5, fontname='Arial')
plt.minorticks_on()
plt.ylabel('Velocity ($\mathrm{km \ s^{-1}}$)', fontsize=6, fontname='Arial')

cbaxes = fig.add_axes([0.999, 0.15, 0.01, 0.6])  # Adjust these values for better positioning
cbaxes.yaxis.set_ticks_position("right")
ds_colorbar = plt.colorbar(ds_imshow, cax=cbaxes, orientation='vertical', extend='max')
ds_colorbar.outline.set_edgecolor('black')
ds_colorbar.ax.tick_params(axis='x', colors='black')
ds_colorbar.ax.tick_params(axis='y', colors='black')
ds_colorbar.set_label(r'$\Delta S (t,\lambda,\mathrm{TR})$', color='black', fontsize=5, fontname='Arial')
# ds_colorbar.set_ticks([-0.005 , 0.000 , 0.005])
ds_colorbar.minorticks_on()
ds_colorbar_labels = ds_colorbar.ax.get_xticklabels() + ds_colorbar.ax.get_yticklabels()
[label.set_fontproperties(FontProperties(size=3, family='Arial')) for label in ds_colorbar_labels]
[label.set_color('black') for label in ds_colorbar_labels]

# Plot the dynamic spectra - end
# ------------------------------

plt.savefig('Delta_S(t,lambda).pdf', bbox_inches='tight')
plt.close()

# TR - end
# ..............................

# ..............................
# TReqFD - begin

# ------------------------------
# Plot the dynamic spectra - begin

os.chdir(output_dir_figures + '/ds_TReqFD')

fig = plt.figure(figsize=(3.3464567+1, 3.3464567 * 0.75), dpi=450)

plt.xlim(x_min, x_max)
plt.ylim(y_min, y_max)
plt.yticks(fontsize=5, fontname='Arial')
plt.minorticks_on()
plt.ylabel('Wavelength ($\mathrm{\AA}$) from $6562.8 \ \mathrm{\AA}$', fontsize=6, fontname='Arial')

ds_imshow = plt.imshow(Delta_S_t_l_interpolated_TReqFD, origin='lower', cmap='PuOr_r', aspect='auto',
                       vmin=-np.max(Delta_S_t_l_interpolated_TReqFD) * 0.75,
                       vmax=np.max(Delta_S_t_l_interpolated_TReqFD) * 0.75, extent=(x_min, x_max, y_min, y_max))
plt.axhline(0, linewidth=1, linestyle='--', color='yellowgreen')

# plt.plot(    start_time_s_array[fit_starttime_index : fit_endtime_index+1] +71/2 , mu_gauss_list , linestyle='-' , linewidth=1 , color='dodgerblue')
# plt.scatter( start_time_s_array[fit_starttime_index : fit_endtime_index+1] +71/2 , mu_gauss_list , s=3 , color='dodgerblue')
# plt.errorbar(start_time_s_array[fit_starttime_index : fit_endtime_index+1] +71/2 , mu_gauss_list , yerr=sigma_gauss_list , color='dodgerblue' , fmt='none' , linewidth=0.5 , capsize=0.5 , capthick=0.5)

twinx = plt.twinx()
plt.ylim(y_min / 6562.8 * 3e5, y_max / 6562.8 * 3e5)
plt.yticks(fontsize=5, fontname='Arial')
plt.minorticks_on()
plt.ylabel('Velocity ($\mathrm{km \ s^{-1}}$)', fontsize=6, fontname='Arial')

cbaxes = fig.add_axes([0.999, 0.15, 0.01, 0.6])  # Adjust these values for better positioning
cbaxes.yaxis.set_ticks_position("right")
ds_colorbar = plt.colorbar(ds_imshow, cax=cbaxes, orientation='vertical', extend='max')
ds_colorbar.outline.set_edgecolor('black')
ds_colorbar.ax.tick_params(axis='x', colors='black')
ds_colorbar.ax.tick_params(axis='y', colors='black')
ds_colorbar.set_label(r'$\Delta S (t,\lambda,\mathrm{FD})$', color='black', fontsize=5, fontname='Arial')
# ds_colorbar.set_ticks([-0.005 , 0.000 , 0.005])
ds_colorbar.minorticks_on()
ds_colorbar_labels = ds_colorbar.ax.get_xticklabels() + ds_colorbar.ax.get_yticklabels()
[label.set_fontproperties(FontProperties(size=3, family='Arial')) for label in ds_colorbar_labels]
[label.set_color('black') for label in ds_colorbar_labels]

# Plot the dynamic spectra - end
# ------------------------------

plt.savefig('Delta_S(t,lambda)_TReqFD.pdf', bbox_inches='tight')
plt.close()

# TReqFD - end
# ..............................
plt.rcParams['font.family'] = 'Arial'   # 可以改为 'Times New Roman' 等期刊常用字体
plt.rcParams['font.size'] = 5         # 统一控制字体大小
plt.rcParams['figure.dpi'] = 300       # 设置分辨率，保证打印/出版质量
os.chdir(output_dir_figures + '/EW')
print(len(start_time_s_array),len(ew_list),len(MF))
fig, ax = plt.subplots(figsize=(7, 4))
min_idx = np.argmax(delta_ew_list)
min_val = delta_ew_list[min_idx]
ax.scatter(start_time_s_array[1:], delta_ew_list[1:], marker='o', edgecolor='black', s=30, alpha=0.8, label='Mass')
ax.set_xlabel('Time (s)', fontsize=12)
ax.set_ylabel(r"$\Delta \ \mathrm{EW}$" , fontsize=12)
ax.set_title(r"$\Delta \ \mathrm{EW}\ \mathrm{VS}\ \mathrm{Time}$" , fontsize=14)
ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
ax.legend(frameon=False)
ax.grid(True, linestyle='--', alpha=0.5)
ax.text(min_idx , min_val, color='red', fontsize=6, fontname='Arial', ha='left', va='bottom')
plt.tight_layout()
plt.savefig('ew.pdf', bbox_inches='tight')
plt.close()
# Plot the dynamic spectra and EW-t - end
# ==============================
# 假设所有列表长度相同，这里用len(mu_cloud2_list)举例
with open('fit_results.txt', 'w', encoding='utf-8') as f:
    # 写入表头，可根据需要自行调整
    f.write("Iteration\tmu_cloud2\tsigma_cloud2\tI_cloud2\tI_voigt\tS_l\tgamma_voigt\tsigma_voigt\tMF\n")
    # 写入每一行数据
    for i in range(len(mu_cloud2_list)):
        # 获取各项数据
        mu_cloud2 = mu_cloud2_list[i]
        sigma_cloud2 = sigma_cloud2_list[i]
        I_cloud2 = I_cloud2_list[i]
        I_voigt = I_voigt_list[i]
        S_l_val = S_l[i]
        gamma_voigt = gamma_voigt_list[i]
        sigma_voigt = sigma_voigt_list[i]
        M_fit = MF[i]
        # 将None值转成字符串"None"，或自行在此做判断处理
        line = f"{i}\t{mu_cloud2}\t{sigma_cloud2}\t{I_cloud2}\t{I_voigt}\t{S_l_val}\t{gamma_voigt}\t{sigma_voigt}\t{MF}\n"
        f.write(line)










