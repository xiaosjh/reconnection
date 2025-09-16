import numpy as np
from scipy.optimize import curve_fit
from astropy.io import fits
import matplotlib.pyplot as plt
import os
from itertools import product
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import pandas as pd

# 假设 lam, I, 和 I0 是你的数据数组
# 这里的 lam 是波长数组，I 是通过暗条的光谱强度，I0 是背景光谱强度

input_dir_ha = r"/Volumes/T7Shield/CHASEdata/CME1/20230313/data00/Ha"
input_dir_fe = r"/Volumes/T7Shield/CHASEdata/CME1/20230313/data00/Fe"
output_dir_figures = r"/Users/yvonne/Documents/sunas/CME1/output0313/Figures3"

image_half_side = 1100  # The half side length of the cropped image. Unit: pixels.
psr = 1.0436  # Pixel spatial resolution. Unit: arcsec per pixel.

TR_center = [-130,
             640]  # The center of the target region. Relative to the center of the solar disk. [relative_x , relative_y]. Unit: arcsec.
TR_half_side = 120  # The half side length of the target region. Unit: arcsec.
TR_width = 140
TR_height = 140

ha_left_index = 30  # The index of the left  border of Ha line. Minimum: 30
ha_right_index = 117  # The index of the right border of Ha line. Maximum: 117

fit_starttime_index = 1  # The index of the start time for fitting. Minimum: 0
fit_endtime_index = 14  # The index of the end   time for fitting. Maximum: 22

dt_format = '%Y-%m-%d %H:%M:%S'
input_filename_ha_list = os.listdir(input_dir_ha)
input_filename_fe_list = os.listdir(input_dir_fe)
input_fitsfilename_ha_list = []
input_fitsfilename_fe_list = []
TR_center_x = TR_center[0]
TR_center_y = TR_center[1]

start_time_str_list = []
tau0_opt0 = []
delta_lambda_opt0 = []
S_opt0 = []
mu_opt0 = []
MF1=[]
length_1 = []
for i in range(len(input_filename_ha_list)):
    if input_filename_ha_list[i][0:3] == 'RSM' and input_filename_ha_list[i][-7:] == 'HA.fits':
        input_fitsfilename_ha_list.append(input_filename_ha_list[i])
input_fitsfilename_ha_list = sorted(input_fitsfilename_ha_list)
for i in range(len(input_filename_fe_list)):
    if input_filename_fe_list[i][0:3] == 'RSM' and input_filename_fe_list[i][-7:] == 'FE.fits':
        input_fitsfilename_fe_list.append(input_filename_fe_list[i])
input_fitsfilename_fe_list = sorted(input_fitsfilename_fe_list)
print('Total number of Ha FITS files:', len(input_fitsfilename_ha_list))
print('Total number of Fe FITS files:', len(input_fitsfilename_fe_list))
for i in range(len(input_fitsfilename_ha_list)):
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
    a = input_fitsfile_ha_1.header[0]
    print(a)
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
    if i == 0:
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

    data_hac_image = input_fitsfile_ha_1.data[68,
                     int(round(cy_ha - image_half_side, 0)):int(round(cy_ha + image_half_side, 0)) + 1,
                     int(round(cx_ha - image_half_side, 0)):int(
                         round(cx_ha + image_half_side)) + 1]  # hac: Ha line center

    data_ha_TR = input_fitsfile_ha_1.data[:, int(round(cy_ha + TR_center_y / psr - TR_height / psr)):int(
        round(cy_ha + TR_center_y / psr + TR_height / psr)) + 1,
                 int(round(cx_ha + TR_center_x / psr - TR_width / psr)):int(
                     round(cx_ha + TR_center_x / psr + TR_width / psr)) + 1]

    data_hac_TR_image = np.array(data_ha_TR[68, :, :], dtype=np.float64) # 复制第 68 层数据
    for i in range(1, 41):  # 从 1 到 10 进行累加
        data_hac_TR_image += data_ha_TR[68 - i, :, :]
        data_hac_TR_image += data_ha_TR[68 + i, :, :]
    data_hac_TR_image /= 120

    plt.figure()
    plt.minorticks_on()
    plt.xlabel('Solar X (arcsec) from center', fontsize=6, fontname='Arial')
    plt.ylabel('Solar Y (arcsec) from center', fontsize=6, fontname='Arial')
    plt.imshow(data_hac_TR_image, origin='lower', cmap='afmhot', vmin=0, vmax=4 * data_hac_image.mean(), extent=(
    TR_center_x - TR_width, TR_center_x + TR_width, TR_center_y - TR_height, TR_center_y + TR_height))

    indices = np.where(input_fitsfile_ha_1.data > 0)
    low_flux_indices = np.where((data_hac_TR_image < 520) & (data_hac_TR_image > 0))
    i_indices, j_indices = low_flux_indices
    high_flux_indices = np.where(data_hac_TR_image >= 520)
    k_indices, l_indices = high_flux_indices

    N_y, N_x = data_hac_TR_image.shape
    x_min = TR_center_x - TR_width
    x_max = TR_center_x + TR_width
    y_min = TR_center_y - TR_height
    y_max = TR_center_y + TR_height
    x_coords = x_min + (j_indices / N_x) * (x_max - x_min)
    y_coords = y_min + (i_indices / N_y) * (y_max - y_min)

    plt.scatter(x_coords, y_coords, color='red', s=1)
    # 显示图像
    plt.show()
    f_d_l_TR = data_ha_TR[:, i_indices, j_indices].mean(axis=1)
    pix_dark = len(i_indices)
    f_o_l_TR = data_ha_TR[:, k_indices, l_indices].mean(axis=1)
    pix_others = len(k_indices)
    f_l_TR = data_ha_TR[:].mean(axis=(1,2))

    x_fit = np.array(wavelength_point_list_ha_t1) - 6562.8

    data_fe_TR = input_fitsfile_fe_1.data[:, int(round(cy_ha + TR_center_y / psr - TR_height / psr)):int(
        round(cy_ha + TR_center_y / psr + TR_height / psr)) + 1,
                 int(round(cx_ha + TR_center_x / psr - TR_width / psr)):int(
                     round(cx_ha + TR_center_x / psr + TR_width / psr)) + 1]
    data_few_TR_image = data_fe_TR[5, :, :]
    data_fe_TR_ave = data_fe_TR.mean(axis=(1, 2))
    data_few_TR_ave = float(
            interp1d(wavelength_point_list_fe_t1, data_fe_TR_ave, kind='cubic', fill_value='extrapolate')(6568))

    def model_spectrum(lam, tau0, delta_lambda, S, mu):
        tau_lam = tau0 * np.exp(-0.5*((lam-mu)**2) / delta_lambda**2)
        return f_o_l_TR* np.exp(-tau_lam) + S * (1 - np.exp(-tau_lam))

    initial_tau0 = 5
    initial_delta_lambda = 0.5
    initial_S = 100
    initial_mu = -1
    popt, pcov = curve_fit(model_spectrum, (np.array(wavelength_point_list_ha_t1) - 6562.8), f_d_l_TR ,
                           p0=[initial_tau0, initial_delta_lambda, initial_S,initial_mu], bounds = ([0.3,0.1, 0.001,-3],[10,1.2,500,-0.01]),maxfev = 2000)

    # 输出最优参数
    tau0_opt, delta_lambda_opt, S_opt, mu_opt = popt
    tau0_opt0.append(tau0_opt)
    delta_lambda_opt0.append(delta_lambda_opt)
    S_opt0.append(S_opt,)
    mu_opt0.append(mu_opt)
    length_1.append(len(i_indices))
    print("Optimal τ₀:", tau0_opt)
    print("Optimal Δλ:", delta_lambda_opt)
    print("Optimal S:", S_opt)
    print("Optimal mu:", mu_opt)

    plt.figure(figsize=(3.3464567, 3.3464567 * 0.75), dpi=450)
    plt.suptitle(start_time_str + ' UT', fontsize=5, fontname='Arial')
    plt.yticks(fontsize=5, fontname='Arial')
    plt.xticks([-3.0,-2.0,-1.0,0.0,1.0,2.0], fontsize=5, fontname='Arial')
    plt.axvline(0, linewidth=1, linestyle='-.', color='red', zorder=0.1)
    plt.axvline(mu_opt, linewidth=1, linestyle='-.', color='blue', zorder=0.1)
    plt.plot(np.array(wavelength_point_list_ha_t1) - 6562.8, f_d_l_TR, label="f_filament",
             linestyle='--', color='green', zorder=0.2)
    plt.plot(np.array(wavelength_point_list_ha_t1) - 6562.8, f_o_l_TR, label="f_background", linestyle='-.',
             color='red', zorder=0.2)
    plt.plot(x_fit, model_spectrum(x_fit, tau0_opt, delta_lambda_opt, S_opt, mu_opt), label="f_fit", linestyle=':',
           color='blue', zorder=0.2)
    plt.legend(loc='upper right', fontsize=5, frameon=False)
    plt.savefig(output_dir_figures + '/mass/' + input_fitsfilename_ha[:-5] + '_mass.pdf', dpi=450, bbox_inches='tight')
    plt.show()

    d = 2000
    N2 = 7.26e7 * (tau0_opt * delta_lambda_opt) / d
    Ne = 3.2e8 * np.sqrt(N2)
    NH = 5e8 * 10**(0.5 * np.log10(N2))
    m_H = 1.67e-24
    M = (NH * m_H + 0.0851 * NH * 3.97 * m_H) * d * 1e5
    MT = len(i_indices) * M * (756.7e5)**2 * 10e-3
    print(f"{NH:.6e}","cm^-3",f"{MT:.6e}""kg")
    MF1.append(MT)

os.chdir(output_dir_figures + '/EW')
plt.figure(figsize=(6, 4))
max_idx = np.argmax(MF1)
max_val = MF1[max_idx]
plt.scatter(range(len(MF1)), MF1, label='Mass')
plt.xlabel('Time (s)', fontsize=6, fontname='Arial')
plt.ylabel('Mass', fontsize=6, fontname='Arial')
plt.title('Mass vs Time(every second opt)', fontsize=7, fontname='Arial')
# 标记最大值点
plt.scatter(max_idx , MF1[max_idx] , color='red', s=20, zorder=5, label='Max')
plt.text(max_idx, MF1[max_idx], f'Max = {max_val:.2e}\ntau = {tau0_opt0[max_idx]}\ndelta_lambda ={delta_lambda_opt0[max_idx]}\nS = {S_opt0[max_idx]}\nmu = {mu_opt0[max_idx]}'
                                        f'\ntaum = {np.mean(tau0_opt0)}\ndelta_lambda ={np.mean(delta_lambda_opt0)}\nS = {np.mean(S_opt0)}\nmu = {np.mean(mu_opt0)}', color='red',
         fontsize=6, fontname='Arial', ha='left', va='bottom')
plt.legend(fontsize=5)
plt.grid(True)
plt.savefig('Mass vs Time(every second opt).pdf')
plt.show()

d_range = np.linspace(1800, 2200, 10)
indices_range = np.linspace(0.9* length_1[max_idx], 1.5* length_1[max_idx] , 10)

results = []
for d ,indices in product(d_range,indices_range):
    try:
        N2 = 7.26e7 * (tau0_opt0[max_idx] * delta_lambda_opt0[max_idx]) / d
        Ne = 3.2e8 * np.sqrt(N2)
        NH = 5e8 * 10 ** (0.5 * np.log10(N2))
        m_H = 1.67e-24
        M = (NH * m_H + 0.0851 * NH * 3.97 * m_H) * d * 1e5
        MT = indices * M * (756.7e5) ** 2 * 10e-3
        results.append((d, indices, float(MT)))
    except Exception as e:
        print(f"Error at d={d}, indices = {indices}: {e}")
        continue

# 转成 DataFrame
df = pd.DataFrame(results, columns=["d","indices", "MT"])
df_sorted = df.sort_values(by="MT", ascending=False).reset_index(drop=True)
avg_mt = df_sorted["MT"].mean()
avg_row = pd.DataFrame({
    "d": ["平均值"],
    "indices": [np.nan],
    "MT": [avg_mt]
})
df_final = pd.concat([df_sorted, avg_row], ignore_index=True)

# 写入 CSV，MT列格式化为科学计数法
df_final["MT"] = df_final["MT"].apply(lambda x: f"{x:.4e}" if pd.notnull(x) else "")
df_final.to_csv("truemass.csv", index=False)
# 显示浮动范围
print("\n--- MT0 浮动范围统计 ---")
print(f"最小值：{float(df['MT'].min()):.3e} kg")
print(f"最大值：{float(df['MT'].max()):.3e} kg")
print(f"平均值：{float(df['MT'].mean()):.3e} kg")