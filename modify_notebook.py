import json

# 读取notebook
with open('draw_image3.ipynb', 'r', encoding='utf-8') as f:
    nb = json.load(f)

# 找到并修改目标cell（Cell 20）
target_cell_idx = None
for idx, cell in enumerate(nb['cells']):
    if cell['cell_type'] == 'code':
        source = ''.join(cell['source'])
        if 'contour=ax3.contour' in source:
            target_cell_idx = idx
            break

if target_cell_idx is not None:
    # 新的源代码
    new_source = """def double_cloud_vel(params, vel, Icont):
    ltauL, ltauU, vL, vU, lWl, lWu, SL, SU = params
    tau0L = np.exp(ltauL)
    tau0U = np.exp(ltauU)
    WL = np.exp(lWl)
    WU = np.exp(lWu)
    
    tauL = tau0L * np.exp(-((vel - vL)**2) / (WL**2))
    tauU = tau0U * np.exp(-((vel - vU)**2) / (WU**2))
    
    I = Icont * np.exp(-(tauL + tauU)) \\
        + SL * (1 - np.exp(-tauL)) * np.exp(-tauU) \\
        + SU * (1 - np.exp(-tauU))
    return Iimport matplotlib as mpl
mpl.rcParams['text.usetex'] = False
#数据准备
dir=r'C:\\Learning\\PHD1st\\magnetic_reconnecion\\data\\CHASE_Ha\\RSM20240618T211711_0000_HA.fits'
rsm=fits.open(dir)
left_index=44
right_index=96
lam00 = rsm[1].header['CRVAL3'] + np.arange(rsm[1].header['NAXIS3']) * rsm[1].header['CDELT3']
lam=lam00[left_index:right_index]#±1.5A
I0=rsm[1].data[left_index:right_index,797:826,1540:1569].mean(axis=(1,2))
print(lam00[39],lam00[44],lam00[54],lam00[96],lam00[101])
lam0 = 6562.8
c_km = 299792.458  # 光速 (km/s)
# 1. wavelength → velocity
lam=lam00[left_index:right_index]#±1.5A
vel = (lam - lam0) / lam0 * c_km  # km/s#画图
hawing = rsm[1].data[44:54, :, :].mean(axis=0)
coord_HIS = SkyCoord(0 * u.arcsec, 0 * u.arcsec, obstime = '2024-06-18 21:17:11', observer = 'earth', \\
                     frame = frames.Helioprojective)
headerwing = sunpy.map.make_fitswcs_header(hawing, coord_HIS,
                                       reference_pixel = \\
                                       [rsm[1].header['CRPIX1'], rsm[1].header['CRPIX2']] * u.pixel,
                                       scale = [0.5218 * 2, 0.5218 * 2] * u.arcsec / u.pixel,
                                       telescope = 'CHASE', instrument = 'RSM')
hawing_map = sunpy.map.Map(hawing, headerwing)left_bottom=SkyCoord(Tx=200*u.arcsec,Ty=-520*u.arcsec,frame=hawing_map.coordinate_frame)
top_right=SkyCoord(Tx=420*u.arcsec,Ty=-300*u.arcsec,frame=hawing_map.coordinate_frame)
sub_hawing_map=hawing_map.submap(left_bottom,top_right=top_right)#aia131图像
dir2=r'C:\\Learning\\PHD1st\\magnetic_reconnecion\\data\\AIA2\\131\\aia.lev1_euv_12s.2024-06-18T211708Z.131.image_lev1.fits'
aia=sunpy.map.Map(dir2)
left_bottom=SkyCoord(Tx=200*u.arcsec,Ty=-520*u.arcsec,frame=aia.coordinate_frame)
top_right=SkyCoord(Tx=420*u.arcsec,Ty=-300*u.arcsec,frame=aia.coordinate_frame)
sub_aia=aia.submap(left_bottom,top_right=top_right)
#----------------------------------------------------------------------------------------------------------------------------
#图1
fig=plt.figure(figsize=(9.5,7))
gs=gridspec.GridSpec(2,1,hspace=0.25,height_ratios=[1,0.9])
gs_top=gridspec.GridSpecFromSubplotSpec(1,2,subplot_spec=gs[0],wspace=0.15,width_ratios=[1,1])
old_norm=sub_aia.plot_settings['norm']
old_norm.vmin=0.01
old_norm.vmax=2000
ax1=fig.add_subplot(gs_top[0],projection=sub_aia)
sub_aia.plot(axes=ax1,norm=old_norm)
#------------------------------------------------------------------------------------------------------------------
#图2
ax2=fig.add_subplot(gs_top[1],projection=sub_hawing_map)
sub_hawing_map.plot(axes=ax2,cmap = 'afmhot', vmin = 0, vmax = 2 * sub_hawing_map.data.mean())

#point1是宁静区左下角
point1=SkyCoord(Tx=380*u.arcsec,Ty=-370*u.arcsec,frame=hawing_map.coordinate_frame)
xpix1,ypix1=hawing_map.world_to_pixel(point1)
print(f'宁静区左下角的位置  y:{ypix1.value:.2f},x:{xpix1.value:.2f}')
#point2是宁静区右上角
point2=SkyCoord(Tx=410*u.arcsec,Ty=-340*u.arcsec,frame=hawing_map.coordinate_frame)
xpix2,ypix2=hawing_map.world_to_pixel(point2)
print(f'宁静区右上角的位置  y:{ypix2.value:.2f},x:{xpix2.value:.2f}')
sub_hawing_map.draw_quadrangle(point1,top_right=point2,axes=ax2,edgecolor='black')

jet_point=(1520,791)
point3=hawing_map.wcs.pixel_to_world(jet_point[0],jet_point[1])
#ax2.plot_coord(point3,marker='o',markersize=3)
xpix3,ypix3=hawing_map.world_to_pixel(point3)
print(f'喷流点选取的位置  y:{ypix3.value:.2f},x:{xpix3.value:.2f}')
#point4是喷流选取左下角
jet_point4=(1515,789)
point4=hawing_map.wcs.pixel_to_world(jet_point4[0],jet_point4[1])
xpix4,ypix4=hawing_map.world_to_pixel(point4)
print(f'喷流选取的位置  y:{ypix4.value:.2f},x:{xpix4.value:.2f}')
#point5是喷流选取右上角
jet_point5=(1525,806)
point5=hawing_map.wcs.pixel_to_world(jet_point5[0],jet_point5[1])
#point5=SkyCoord(Tx=364*u.arcsec,Ty=-363*u.arcsec,frame=hawing_map.coordinate_frame)
xpix5,ypix5=hawing_map.world_to_pixel(point5)
print(f'喷流选取的位置  y:{ypix5.value:.2f},x:{xpix5.value:.2f}')

sub_hawing_map.draw_quadrangle(point4,top_right=point5,axes=ax2,edgecolor='white')
sub_aia.draw_quadrangle(point4,top_right=point5,axes=ax1,edgecolor='white')

#在ax2上添加喷流区域的contour（level=800）
jet_region=sub_hawing_map.submap(point4,top_right=point5)
contour=ax2.contour(jet_region.data,levels=[800],colors=['blue'],transform=ax2.get_transform(jet_region.wcs))

#----------------------------------------------------------------------------------------------
#图3（原图4，现在变成图c）
gs_down=gridspec.GridSpecFromSubplotSpec(1,1,subplot_spec=gs[1])
ax3=fig.add_subplot(gs_down[0])
#获取待拟合光谱
#由contour生成mask
cs = contour
from matplotlib.path import Path
# 初始化 mask
mask = np.zeros_like(jet_region.data, dtype=bool)
# 构造像素坐标网格
ny, nx = jet_region.data.shape
X, Y = np.meshgrid(np.arange(nx), np.arange(ny))
points = np.vstack((X.flatten(), Y.flatten())).T  # (N, 2)
# 遍历所有封闭区域
for path in cs.get_paths():
    poly = Path(path.vertices)
    mask |= poly.contains_points(points).reshape(ny, nx)
data_jet=rsm[1].data[left_index:right_index,789:807,1515:1526]
data_mask=np.tile(~mask,(right_index-left_index,1,1))
jet_masked=np.ma.array(data_jet,mask=data_mask)
I=jet_masked.mean(axis=(1,2))

#获取拟合参数
x_opt=np.load('x_opt.npy')
ax3.errorbar(lam,I,yerr=np.sqrt(I),fmt='o',capsize=2,markersize=2,label='raw data')
y_fit = double_cloud_vel(x_opt, vel, I0)
ax3.plot(lam, y_fit, label='fitting curve')
ax3.axvline(x=6562.8,linestyle='--',color='black')
ax3.set_xlabel(r'Wavelength $(\mathrm{\\AA})$')
ax3.set_ylabel('Intensity (Count)')
ax3.legend(bbox_to_anchor=(0.45, 0.77))

#---------------------------------------------------------------------------------------------------
#画图乱七八糟的部分
ax_all=[ax1,ax2]
for ax in ax_all:
    ax.set_xlabel('X (arcsec)',labelpad=1.2)
    ax.set_ylabel('Y (arcsec)')
    ax.tick_params(axis='x',which='both',bottom=True,labelbottom=True,top=False,labeltop=False)
    ax.tick_params(axis='y',which='both',left=True,labelleft=True,right=False,labelright=False)
    ax.coords[0].grid(draw_grid=False)
    ax.coords[1].grid(draw_grid=False)
    ax.coords[0].set_major_formatter('s',show_decimal_unit=False)
    ax.coords[1].set_major_formatter('s',show_decimal_unit=False)
    ax.set_title('')
ax_all[0].text(0.02,0.93,'(a)  AIA 131  21:17:06 UT',color='white',transform=ax_all[0].transAxes,fontsize=12)
ax_all[1].text(0.02,0.85,r'(b)  CHASE H$\\mathrm{\\alpha}$ 6561.55-6562.03$\\mathrm{\\AA}$'+'\n       21:17:11 UT',color='white',transform=ax_all[1].transAxes,fontsize=12)
ax3.text(0.02,0.9,'(c)',color='black',transform=ax3.transAxes,fontsize=12)
ax3.text(0.69,0.05,r'$\\mathrm{v_U}$'+'=-15.65 '+r'$\\mathrm{km \\ s^{-1}}$'+'\n'+r'$\\mathrm{v_L} $'+'= 20.78 '+r'$\\mathrm{km \\ s^{-1}}$',color='black',
         fontsize=10,transform=ax3.transAxes)
#plt.savefig(r'C:\\Learning\\PHD1st\\magnetic_reconnecion\\paper\\fig5\\double_cloud.pdf',format='pdf',bbox_inches='tight', pad_inches=0.1)"""
    
    # 转换为列表格式（notebook中source是列表）
    new_source_list = new_source.split('\n')
    
    # 更新cell的source
    nb['cells'][target_cell_idx]['source'] = new_source_list
    
    print(f"Modified cell {target_cell_idx}")
    
    # 保存修改后的notebook
    with open('draw_image3.ipynb', 'w', encoding='utf-8') as f:
        json.dump(nb, f, ensure_ascii=False, indent=1)
    
    print("Notebook saved successfully!")
else:
    print("Target cell not found!")
