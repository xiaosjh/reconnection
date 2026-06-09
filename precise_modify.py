import json

# 读取备份
with open('draw_image3.ipynb.bak', 'r', encoding='utf-8') as f:
    nb = json.load(f)

cell_source_list = nb['cells'][20]['source']
cell_text = ''.join(cell_source_list)

# 分割关键部分进行精确编辑
# 1. 找到"sub_aia.draw_quadrangle(point4,top_right=point5,axes=ax1,edgecolor='white')"之后的部分

# 用marker来帮助精确定位
MARKER_START = "#------------------------------"
MARKER_END = "#plt.savefig"

# 找到关键位置
marker_pos_start = cell_text.find(MARKER_START)
marker_pos_end = cell_text.find(MARKER_END)

if marker_pos_start == -1 or marker_pos_end == -1:
    print("ERROR: markers not found!")
    exit(1)

# 分割文本
before_marker = cell_text[:marker_pos_start]
after_marker = cell_text[marker_pos_end:]

# 构建新的中间部分
new_middle = """#----------------------------------------------------------------------------------------------
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
ax3.set_xlabel(r'Wavelength $(\mathrm{\AA})$')
ax3.set_ylabel('Intensity (Count)')
ax3.legend(bbox_to_anchor=(0.45, 0.77))

#---------------------------------------------------------------------------------------------------
#画图乱七八糟的部分
ax_all=[ax1,ax2]
for ax in ax_all:
    #ax.set_xlabel('X (arcsecs)')
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
ax_all[1].text(0.02,0.85,r'(b)  CHASE H$\mathrm{\alpha}$ 6561.55-6562.03$\mathrm{\AA}$'+'\n       21:17:11 UT',color='white',transform=ax_all[1].transAxes,fontsize=12)
ax3.text(0.02,0.9,'(c)',color='black',transform=ax3.transAxes,fontsize=12)
ax3.text(0.69,0.05,r'$\mathrm{v_U}$'+'=-15.65 '+r'$\mathrm{km \ s^{-1}}$'+'\n'+r'$\mathrm{v_L} $'+'= 20.78 '+r'$\mathrm{km \ s^{-1}}$',color='black',
         fontsize=10,transform=ax3.transAxes)
"""

# 也需要修改figsize
new_text = before_marker + new_middle + after_marker

# 替换figsize
new_text = new_text.replace('figsize=(9.5,8)', 'figsize=(9.5,7)')

# 转换为source list格式（保持原有的换行结构）
new_source_list = new_text.split('\n')

# 更新cell
nb['cells'][20]['source'] = new_source_list

# 保存
with open('draw_image3.ipynb', 'w', encoding='utf-8') as f:
    json.dump(nb, f, ensure_ascii=False, indent=1)

print("修改成功!")
print(f"新文本长度: {len(new_text)} characters")
