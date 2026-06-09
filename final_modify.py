import json

# 读取备份
with open('draw_image3.ipynb.bak', 'r', encoding='utf-8') as f:
    nb = json.load(f)

# 获取原始source list
original_source_list = nb['cells'][20]['source']

# 将list转为单个字符串便于处理
cell_text = ''.join(original_source_list)

# 执行替换
# 1. 改变figsize
cell_text = cell_text.replace('figsize=(9.5,8)', 'figsize=(9.5,7)')

# 2. 查找需要添加contour的位置，然后在那儿加上代码
old_section = """sub_hawing_map.draw_quadrangle(point4,top_right=point5,axes=ax2,edgecolor='white')
sub_aia.draw_quadrangle(point4,top_right=point5,axes=ax1,edgecolor='white')
#------------------------------"""

new_section = """sub_hawing_map.draw_quadrangle(point4,top_right=point5,axes=ax2,edgecolor='white')
sub_aia.draw_quadrangle(point4,top_right=point5,axes=ax1,edgecolor='white')

#在ax2上添加喷流区域的contour（level=800）
jet_region=sub_hawing_map.submap(point4,top_right=point5)
contour=ax2.contour(jet_region.data,levels=[800],colors=['blue'],transform=ax2.get_transform(jet_region.wcs))
#------------------------------"""

cell_text = cell_text.replace(old_section, new_section)

# 3. 删除原始的图3代码段并替换
old_section3 = """----------------------------------------------------------------------------------------------
#图3
gs_down=gridspec.GridSpecFromSubplotSpec(1,2,subplot_spec=gs[1],width_ratios=[0.6,1])
jet_region=sub_hawing_map.submap(point4,top_right=point5)
ax3=fig.add_subplot(gs_down[0],projection=jet_region)
im=jet_region.plot(axes=ax3,cmap='afmhot',vmin=0, vmax = 2 * sub_hawing_map.data.mean())
ax3.coords[0].set_ticks(number=4)

contour=ax3.contour(jet_region.data,levels=[800],colors=['blue'],transform=ax3.get_transform(jet_region.wcs))

#-------------------------------------------------------------------------------------------------
#图4
ax4=fig.add_subplot(gs_down[1])"""

new_section3 = """----------------------------------------------------------------------------------------------
#图3（原图4，现在变成图c）
gs_down=gridspec.GridSpecFromSubplotSpec(1,1,subplot_spec=gs[1])
ax3=fig.add_subplot(gs_down[0])"""

cell_text = cell_text.replace(old_section3, new_section3)

# 4. 替换ax_all列表
cell_text = cell_text.replace('ax_all=[ax1,ax2,ax3]', 'ax_all=[ax1,ax2]')

# 5. 替换底部的ax4标签和代码为ax3
# 删除中间的ax3循环部分里面关于ax3的坐标处理
cell_text = cell_text.replace(
    """ax_all[0].text(0.02,0.93,'(a)  AIA 131  21:17:06 UT',color='white',transform=ax_all[0].transAxes,fontsize=12)
ax_all[1].text(0.02,0.85,r'(b)  CHASE H$\\mathrm{\\alpha}$ 6561.55-6562.03$\\mathrm{\\AA}$'+'\n       21:17:11 UT',color='white',transform=ax_all[1].transAxes,fontsize=12)
ax_all[2].text(0.02,0.93,'(c)',color='white',transform=ax_all[2].transAxes,fontsize=12)
ax4.text(0.02,0.9,'(d)',color='black',transform=ax4.transAxes,fontsize=12)
ax4.text(0.69,0.05,r'$\\mathrm{v_U}$'+'=-15.65 '+r'$\\mathrm{km \\ s^{-1}}$'+'\\n'+r'$\\mathrm{v_L} $'+'= 20.78 '+r'$\\mathrm{km \\ s^{-1}}$',color='black',
         fontsize=10,transform=ax4.transAxes)""",
    """ax_all[0].text(0.02,0.93,'(a)  AIA 131  21:17:06 UT',color='white',transform=ax_all[0].transAxes,fontsize=12)
ax_all[1].text(0.02,0.85,r'(b)  CHASE H$\\mathrm{\\alpha}$ 6561.55-6562.03$\\mathrm{\\AA}$'+'\n       21:17:11 UT',color='white',transform=ax_all[1].transAxes,fontsize=12)
ax3.text(0.02,0.9,'(c)',color='black',transform=ax3.transAxes,fontsize=12)
ax3.text(0.69,0.05,r'$\\mathrm{v_U}$'+'=-15.65 '+r'$\\mathrm{km \\ s^{-1}}$'+'\\n'+r'$\\mathrm{v_L} $'+'= 20.78 '+r'$\\mathrm{km \\ s^{-1}}$',color='black',
         fontsize=10,transform=ax3.transAxes)"""
)

# 6. 替换其他ax4的引用为ax3
cell_text = cell_text.replace('ax4.errorbar(lam,I,', 'ax3.errorbar(lam,I,')
cell_text = cell_text.replace('ax4.plot(lam, y_fit', 'ax3.plot(lam, y_fit')
cell_text = cell_text.replace('ax4.axvline(x=6562.8', 'ax3.axvline(x=6562.8')
cell_text = cell_text.replace('ax4.set_xlabel(r', 'ax3.set_xlabel(r')
cell_text = cell_text.replace('ax4.set_ylabel(', 'ax3.set_ylabel(')
cell_text = cell_text.replace('ax4.legend(', 'ax3.legend(')

# 转回source list格式（保持每一行的'\n'）
# 按照原始的格式，每一行都应该以'\n'结尾
new_source_list = []
for line in cell_text.split('\n')[:-1]:  # 排除最后空行
    new_source_list.append(line + '\n')
# 最后一行不加'\n'
if cell_text.split('\n')[-1]:
    new_source_list.append(cell_text.split('\n')[-1])

# 更新cell
nb['cells'][20]['source'] = new_source_list

# 保存
with open('draw_image3.ipynb', 'w', encoding='utf-8') as f:
    json.dump(nb, f, ensure_ascii=False, indent=1)

print(f"修改成功!")
print(f"New source list length: {len(new_source_list)}")
print(f"First 10 items:\n")
for i, item in enumerate(new_source_list[:10]):
    print(f"  [{i}]: {repr(item)}")
