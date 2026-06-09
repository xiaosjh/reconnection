import json
import os

# 备份原文件
backup_path = 'draw_image3.ipynb.bak'
if not os.path.exists(backup_path):
    os.system(f'copy draw_image3.ipynb {backup_path}')

# 读取notebook
with open(backup_path, 'r', encoding='utf-8') as f:
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
    print(f"Found target cell at index {target_cell_idx}")
    
    # 获取原始source并分析结构
    original_source = nb['cells'][target_cell_idx]['source']
    if isinstance(original_source, list):
        original_text = ''.join(original_source)
    else:
        original_text = original_source
    
    print(f"Original source length: {len(original_text)} characters")
    
    # 执行替换操作：
    # 1. 改变figsize从(9.5,8)到(9.5,7)
    modified_text = original_text.replace('figsize=(9.5,8)', 'figsize=(9.5,7)')
    
    # 2. 在ax2上添加contour（在sub_aia.draw_quadrangle之后添加）
    contour_code = """
#在ax2上添加喷流区域的contour（level=800）
jet_region=sub_hawing_map.submap(point4,top_right=point5)
contour=ax2.contour(jet_region.data,levels=[800],colors=['blue'],transform=ax2.get_transform(jet_region.wcs))
"""
    
    # 找到"sub_aia.draw_quadrangle"的第二个出现位置后插入
    insert_pos = modified_text.find('sub_aia.draw_quadrangle(point4,top_right=point5,axes=ax1,edgecolor=\'white\')')
    if insert_pos != -1:
        insert_pos = modified_text.find('\n', insert_pos) + 1
        modified_text = modified_text[:insert_pos] + contour_code + modified_text[insert_pos:]
    
    # 3. 删除原来的图3代码块（从"#-----"到"ax3.coords[0].set_ticks"的部分）
    # 查找需要删除的部分
    start_delete = modified_text.find('#----------------------------------------------------------------------------------------------\n#图3')
    end_delete = modified_text.find('contour=ax3.contour(jet_region.data,levels=[800]')
    if start_delete != -1 and end_delete != -1:
        # 找到这行的末尾
        end_delete = modified_text.find(')\n', end_delete) + 2
        # 删除
        deleted_section = modified_text[start_delete:end_delete]
        modified_text = modified_text[:start_delete] + modified_text[end_delete:]
        print(f"Deleted original figure C section ({end_delete - start_delete} chars)")
    
    # 4. 改变gridspec_down从(1,2)到(1,1)
    modified_text = modified_text.replace(
        'gs_down=gridspec.GridSpecFromSubplotSpec(1,2,subplot_spec=gs[1],width_ratios=[0.6,1])',
        'gs_down=gridspec.GridSpecFromSubplotSpec(1,1,subplot_spec=gs[1])'
    )
    
    # 5. 改变ax3的subplot定义
    modified_text = modified_text.replace(
        'ax3=fig.add_subplot(gs_down[0],projection=jet_region)',
        'ax3=fig.add_subplot(gs_down[0])'
    )
    
    # 6. 删除"im=jet_region.plot"和"ax3.coords"那几行
    lines = modified_text.split('\n')
    new_lines = []
    skip_next_lines = False
    for i, line in enumerate(lines):
        if 'ax3=fig.add_subplot(gs_down[0])' in line:
            new_lines.append(line)
        elif 'im=jet_region.plot(axes=ax3' in line:
            continue  # skip this line
        elif 'ax3.coords[0].set_ticks(number=4)' in line:
            continue  # skip this line
        else:
            new_lines.append(line)
    
    modified_text = '\n'.join(new_lines)
    
    # 7. 更新ax_all列表（删除ax3）
    modified_text = modified_text.replace(
        'ax_all=[ax1,ax2,ax3]',
        'ax_all=[ax1,ax2]'
    )
    
    # 8. 更新标签（ax3变成(c)，ax4变成(c)）
    # 先改ax3的标签（这个应该已经被删除了）
    # 改ax4的标签为ax3，并改标签为(c)
    modified_text = modified_text.replace(
        "ax4.text(0.02,0.9,'(d)',color='black'",
        "ax3.text(0.02,0.9,'(c)',color='black'"
    )
    modified_text = modified_text.replace(
        'ax4.text(0.69,0.05,r\'$\\mathrm{v_U}$',
        'ax3.text(0.69,0.05,r\'$\\mathrm{v_U}$'
    )
    modified_text = modified_text.replace(
        'ax4.legend(bbox_to_anchor=(0.45, 0.77))',
        'ax3.legend(bbox_to_anchor=(0.45, 0.77))'
    )
    modified_text = modified_text.replace(
        'ax4.errorbar(lam,I,',
        'ax3.errorbar(lam,I,'
    )
    modified_text = modified_text.replace(
        'ax4.plot(lam, y_fit',
        'ax3.plot(lam, y_fit'
    )
    modified_text = modified_text.replace(
        'ax4.axvline(x=6562.8',
        'ax3.axvline(x=6562.8'
    )
    modified_text = modified_text.replace(
        'ax4.set_xlabel(r\'Wavelength',
        'ax3.set_xlabel(r\'Wavelength'
    )
    modified_text = modified_text.replace(
        'ax4.set_ylabel(\'Intensity (Count)\')',
        'ax3.set_ylabel(\'Intensity (Count)\')'
    )
    
    # 转换为列表格式
    modified_source_list = modified_text.split('\n')
    
    # 更新cell的source
    nb['cells'][target_cell_idx]['source'] = modified_source_list
    
    # 写入新的notebook
    # 需要处理文件锁的问题
    temp_path = 'draw_image3_new.ipynb'
    with open(temp_path, 'w', encoding='utf-8') as f:
        json.dump(nb, f, ensure_ascii=False, indent=1)
    
    print(f"New notebook saved to {temp_path}")
    print(f"Modified source length: {len(modified_text)} characters")
    
else:
    print("Target cell not found!")
