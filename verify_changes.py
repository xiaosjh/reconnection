import json

# 验证修改
with open('draw_image3_new.ipynb', 'r', encoding='utf-8') as f:
    nb = json.load(f)

# 检查Cell 20
cell = nb['cells'][20]
source = ''.join(cell['source'])

# 检查关键代码
checks = [
    ('figsize=(9.5,7)', 'figsize改正'),
    ('jet_region=sub_hawing_map.submap(point4,top_right=point5)', 'jet_region定义'),
    ('contour=ax2.contour(jet_region.data,levels=[800]', 'contour在ax2上'),
    ('gs_down=gridspec.GridSpecFromSubplotSpec(1,1,subplot_spec=gs[1])', 'gridspec改为1列'),
    ('ax3=fig.add_subplot(gs_down[0])', 'ax3变为非投影subplot'),
    ('ax_all=[ax1,ax2]', 'ax_all移除了ax3'),
    ('ax3.text(0.02,0.9,\'(c)\'', '标签改为(c)'),
]

print("验证修改:")
for check_str, desc in checks:
    if check_str in source:
        print(f"✓ {desc}")
    else:
        print(f"✗ {desc} - NOT FOUND")
        
# 检查不应该出现的代码
unwanted = [
    ('ax3=fig.add_subplot(gs_down[0],projection=jet_region)', '不应该有投影ax3'),
    ('ax3.coords[0].set_ticks(number=4)', '不应该有ax3坐标设置'),
    ('ax4.text(0.02,0.9,\'(d)\'', '不应该有ax4标签'),
]

print("\n检查不应该出现的代码:")
for check_str, desc in unwanted:
    if check_str not in source:
        print(f"✓ {desc}")
    else:
        print(f"✗ {desc} - STILL PRESENT")
