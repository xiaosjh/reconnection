import json

print("\n" + "="*80)
print("✅ 修改完成验证 - draw_image3.ipynb")
print("="*80 + "\n")

with open('draw_image3.ipynb', 'r', encoding='utf-8') as f:
    nb = json.load(f)

source = ''.join(nb['cells'][20]['source'])

# 所有验证
checks = [
    ("图表尺寸 (9.5,7)", "figsize=(9.5,7)" in source),
    ("在图b (ax2)添加contour", "contour=ax2.contour(jet_region.data,levels=[800]" in source),
    ("蓝色线 (blue)", "colors=['blue']" in source),
    ("gridspec为1列", "GridSpecFromSubplotSpec(1,1,subplot_spec=gs[1])" in source),
    ("ax3非投影subplot", "ax3=fig.add_subplot(gs_down[0])" in source),
    ("删除投影ax3", "projection=jet_region" not in source),
    ("ax_all只有[ax1,ax2]", "ax_all=[ax1,ax2]" in source),
    ("原图d变成图c", "ax3.text(0.02,0.9,'(c)'" in source),
    ("移除ax_all[2]", "ax_all[2]" not in source),
    ("移除ax4标签", "ax4.text(0.02,0.9," not in source),
    ("速度文本改为ax3", "transform=ax3.transAxes)" in source and "v_U" in source),
]

print("修改项目:")
print("-" * 80)
for i, (desc, result) in enumerate(checks, 1):
    status = "✅" if result else "❌"
    print(f"{status} {i:2d}. {desc}")

# 统计结果
passed = sum(1 for _, result in checks if result)
total = len(checks)

print("-" * 80)
print(f"\n总体状态: {passed}/{total} 通过\n")

if passed == total:
    print("🎉 所有修改完成且验证通过！")
    print("\n修改摘要:")
    print("  • 删除原图(c)（喷流区域图像）")
    print("  • 在图(b)上添加level=800蓝色等高线（仅喷流区域）")
    print("  • 原图(d)变成新图(c)")
    print("  • 布局调整为2行（顶部两个图，底部一个图）")
else:
    print("⚠️ 有项目未通过，请检查!")

print("\n" + "="*80 + "\n")
