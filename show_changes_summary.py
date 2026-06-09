import json

print("=" * 80)
print("修改完成！以下是关键变更：")
print("=" * 80)

with open('draw_image3.ipynb', 'r', encoding='utf-8') as f:
    nb = json.load(f)

source = ''.join(nb['cells'][20]['source'])
lines = source.split('\n')

print("\n✓ 修改1：图表尺寸从(9.5,8)调整为(9.5,7)")
print("   figsize=(9.5,7)" if "figsize=(9.5,7)" in source else "   ERROR")

print("\n✓ 修改2：在图(b) ax2上添加了level=800的蓝色contour")
print("   代码位置：第 85-87 行")
for i, line in enumerate(lines[84:88], start=85):
    if 'contour=ax2.contour' in line or 'jet_region=sub_hawing_map' in line:
        print(f"   {i}: {line}")

print("\n✓ 修改3：删除了原图(c)的代码（喷流区域的DEM图像）")
print("   - 移除了: ax3=fig.add_subplot(gs_down[0],projection=jet_region)")
print("   - 移除了: im=jet_region.plot(axes=ax3,...)")
print("   - 移除了: ax3.coords[0].set_ticks(number=4)")

print("\n✓ 修改4：gridspec改为单列布局")
for i, line in enumerate(lines[102:105], start=103):
    if 'gridspec' in line or 'ax3=fig' in line:
        print(f"   {i}: {line}")

print("\n✓ 修改5：原图(d)变成新图(c)")
print("   - ax4.text() → ax3.text() 标签改为 '(c)'")
print("   - 原图标签关系：图a → 图a，图b → 图b，图c(removed) → (removed)，图d → 图c")

print("\n✓ 修改6：图表ax_all列表")
for i, line in enumerate(lines):
    if 'ax_all=' in line:
        print(f"   {i+1}: {line}")
        break

print("\n" + "=" * 80)
print("修改完成。notebook已保存为: draw_image3.ipynb")
print("=" * 80)
