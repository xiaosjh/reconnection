import json

print("\n📊 最终验证 - 修改细节\n")

with open('draw_image3.ipynb', 'r', encoding='utf-8') as f:
    nb = json.load(f)

source = ''.join(nb['cells'][20]['source'])
lines = source.split('\n')

# 显示关键代码片段
print("1️⃣  在图b (ax2)上添加contour的代码:")
print("-" * 70)
for i, line in enumerate(lines[80:92], start=81):
    print(f"{i:3d}: {line}")

print("\n2️⃣  gridspec和ax3的新定义:")
print("-" * 70)
for i, line in enumerate(lines[100:110], start=101):
    print(f"{i:3d}: {line}")

print("\n3️⃣  图表标签和ax_all列表:")
print("-" * 70)
for i, line in enumerate(lines[135:155], start=136):
    print(f"{i:3d}: {line}")

print("\n✅ 验证总结:")
print("-" * 70)
checks = [
    ("figsize调整为(9.5,7)", "figsize=(9.5,7)" in source),
    ("在ax2上添加了contour", "contour=ax2.contour" in source),
    ("level=800的蓝色线", "levels=[800],colors=['blue']" in source),
    ("gridspec改为1列", "GridSpecFromSubplotSpec(1,1,subplot_spec=gs[1])" in source),
    ("ax3变为非投影subplot", "ax3=fig.add_subplot(gs_down[0])" in source),
    ("ax_all只包含ax1和ax2", "ax_all=[ax1,ax2]" in source),
    ("原图d现在是图c", "ax3.text(0.02,0.9,'(c)',color='black'" in source),
    ("没有投影ax3", "projection=jet_region" not in source),
]

for desc, result in checks:
    status = "✓" if result else "✗"
    print(f"{status} {desc}")

print("\n" + "=" * 70)
