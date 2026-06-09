import json

with open('draw_image3.ipynb', 'r', encoding='utf-8') as f:
    nb = json.load(f)

source = ''.join(nb['cells'][20]['source'])

# 显示关键行
print("gridspec lines:")
for i, line in enumerate(source.split('\n')):
    if 'gridspec.GridSpecFromSubplotSpec' in line and 'gs[1]' in line:
        print(f"  Line: {line}")

print("\nax3 subplot lines:")
for i, line in enumerate(source.split('\n')):
    if 'ax3=fig.add_subplot' in line and 'gs_down' in line:
        print(f"  Line: {line}")
