import json

with open('draw_image3.ipynb', 'r', encoding='utf-8') as f:
    nb = json.load(f)

source = ''.join(nb['cells'][20]['source'])

# 查找所有gridspec相关的行
print("Looking for gridspec and ax3=fig.add_subplot:")
lines = source.split('\n')
for i, line in enumerate(lines):
    if 'gridspec.GridSpecFromSubplotSpec' in line and 'gs_down' in line:
        print(f"Line {i}: {line}")
    elif 'ax3=fig.add_subplot(gs_down' in line:
        print(f"Line {i}: {line}")
        
# 验证关键替换
print("\nKey replacements:")
print(f"✓ figsize=(9.5,7): {'figsize=(9.5,7)' in source}")
print(f"✓ gs_down and (1,1): {'GridSpecFromSubplotSpec(1,1,subplot_spec=gs[1])' in source}")
print(f"✓ ax3=fig.add_subplot(gs_down[0]): {'ax3=fig.add_subplot(gs_down[0])' in source}")
