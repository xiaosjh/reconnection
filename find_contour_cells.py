import json

with open('draw_image3.ipynb', 'r', encoding='utf-8') as f:
    nb = json.load(f)

# 查找包含 contour 和 level 的代码单元
for idx, cell in enumerate(nb['cells']):
    if cell['cell_type'] == 'code':
        source = ''.join(cell['source'])
        if 'contour' in source and 'level' in source:
            print(f"\n{'='*80}\nCell {idx}:\n{'='*80}")
            print(source)
