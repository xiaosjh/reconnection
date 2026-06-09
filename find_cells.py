import json

with open('draw_image3.ipynb', 'r', encoding='utf-8') as f:
    nb = json.load(f)

# 查找包含 ax3 或 ax4 的代码单元
for idx, cell in enumerate(nb['cells']):
    if cell['cell_type'] == 'code':
        source = ''.join(cell['source'])
        if 'ax3' in source or 'ax4' in source or ('contour' in source and 'level' in source):
            print(f"\n{'='*60}\nCell {idx}:\n{'='*60}")
            print(source)
            print(f"\nCell {idx} length: {len(source)}")
