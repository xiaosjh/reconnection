import json

with open('draw_image3.ipynb', 'r', encoding='utf-8') as f:
    nb = json.load(f)

# 找到Cell 20并获取其完整源代码
for idx, cell in enumerate(nb['cells']):
    if cell['cell_type'] == 'code':
        source = ''.join(cell['source'])
        if 'contour=ax3.contour' in source:
            print(f"Found target cell at index {idx}")
            print("Full source code:")
            print(source)
            break
