import json

# 读取备份
with open('draw_image3.ipynb.bak', 'r', encoding='utf-8') as f:
    nb = json.load(f)

# 获取原始source list（保持原始格式）
original_source_list = nb['cells'][20]['source']

# 将source list转换为单个字符串，处理逗号和转义
if isinstance(original_source_list, list):
    cell_text = ''.join(original_source_list)
else:
    cell_text = original_source_list

print(f"Original is list: {isinstance(original_source_list, list)}")
print(f"First 5 items: {original_source_list[:5]}")
print(f"Total items: {len(original_source_list)}")
