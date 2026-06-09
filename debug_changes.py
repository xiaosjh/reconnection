import json

with open('draw_image3_new.ipynb', 'r', encoding='utf-8') as f:
    nb = json.load(f)

source = ''.join(nb['cells'][20]['source'])

# 查找gridspec和ax3相关的行
import re

print("Search for gridspec_down:")
matches = re.findall(r'gs_down.*', source)
for m in matches:
    print(f"  {m}")

print("\nSearch for ax3=fig.add_subplot:")
matches = re.findall(r'ax3=fig\.add_subplot\([^)]*\)', source)
for m in matches:
    print(f"  {m}")
