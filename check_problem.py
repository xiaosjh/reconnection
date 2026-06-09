import json

with open('draw_image3.ipynb', 'r', encoding='utf-8') as f:
    nb = json.load(f)

source = ''.join(nb['cells'][20]['source'])

# 查找问题
print("查找问题代码:")
if "ax_all[2].text" in source:
    print("✗ 问题：ax_all[2]仍然存在（ax_all只应该有ax1和ax2）")
    
if "ax4.text(0.02,0.9,'(d)'" in source:
    print("✗ 问题：ax4.text标签仍然存在")
    
if "ax4.text(0.69,0.05" in source:
    print("✗ 问题：ax4的速度文本仍然存在")

# 显示源代码的这一部分
lines = source.split('\n')
print("\n问题行：")
for i, line in enumerate(lines[149:156], start=150):
    print(f"{i}: {line}")
