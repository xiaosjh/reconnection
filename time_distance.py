import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.interpolate import interp1d, RegularGridInterpolator
from scipy.ndimage import zoom
import glob
from scipy.interpolate import splprep, splev
data=np.load(r'D:\Learning\PHD1st\magnetic_reconnecion\program\temp.npy')
#——————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
# 生成示例数据 (384*384*34 矩阵)
#这里的34个数据，在时间上对应的是 D:\Learning\PHD1st\magnetic_reconnecion\data\AIA2\DEM10
#里面从32到98的34个时刻的数据
fig, ax = plt.subplots()
ax.imshow(np.log10(data[:, :, 18]),cmap='jet',origin='lower',vmin=6.0,vmax=7.3)  # 参考图，方便画曲线
clicked_points = []

def on_click(event):
    if event.button == 1:  # 左键点击
        x, y = event.xdata, event.ydata  # 记录浮点数坐标
        if x is not None and y is not None:
            clicked_points.append((x, y))
            ax.plot(x, y, 'ro')  # 画出点
            fig.canvas.draw()

def on_key(event):
    if event.key == 'enter':  # 按 Enter 结束
        fig.canvas.mpl_disconnect(cid_click)
        fig.canvas.mpl_disconnect(cid_key)
        plt.close()

cid_click = fig.canvas.mpl_connect('button_press_event', on_click)
cid_key = fig.canvas.mpl_connect('key_press_event', on_key)
plt.show()


#需要注意的是，clicked_points第一个坐标返回的是x
#这里我没有按照x排序，因为那个曲线可能会往回勾
#clicked_points = np.array(sorted(clicked_points, key=lambda p: p[0]))
clicked_points = np.array(clicked_points)
# 计算插值点数：按点击点间距的平均值来决定
total_length = np.sum(np.sqrt(np.diff(clicked_points[:, 0])**2 + np.diff(clicked_points[:, 1])**2))
num_interp = int(total_length / 5)  # 经验值：每 5 个像素插值一个点


#插值会导致最后的曲线弯弯绕绕不太美观，直接拟合则存在x不是单调的问题

from scipy.interpolate import make_smoothing_spline
#我用这个B样条插值，通过参数U的方式，绕过x不是单调的问题，因为x不是单调但u单调就可以了（机智！！！！！！！！！！！！）
x=clicked_points[:,0]
y=clicked_points[:,1]
tck,u=splprep([x,y],s=0)
new_u=np.linspace(0,1,num_interp)
smooth_x,smooth_y=splev(new_u,tck)

#这个数组记录的是插值后的点的坐标
fit_points = np.column_stack((smooth_y, smooth_x))
#interp_points2=interp_points-2

A_x = smooth_x
A_y = smooth_y

# 计算切向量（dx/du, dy/du）
dx, dy = splev(new_u, tck, der=1)

# 计算法向量（旋转 90°）,这个是向里走
norms = np.sqrt(dx**2 + dy**2)  # 计算切向量长度
nx, ny = dy / norms, -dx / norms  # 法向量归一化（向外方向）

# 沿法向量方向收缩 5，10 个像素
shrink_distance = 5
B_y = A_y - shrink_distance * ny
B_x = A_x - shrink_distance * nx
C_y = A_y - 2 * shrink_distance * ny
C_x = A_x - 2 * shrink_distance * nx

fit_points2 = np.column_stack((B_y, B_x))
fit_points3 = np.column_stack((C_y, C_x))


#这里画一下图，看物质流的运动和那条线的关系，方便修正线的位置
for t in range(data.shape[2]):
    plt.imshow(np.log10(data[:, :, t]),cmap='jet',origin='lower',vmin=6.0,vmax=7.3)
    plt.plot(fit_points[:,1],fit_points[:,0],color='red')
    plt.plot(B_x, B_y, color='green')
    plt.plot(C_x, C_y, color='blue')
    #plt.plot(interp_points2[:,1],interp_points2[:,0],color='blue')
    plt.axis('off')
    plt.savefig(r'D:\Learning\PHD1st\magnetic_reconnecion\data\AIA2\DEM_slice\\'+str(t)+'.png')
    plt.close()


y_grid = np.arange(data.shape[0])
x_grid = np.arange(data.shape[1])
# 使用双线性插值获取 20 个时间步上的像素值
result_matrix = np.zeros((num_interp, data.shape[2]))
for t in range(data.shape[2]):
    interp_func = RegularGridInterpolator((y_grid,x_grid), data[:, :, t], method='linear')
    #这个矩阵记录的是插值点对应的值
    result_matrix[:, t] = (interp_func(fit_points)+interp_func(fit_points2)+interp_func(fit_points3))/3
np.savez(r'D:\Learning\PHD1st\magnetic_reconnecion\data\AIA2\fit_points.npz',
         fit_points=fit_points, fit_points2=fit_points2, fit_points3=fit_points3)

np.save(r'D:\Learning\PHD1st\magnetic_reconnecion\data\AIA2\result_matrix.npy',result_matrix)

print("插值后点的数量:", num_interp)
print("最终的结果矩阵形状：", result_matrix.shape)
