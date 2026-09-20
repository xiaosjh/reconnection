from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Arc

plt.rcParams.update({'font.family': 'Microsoft YaHei', 'font.size': 12,
                     'axes.unicode_minus': False, 'mathtext.fontset': 'dejavusans'})
blue, green, orange, dark = '#1261a0', '#16806a', '#d27712', '#243447'
O = np.array([0., 0.])
V = np.array([3., 0.])
phi = np.deg2rad(30)
P = np.array([np.cos(phi), np.sin(phi)])
T = np.array([1/3, np.sqrt(1-1/9)])
rho = np.degrees(np.arctan2(P[1], V[0]-P[0]))
alpha = np.degrees(np.arcsin(1/3))

fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
fig.subplots_adjust(left=.03, right=.99, bottom=.12, top=.82, wspace=.12)
fig.suptitle('太阳—观测者的二维截面：角度的顶点在哪里？', fontsize=20, y=.97)

def line(ax, a, b, **kwargs):
    ax.plot([a[0], b[0]], [a[1], b[1]], **kwargs)

def base(ax):
    ax.add_patch(Circle(O, 1, facecolor='#fff4dc', edgecolor=orange, lw=2))
    line(ax, O, V, color=dark, lw=1.8)
    ax.scatter([0, 3], [0, 0], s=30, color=dark, zorder=5)
    ax.text(-.13, -.22, r'$O$', fontsize=17)
    ax.text(-.8, -.63, '太阳截面', color=orange)
    ax.text(2.91, -.25, r'$V$', fontsize=17)
    ax.text(3.02, -.49, '观测者', ha='center')
    ax.text(1.55, -.22, r'$OV=D$', ha='center', color=dark)
    ax.set(xlim=(-1.15, 3.55), ylim=(-1.17, 1.35), aspect='equal')
    ax.axis('off')

ax = axes[0]
base(ax)
ax.set_title(r'① 目标点：$\rho$ 在观测者处，$\theta$ 在表面处', fontsize=14, pad=13)
line(ax, O, P, color=dark, lw=2)
line(ax, V, P, color=blue, lw=2.3)
ax.annotate('', xy=1.75*P, xytext=P,
            arrowprops={'arrowstyle':'->', 'color':green, 'lw':2})
ax.scatter(*P, color=dark, s=35, zorder=6)
ax.text(P[0]-.23, P[1]+.16, r'$P$', fontsize=17)
ax.text(.18, .37, r'$OP=R$', fontsize=13)
ax.text(1.57, 1.01, '向外法线', color=green)
ax.text(1.69, .51, '指向目标的视线', color=blue, rotation=-rho, fontsize=11)
ax.add_patch(Arc(V, 1.35, 1.35, theta1=180-rho, theta2=180, color=blue, lw=2.3))
a=np.deg2rad(180-rho/2)
ax.text(*(V+.88*np.array([np.cos(a),np.sin(a)])), r'$\rho$', color=blue, fontsize=20,
        ha='center', va='center')
ax.add_patch(Arc(P, .64, .64, theta1=-rho, theta2=30, color=green, lw=2.3))
a=np.deg2rad((30-rho)/2)
ax.text(*(P+.46*np.array([np.cos(a),np.sin(a)])), r'$\theta$', color=green,
        fontsize=20, ha='center', va='center')

ax=axes[1]
base(ax)
ax.set_title(r'② 圆盘边缘：$\alpha$ 是太阳角半径', fontsize=14, pad=13)
line(ax, O, T, color=dark, lw=2)
line(ax, V, T, color=orange, lw=2.3)
ax.scatter(*T, color=dark, s=35, zorder=6)
ax.text(T[0]-.22, T[1]+.08, r'$T$', fontsize=17)
ax.text(.55, 1.13, '边缘切点', fontsize=12, color=orange)
ax.text(-.24, .48, r'$OT=R$', fontsize=13)
u=(O-T)/np.linalg.norm(O-T)
v=(V-T)/np.linalg.norm(V-T)
q=np.array([T+.13*u, T+.13*(u+v), T+.13*v])
ax.plot(q[:,0],q[:,1],color=dark,lw=1.4)
ax.add_patch(Arc(V, 1.4, 1.4, theta1=180-alpha, theta2=180, color=orange, lw=2.3))
a=np.deg2rad(180-alpha/2)
ax.text(*(V+.92*np.array([np.cos(a),np.sin(a)])), r'$\alpha$', color=orange,
        fontsize=20, ha='center', va='center')
ax.text(1.64, .77, r'$OT\perp VT$', color=dark, fontsize=16)

fig.text(.5, .035, '示意图未按真实日地距离比例绘制。圆盘中心方向是 V → O；P 是太阳表面上的目标点。',
         ha='center', color='#52606d', fontsize=12)
out=Path(__file__).resolve().parent
fig.savefig(out/'rho_geometry_explanation.png', dpi=180, facecolor='white')
fig.savefig(out/'rho_geometry_explanation.svg', facecolor='white')
plt.close(fig)
print(out/'rho_geometry_explanation.png')
