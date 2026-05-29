"""
Vizualizace 3D výpočetní sítě: průřezy XY a XZ, 3D schéma.
Spuštění: python plot_mesh.py
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import sys, os

DAT = "mesh.dat"
if not os.path.exists(DAT):
    print(f"'{DAT}' nenalezen – spusť nejdříve ./channelCube")
    sys.exit(1)

data  = np.loadtxt(DAT, comments="#")
i_arr = data[:,0].astype(int)
j_arr = data[:,1].astype(int)
k_arr = data[:,2].astype(int)
x     = data[:,3]; y = data[:,4]; z = data[:,5]
solid = data[:,6].astype(int)

nx = i_arr.max()+1; ny = j_arr.max()+1; nz = k_arr.max()+1
dx = x[1]-x[0]; dy = y[nx]-y[0]; dz = z[nx*ny]-z[0]
Lx = nx*dx; Ly = ny*dy; Lz = nz*dz

print(f"3D síť: {nx}×{ny}×{nz}  dx={dx:.4f}  dy={dy:.4f}  dz={dz:.4f}")
print(f"Pevných buněk: {solid.sum()}")

# Pozice pevné krychle
sx = x[solid==1]; sy = y[solid==1]; sz = z[solid==1]
cx0,cx1 = sx.min()-dx/2, sx.max()+dx/2
cy0,cy1 = sy.min()-dy/2, sy.max()+dy/2
cz0,cz1 = sz.min()-dz/2, sz.max()+dz/2

fig = plt.figure(figsize=(16,10))

# ── 1. 3D schéma domény ──────────────────────────────────────────────────────
ax3d = fig.add_subplot(1,2,1, projection='3d')

def draw_box(ax, x0,x1,y0,y1,z0,z1, color, alpha=0.15, lw=1):
    verts = [
        [(x0,y0,z0),(x1,y0,z0),(x1,y1,z0),(x0,y1,z0)],
        [(x0,y0,z1),(x1,y0,z1),(x1,y1,z1),(x0,y1,z1)],
        [(x0,y0,z0),(x0,y1,z0),(x0,y1,z1),(x0,y0,z1)],
        [(x1,y0,z0),(x1,y1,z0),(x1,y1,z1),(x1,y0,z1)],
        [(x0,y0,z0),(x1,y0,z0),(x1,y0,z1),(x0,y0,z1)],
        [(x0,y1,z0),(x1,y1,z0),(x1,y1,z1),(x0,y1,z1)],
    ]
    pc = Poly3DCollection(verts, alpha=alpha, facecolor=color, edgecolor='k', lw=lw)
    ax.add_collection3d(pc)

draw_box(ax3d, 0,Lx, 0,Ly, 0,Lz, 'lightcyan', alpha=0.08, lw=0.5)
draw_box(ax3d, cx0,cx1, cy0,cy1, cz0,cz1, 'steelblue', alpha=0.7, lw=1.5)

# Průřez XY (z=Lz/2, červený rám)
zmid = Lz/2
rect_verts = [[(0,0,zmid),(Lx,0,zmid),(Lx,Ly,zmid),(0,Ly,zmid)]]
pc2 = Poly3DCollection(rect_verts, alpha=0.2, facecolor='red', edgecolor='red', lw=1.5)
ax3d.add_collection3d(pc2)

# Průřez XZ (y=Ly/2, zelený rám)
ymid = Ly/2
rect_verts2 = [[(0,ymid,0),(Lx,ymid,0),(Lx,ymid,Lz),(0,ymid,Lz)]]
pc3 = Poly3DCollection(rect_verts2, alpha=0.2, facecolor='green', edgecolor='green', lw=1.5)
ax3d.add_collection3d(pc3)

# Šipka toku
ax3d.quiver(-0.3,Ly/2,Lz/2, 0.5,0,0, color='blue', lw=2, arrow_length_ratio=0.3)
ax3d.text(-0.4,Ly/2,Lz/2*1.4, "U_in=1 m/s", color='blue', fontsize=9)

ax3d.set_xlim(-0.2,Lx+0.1); ax3d.set_ylim(-0.1,Ly+0.1); ax3d.set_zlim(-0.1,Lz+0.1)
ax3d.set_xlabel("x [m]"); ax3d.set_ylabel("y [m]"); ax3d.set_zlabel("z [m]")
ax3d.set_title(f"3D tunel  {nx}×{ny}×{nz} buněk\n(modrá=krychle, červený řez=XY, zelený=XZ)")
ax3d.view_init(elev=20, azim=-50)

# ── 2. Průřezy sítě XY a XZ ──────────────────────────────────────────────────
ax2 = fig.add_subplot(2,2,2)
ax3 = fig.add_subplot(2,2,4)

# XY řez (k = nz//2)
kmid = nz//2
mask_xy = (k_arr == kmid)
xi = x[mask_xy]; yi = y[mask_xy]; si = solid[mask_xy]

# Mřížkové čáry
step = max(1, nx//20)
for jj in range(0, ny+1, step):
    ax2.axhline(jj*dy, color='lightgray', lw=0.4)
for ii in range(0, nx+1, step):
    ax2.axvline(ii*dx, color='lightgray', lw=0.4)

# Pevné buňky
sxi = xi[si==1]; syi = yi[si==1]
if len(sxi):
    rect = mpatches.Rectangle((sxi.min()-dx/2, syi.min()-dy/2),
                               sxi.max()-sxi.min()+dx, syi.max()-syi.min()+dy,
                               fc='steelblue', ec='navy', lw=1.5)
    ax2.add_patch(rect)

ax2.add_patch(mpatches.Rectangle((0,0),Lx,Ly,fc='none',ec='black',lw=2))
ax2.set_xlim(-0.1,Lx+0.1); ax2.set_ylim(-0.1,Ly+0.1)
ax2.set_aspect('equal')
ax2.set_xlabel("x [m]"); ax2.set_ylabel("y [m]")
ax2.set_title(f"Řez XY  (z = {(kmid+0.5)*dz:.2f} m)", color='red')
ax2.text(0.02,0.95,f"Inlet → U=1\nOutlet → p=0",
         transform=ax2.transAxes, fontsize=8, va='top',
         bbox=dict(fc='lightyellow',ec='gray'))

# XZ řez (j = ny//2)
jmid = ny//2
mask_xz = (j_arr == jmid)
xi2 = x[mask_xz]; zi2 = z[mask_xz]; si2 = solid[mask_xz]

stepz = max(1, nz//20)
for kk in range(0, nz+1, stepz):
    ax3.axhline(kk*dz, color='lightgray', lw=0.4)
for ii in range(0, nx+1, step):
    ax3.axvline(ii*dx, color='lightgray', lw=0.4)

sxi2 = xi2[si2==1]; szi2 = zi2[si2==1]
if len(sxi2):
    rect2 = mpatches.Rectangle((sxi2.min()-dx/2, szi2.min()-dz/2),
                                sxi2.max()-sxi2.min()+dx, szi2.max()-szi2.min()+dz,
                                fc='steelblue', ec='navy', lw=1.5)
    ax3.add_patch(rect2)

ax3.add_patch(mpatches.Rectangle((0,0),Lx,Lz,fc='none',ec='black',lw=2))
ax3.set_xlim(-0.1,Lx+0.1); ax3.set_ylim(-0.1,Lz+0.1)
ax3.set_aspect('equal')
ax3.set_xlabel("x [m]"); ax3.set_ylabel("z [m]")
ax3.set_title(f"Řez XZ  (y = {(jmid+0.5)*dy:.2f} m)", color='green')

plt.suptitle(f"Výpočetní síť 3D tunelu  {nx}×{ny}×{nz}  Re=20", fontsize=13)
plt.tight_layout()
out = "mesh_visualization.png"
plt.savefig(out, dpi=150, bbox_inches='tight')
print(f"Uloženo: {out}")
plt.show()
