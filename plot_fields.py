"""
Vizualizace 3D výsledků: řezy XY (z=střed), XZ (y=střed), YZ (x=za krychlí).
Spuštění: python plot_fields.py
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import TwoSlopeNorm
import sys, os

DAT = "fields.dat"
if not os.path.exists(DAT):
    print(f"'{DAT}' nenalezen – spusť nejdříve ./channelCube")
    sys.exit(1)

data  = np.loadtxt(DAT, comments="#")
i_arr = data[:,0].astype(int); j_arr = data[:,1].astype(int); k_arr = data[:,2].astype(int)
x=data[:,3]; y=data[:,4]; z=data[:,5]
ux=data[:,6]; uy=data[:,7]; uz=data[:,8]
p=data[:,9]; solid=data[:,11].astype(int)

nx=i_arr.max()+1; ny=j_arr.max()+1; nz=k_arr.max()+1
dx=x[1]-x[0]; dy=y[nx]-y[0]; dz=z[nx*ny]-z[0]
Lx=nx*dx; Ly=ny*dy; Lz=nz*dz

spd = np.sqrt(ux**2+uy**2+uz**2)
spd[solid==1] = np.nan
ux_m=np.where(solid==1,np.nan,ux)
uy_m=np.where(solid==1,np.nan,uy)
p_m =np.where(solid==1,np.nan,p)

def get_slice_xy(k_idx):
    """Vrátí 2D pole pro řez z=const (k=k_idx)."""
    mask = (k_arr==k_idx)
    X2 = x[mask].reshape(ny,nx)
    Y2 = y[mask].reshape(ny,nx)
    UX = ux_m[mask].reshape(ny,nx)
    UY = uy_m[mask].reshape(ny,nx)
    P2 = p_m[mask].reshape(ny,nx)
    SP = spd[mask].reshape(ny,nx)
    S2 = solid[mask].reshape(ny,nx)
    return X2, Y2, UX, UY, P2, SP, S2

def get_slice_xz(j_idx):
    mask = (j_arr==j_idx)
    X2 = x[mask].reshape(nz,nx)
    Z2 = z[mask].reshape(nz,nx)
    UX = ux_m[mask].reshape(nz,nx)
    UZ = uz[mask].reshape(nz,nx)
    UZ = np.where(solid[mask].reshape(nz,nx)==1, np.nan, UZ)
    P2 = p_m[mask].reshape(nz,nx)
    SP = spd[mask].reshape(nz,nx)
    S2 = solid[mask].reshape(nz,nx)
    return X2, Z2, UX, UZ, P2, SP, S2

def add_cube_xy(ax, sx, sy, s2d):
    sxf=sx[s2d==1]; syf=sy[s2d==1]
    if len(sxf):
        r=mpatches.Rectangle((sxf.min()-dx/2,syf.min()-dy/2),
                              sxf.max()-sxf.min()+dx,syf.max()-syf.min()+dy,
                              fc='dimgray',ec='black',lw=1.5,zorder=10)
        ax.add_patch(r)

def add_cube_xz(ax, sx, sz, s2d):
    sxf=sx[s2d==1]; szf=sz[s2d==1]
    if len(sxf):
        r=mpatches.Rectangle((sxf.min()-dx/2,szf.min()-dz/2),
                              sxf.max()-sxf.min()+dx,szf.max()-szf.min()+dz,
                              fc='dimgray',ec='black',lw=1.5,zorder=10)
        ax.add_patch(r)

kmid=nz//2; jmid=ny//2
X_xy,Y_xy,UX_xy,UY_xy,P_xy,SP_xy,S_xy = get_slice_xy(kmid)
X_xz,Z_xz,UX_xz,UZ_xz,P_xz,SP_xz,S_xz = get_slice_xz(jmid)

fig, axes = plt.subplots(3,2, figsize=(16,13))

sk=max(1,nx//20)

# ── Rychlost XY ──────────────────────────────────────────────────────────────
ax=axes[0,0]
cm=ax.pcolormesh(X_xy-dx/2,Y_xy-dy/2,SP_xy,cmap='turbo',shading='auto')
plt.colorbar(cm,ax=ax,label="|U| [m/s]")
ax.quiver(X_xy[::2,::sk],Y_xy[::2,::sk],UX_xy[::2,::sk],UY_xy[::2,::sk],
          scale=12,width=0.002,alpha=0.7,color='white')
add_cube_xy(ax,X_xy,Y_xy,S_xy)
ax.set_xlim(0,Lx); ax.set_ylim(0,Ly); ax.set_aspect('equal')
ax.set_title(f"Rychlost |U|  –  řez XY  (z={( kmid+0.5)*dz:.2f})")
ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]")

# ── Rychlost XZ ──────────────────────────────────────────────────────────────
ax=axes[0,1]
cm=ax.pcolormesh(X_xz-dx/2,Z_xz-dz/2,SP_xz,cmap='turbo',shading='auto')
plt.colorbar(cm,ax=ax,label="|U| [m/s]")
sk2=max(1,nz//10)
ax.quiver(X_xz[::sk2,::sk],Z_xz[::sk2,::sk],UX_xz[::sk2,::sk],UZ_xz[::sk2,::sk],
          scale=12,width=0.002,alpha=0.7,color='white')
add_cube_xz(ax,X_xz,Z_xz,S_xz)
ax.set_xlim(0,Lx); ax.set_ylim(0,Lz); ax.set_aspect('equal')
ax.set_title(f"Rychlost |U|  –  řez XZ  (y={( jmid+0.5)*dy:.2f})")
ax.set_xlabel("x [m]"); ax.set_ylabel("z [m]")

# ── Tlak XY ───────────────────────────────────────────────────────────────────
ax=axes[1,0]
pmax=np.nanmax(np.abs(P_xy))+1e-10
cm2=ax.pcolormesh(X_xy-dx/2,Y_xy-dy/2,P_xy,cmap='RdBu_r',
                   norm=TwoSlopeNorm(0,-pmax,pmax),shading='auto')
plt.colorbar(cm2,ax=ax,label="p [Pa]")
add_cube_xy(ax,X_xy,Y_xy,S_xy)
ax.set_xlim(0,Lx); ax.set_ylim(0,Ly); ax.set_aspect('equal')
ax.set_title("Tlak p  –  řez XY"); ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]")

# ── Tlak XZ ───────────────────────────────────────────────────────────────────
ax=axes[1,1]
pmax2=np.nanmax(np.abs(P_xz))+1e-10
cm3=ax.pcolormesh(X_xz-dx/2,Z_xz-dz/2,P_xz,cmap='RdBu_r',
                   norm=TwoSlopeNorm(0,-pmax2,pmax2),shading='auto')
plt.colorbar(cm3,ax=ax,label="p [Pa]")
add_cube_xz(ax,X_xz,Z_xz,S_xz)
ax.set_xlim(0,Lx); ax.set_ylim(0,Lz); ax.set_aspect('equal')
ax.set_title("Tlak p  –  řez XZ"); ax.set_xlabel("x [m]"); ax.set_ylabel("z [m]")

# ── Proudnice XY ──────────────────────────────────────────────────────────────
ax=axes[2,0]
UX_s=np.where(S_xy==1,0.0,UX_xy); UY_s=np.where(S_xy==1,0.0,UY_xy)
cm4=ax.pcolormesh(X_xy-dx/2,Y_xy-dy/2,UX_xy,cmap='coolwarm',shading='auto',
                   vmin=np.nanpercentile(UX_xy,2),vmax=np.nanpercentile(UX_xy,98))
plt.colorbar(cm4,ax=ax,label="u_x [m/s]")
xlin=(np.arange(nx)+0.5)*dx; ylin=(np.arange(ny)+0.5)*dy
try:
    ax.streamplot(xlin,ylin,UX_s,UY_s,density=1.8,linewidth=0.8,
                  color='k',arrowsize=0.8,broken_streamlines=True)
except Exception as e:
    print(f"Streamplot XY: {e}")
add_cube_xy(ax,X_xy,Y_xy,S_xy)
ax.set_xlim(0,Lx); ax.set_ylim(0,Ly); ax.set_aspect('equal')
ax.set_title("Proudnice  –  řez XY"); ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]")

# ── Proudnice XZ ──────────────────────────────────────────────────────────────
ax=axes[2,1]
UX_s2=np.where(S_xz==1,0.0,UX_xz); UZ_s=np.where(S_xz==1,0.0,UZ_xz)
cm5=ax.pcolormesh(X_xz-dx/2,Z_xz-dz/2,UX_xz,cmap='coolwarm',shading='auto',
                   vmin=np.nanpercentile(UX_xz,2),vmax=np.nanpercentile(UX_xz,98))
plt.colorbar(cm5,ax=ax,label="u_x [m/s]")
zlin=(np.arange(nz)+0.5)*dz
try:
    ax.streamplot(xlin,zlin,UX_s2,UZ_s,density=1.8,linewidth=0.8,
                  color='k',arrowsize=0.8,broken_streamlines=True)
except Exception as e:
    print(f"Streamplot XZ: {e}")
add_cube_xz(ax,X_xz,Z_xz,S_xz)
ax.set_xlim(0,Lx); ax.set_ylim(0,Lz); ax.set_aspect('equal')
ax.set_title("Proudnice  –  řez XZ"); ax.set_xlabel("x [m]"); ax.set_ylabel("z [m]")

plt.suptitle(f"3D kanálové proudění kolem krychle   Re=20  ({nx}×{ny}×{nz} buněk)",
             fontsize=13, y=1.01)
plt.tight_layout()
out="flow_visualization.png"
plt.savefig(out,dpi=150,bbox_inches='tight')
print(f"Uloženo: {out}")
plt.show()
