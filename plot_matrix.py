"""
Vizualizace FVM matic: tlaková (pEqn) a hybnostní (UEqn).
Spuštění: python plot_matrix.py
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Rectangle
import scipy.sparse as sp
import sys, os

def load_mat(path, ncols_coef):
    d = np.loadtxt(path, comments="#")
    return d

if not (os.path.exists("pmatrix.dat") and os.path.exists("umatrix.dat")):
    print("pmatrix.dat / umatrix.dat nenalezeny"); sys.exit(1)

# ── Načtení ───────────────────────────────────────────────────────────────────
# pmatrix: i j k x y z aP aE aW aN aS aF aB src solid
pm = np.loadtxt("pmatrix.dat", comments="#")
# umatrix: i j k x y z aP aE aW aN aS aF aB src_x src_y src_z solid
um = np.loadtxt("umatrix.dat", comments="#")

ia = pm[:,0].astype(int); ja = pm[:,1].astype(int); ka = pm[:,2].astype(int)
nx = ia.max()+1; ny = ja.max()+1; nz = ka.max()+1
dx = pm[1,3]-pm[0,3]
dy = pm[nx,4]-pm[0,4]
dz = pm[nx*ny,5]-pm[0,5]

solid = pm[:,-1].astype(int)

# Sloupce koeficientů
def coefs(d):
    return d[:,6],d[:,7],d[:,8],d[:,9],d[:,10],d[:,11],d[:,12]

paP,paE,paW,paN,paS,paF,paB = coefs(pm)
uaP,uaE,uaW,uaN,uaS,uaF,uaB = coefs(um)

print(f"Síť {nx}×{ny}×{nz}  ({nx*ny*nz} buněk)")

# ─────────────────────────────────────────────────────────────────────────────
# Pomocná funkce: rekonstrukce řídké matice pro 2D řez (k = kmid)
# ─────────────────────────────────────────────────────────────────────────────
def slice_sparse(aP_g, aE_g, aW_g, aN_g, aS_g, k_slice):
    """Sestaví řídkou matici pro XY řez (k=k_slice)."""
    N2 = nx * ny
    rows, cols, vals = [], [], []
    base = k_slice * N2   # globální index prvního prvku v řezu

    for jj in range(ny):
        for ii in range(nx):
            c_loc = ii + jj*nx          # lokální index v řezu
            c_glo = base + c_loc        # globální index

            # Diagonála
            rows.append(c_loc); cols.append(c_loc); vals.append(aP_g[c_glo])

            # Východ (i+1)
            if ii < nx-1:
                e_loc = c_loc+1
                rows.append(c_loc); cols.append(e_loc); vals.append(aE_g[c_glo])
            # Západ (i-1)
            if ii > 0:
                w_loc = c_loc-1
                rows.append(c_loc); cols.append(w_loc); vals.append(aW_g[c_glo])
            # Sever (j+1)
            if jj < ny-1:
                n_loc = c_loc+nx
                rows.append(c_loc); cols.append(n_loc); vals.append(aN_g[c_glo])
            # Jih (j-1)
            if jj > 0:
                s_loc = c_loc-nx
                rows.append(c_loc); cols.append(s_loc); vals.append(aS_g[c_glo])

    return sp.csr_matrix((vals, (rows, cols)), shape=(N2, N2))

kmid = nz // 2
Sp_p = slice_sparse(paP, paE, paW, paN, paS, kmid)
Sp_u = slice_sparse(uaP, uaE, uaW, uaN, uaS, kmid)

# ─────────────────────────────────────────────────────────────────────────────
# Pomocná funkce: 2D pole pro XY řez
# ─────────────────────────────────────────────────────────────────────────────
def xy_slice(arr, k_slice):
    mask = (ka == k_slice)
    return arr[mask].reshape(ny, nx)

# ═════════════════════════════════════════════════════════════════════════════
# FIGURA 1 – Sparsity pattern 2D řezu (pEqn vlevo, UEqn vpravo)
# ═════════════════════════════════════════════════════════════════════════════
fig1, axes = plt.subplots(1, 2, figsize=(16, 7))
fig1.suptitle(f"Sparsity pattern – XY řez (k={kmid})  [{nx*ny}×{nx*ny} podmatice]",
              fontsize=13)

for ax, Sp, title in [
    (axes[0], Sp_p, "Tlaková matice pEqn\n(symetrická Laplaceova)"),
    (axes[1], Sp_u, "Hybnostní matice UEqn\n(nesymetrická – upwind konvekce)"),
]:
    ax.spy(Sp, markersize=0.3, color='steelblue', aspect='auto')

    # Zvýrazni identity řádky (pevné buňky) červeně
    sol_2d = xy_slice(solid, kmid).flatten()
    solid_idx = np.where(sol_2d == 1)[0]
    ax.hlines(solid_idx, -0.5, nx*ny-0.5, colors='red', lw=0.8, alpha=0.7)
    ax.vlines(solid_idx, -0.5, nx*ny-0.5, colors='red', lw=0.8, alpha=0.7)

    # Popisky pásů
    ax.axvline(nx, color='orange', lw=0.6, ls='--', alpha=0.6,
               label=f'±{nx} (N/S sousedé)')
    ax.axvline(nx*ny-nx, color='orange', lw=0.6, ls='--', alpha=0.6)
    ax.legend(fontsize=8, loc='lower right')
    ax.set_title(title + f"\n{Sp.nnz} nenulových prvků  ({Sp.nnz/(nx*ny):.1f}/řádek)",
                 fontsize=10)
    ax.set_xlabel("sloupec (buňka)"); ax.set_ylabel("řádek (buňka)")

plt.tight_layout()
plt.savefig("matrix_sparsity.png", dpi=150, bbox_inches='tight')
print("Uloženo: matrix_sparsity.png")

# ═════════════════════════════════════════════════════════════════════════════
# FIGURA 2 – Prostorové rozložení koeficientů (XY řez)
# ═════════════════════════════════════════════════════════════════════════════
fig2, axes = plt.subplots(3, 4, figsize=(18, 11))
fig2.suptitle(f"Koeficienty matic na XY řezu (k={kmid})", fontsize=13)

def plot_coef(ax, data_2d, title, cmap='RdBu', symm=True):
    sx = xy_slice(solid, kmid)
    d  = np.where(sx == 1, np.nan, data_2d)
    v  = np.nanmax(np.abs(d)) if symm else None
    kw = dict(cmap=cmap, shading='auto')
    if symm and v:
        kw['vmin']=-v; kw['vmax']=v
    xg = (np.arange(nx)+0.5)*dx
    yg = (np.arange(ny)+0.5)*dy
    X,Y = np.meshgrid(xg, yg)
    im = ax.pcolormesh(X, Y, d, **kw)
    plt.colorbar(im, ax=ax, shrink=0.9)
    # Překážka
    sol2 = xy_slice(solid, kmid)
    sx_c = xg[sol2.any(axis=0)]
    sy_c = yg[sol2.any(axis=1)]
    if len(sx_c) and len(sy_c):
        ax.add_patch(Rectangle((sx_c.min()-dx/2, sy_c.min()-dy/2),
                                sx_c.max()-sx_c.min()+dx,
                                sy_c.max()-sy_c.min()+dy,
                                fc='dimgray', ec='k', lw=1.5, zorder=5))
    ax.set_aspect('equal'); ax.set_title(title, fontsize=9)
    ax.set_xlabel("x"); ax.set_ylabel("y")

labels_coef = ["aP (diagonal.)", "aE (východ)", "aN (sever)", "aF (přední z+)"]
p_fields = [xy_slice(paP,kmid), xy_slice(paE,kmid),
            xy_slice(paN,kmid), xy_slice(paF,kmid)]
u_fields = [xy_slice(uaP,kmid), xy_slice(uaE,kmid),
            xy_slice(uaN,kmid), xy_slice(uaF,kmid)]

for col, (lbl, pf, uf) in enumerate(zip(labels_coef, p_fields, u_fields)):
    plot_coef(axes[0,col], pf, f"pEqn – {lbl}")
    plot_coef(axes[1,col], uf, f"UEqn – {lbl}")

# Řádek 3: Diagonální dominance obou matic
def diag_dom(aP_g, aE_g, aW_g, aN_g, aS_g, aF_g, aB_g, k_slice):
    aP2 = xy_slice(aP_g, k_slice)
    off = (np.abs(xy_slice(aE_g,k_slice)) + np.abs(xy_slice(aW_g,k_slice)) +
           np.abs(xy_slice(aN_g,k_slice)) + np.abs(xy_slice(aS_g,k_slice)) +
           np.abs(xy_slice(aF_g,k_slice)) + np.abs(xy_slice(aB_g,k_slice)))
    return np.abs(aP2) / (off + 1e-30)

dd_p = diag_dom(paP,paE,paW,paN,paS,paF,paB, kmid)
dd_u = diag_dom(uaP,uaE,uaW,uaN,uaS,uaF,uaB, kmid)

plot_coef(axes[2,0], dd_p, "pEqn – diag. dominance |aP|/Σ|off|", cmap='viridis', symm=False)
plot_coef(axes[2,1], dd_u, "UEqn – diag. dominance |aP|/Σ|off|", cmap='viridis', symm=False)
axes[2,2].axis('off'); axes[2,3].axis('off')

plt.tight_layout()
plt.savefig("matrix_coefficients.png", dpi=150, bbox_inches='tight')
print("Uloženo: matrix_coefficients.png")

# ═════════════════════════════════════════════════════════════════════════════
# FIGURA 3 – Hustý blok kolem krychle (přímé hodnoty stencilu)
# ═════════════════════════════════════════════════════════════════════════════
# Vezmi oblast kolem přední stěny krychle: i ∈ [i_cube-3, i_cube+5], j ∈ [j_cube-2, j_cube+5]
sol2d = xy_slice(solid, kmid)
ci = np.where(sol2d.any(axis=0))[0]
cj = np.where(sol2d.any(axis=1))[0]

if len(ci) and len(cj):
    i0 = max(0, ci.min()-3); i1 = min(nx, ci.max()+4)
    j0 = max(0, cj.min()-2); j1 = min(ny, cj.max()+3)

    N2 = nx*ny
    base = kmid * N2
    rows_sub = [jj*nx+ii for jj in range(j0,j1) for ii in range(i0,i1)]
    nr = j1-j0; nc_b = i1-i0
    dense_p = np.zeros((nr*nc_b, nr*nc_b))
    dense_u = np.zeros_like(dense_p)

    Sp_p2 = Sp_p.toarray()
    Sp_u2 = Sp_u.toarray()
    for ri, r in enumerate(rows_sub):
        for ci2, c2 in enumerate(rows_sub):
            dense_p[ri, ci2] = Sp_p2[r, c2]
            dense_u[ri, ci2] = Sp_u2[r, c2]

    fig3, axes3 = plt.subplots(1, 2, figsize=(16, 7))
    fig3.suptitle(f"Hustý blok matic kolem přední stěny krychle\n"
                  f"Řádky/sloupce = buňky i∈[{i0},{i1-1}] × j∈[{j0},{j1-1}]  (k={kmid})",
                  fontsize=12)

    ticklabels = [f"({ii},{jj})" for jj in range(j0,j1) for ii in range(i0,i1)]
    short_tl   = [f"({ii},{jj})" if (ii in [i0,ci.min()-1,ci.min(),ci.max(),i1-1] or
                                       jj in [j0,cj.min()-1,cj.min(),cj.max(),j1-1])
                  else "" for jj in range(j0,j1) for ii in range(i0,i1)]

    for ax3, dense, title3 in [
        (axes3[0], dense_p, "Tlaková matice pEqn"),
        (axes3[1], dense_u, "Hybnostní matice UEqn"),
    ]:
        vmax3 = np.abs(dense).max()
        im3 = ax3.imshow(dense, cmap='RdBu', vmin=-vmax3, vmax=vmax3, aspect='auto')
        plt.colorbar(im3, ax=ax3, shrink=0.8)

        # Vyznač solid buňky (mají identity řádek)
        for ri, r in enumerate(rows_sub):
            sol_loc = sol2d.flatten()[r]
            if sol_loc:
                ax3.axhline(ri-0.5, color='lime', lw=0.5, alpha=0.5)
                ax3.axhline(ri+0.5, color='lime', lw=0.5, alpha=0.5)
                ax3.axvline(ri-0.5, color='lime', lw=0.5, alpha=0.5)
                ax3.axvline(ri+0.5, color='lime', lw=0.5, alpha=0.5)
                ax3.add_patch(Rectangle((ri-0.5,ri-0.5),1,1,
                                         fc='none',ec='lime',lw=2,zorder=5))

        n_blk = len(rows_sub)
        ax3.set_xticks(range(n_blk)); ax3.set_xticklabels(short_tl, rotation=90, fontsize=6)
        ax3.set_yticks(range(n_blk)); ax3.set_yticklabels(short_tl, fontsize=6)
        ax3.set_title(title3, fontsize=11)

    plt.tight_layout()
    plt.savefig("matrix_dense_block.png", dpi=150, bbox_inches='tight')
    print("Uloženo: matrix_dense_block.png")

plt.show()
