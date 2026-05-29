"""
Ověření konzistence CFD výsledků.
Spuštění: python verify.py
"""
import numpy as np
import sys, os

DAT = "fields.dat"
if not os.path.exists(DAT):
    print(f"'{DAT}' nenalezen"); sys.exit(1)

data  = np.loadtxt(DAT, comments="#")
i_arr = data[:,0].astype(int)
j_arr = data[:,1].astype(int)
k_arr = data[:,2].astype(int)
x=data[:,3]; y=data[:,4]; z=data[:,5]
ux=data[:,6]; uy=data[:,7]; uz=data[:,8]
p=data[:,9]; res=data[:,10]; solid=data[:,11].astype(int)

nx=i_arr.max()+1; ny=j_arr.max()+1; nz=k_arr.max()+1
dx=x[1]-x[0]; dy=y[nx]-y[0]; dz=z[nx*ny]-z[0]
Lx=nx*dx; Ly=ny*dy; Lz=nz*dz
Ax=dy*dz; V=dx*dy*dz

print(f"Síť: {nx}×{ny}×{nz}   Lx={Lx:.2f} Ly={Ly:.2f} Lz={Lz:.2f}")
print(f"Pevných buněk: {solid.sum()}\n")

PASS = "  ✓"; FAIL = "  ✗"
ok_all = True

# ─── 1. BILANCE HMOTY ────────────────────────────────────────────────────────
print("═" * 60)
print("1. BILANCE HMOTY")

# Správná kontrola: globální suma div(U) přes všechny buňky.
# Lokální div(U)_c = phiE[c] - phiW[c] + ... Pro uzavřený kanál
# se všechny vnitřní stěny vyruší a zůstane: Q_out - Q_in.
# Kontinuitní rezidua v res[] jsou právě tato lokální div(U).
sum_div_all   = res.sum()
sum_div_fluid = res[solid==0].sum()

# Průměrná rychlost na výstupu (cell-centered proxy)
mask_out = (i_arr == nx-1) & (solid == 0)
Q_out_cc = ux[mask_out].sum() * Ax
Q_in_bc  = 1.0 * (ny*nz) * Ax    # U_in=1 × celý průřez (bez solid na inlet)

ok = abs(sum_div_fluid) < 1e-6
ok_all &= ok
print(f"  Globální Σ div(U) přes fluid buňky: {sum_div_fluid:.3e}")
print(f"  (= Q_out − Q_in na stěnách, má být ≈ 0)")
print(f"  |Σ div| < 1e-6{PASS if ok else FAIL}")
print(f"  Q_outlet (cell-centered proxy): {Q_out_cc:.6f}  (BC inlet: {Q_in_bc:.6f})")
print(f"  Pozn.: cell-centered ≠ face flux – rozdíl je normální")

# ─── 2. KONTINUITNÍ REZIDUUM ─────────────────────────────────────────────────
print("\n" + "═" * 60)
print("2. KONTINUITNÍ REZIDUUM (div U)")

res_fluid = np.abs(res[solid==0])
print(f"  Max |div U| ve fluidních buňkách : {res_fluid.max():.3e}")
print(f"  Průměr |div U|                   : {res_fluid.mean():.3e}")
print(f"  L2 norma                         : {np.sqrt((res_fluid**2).mean()):.3e}")
ok2 = res_fluid.max() < 1e-4
ok_all &= ok2
print(f"  Max < 1e-4 {PASS if ok2 else FAIL}")

# ─── 3. SYMETRIE KOLEM y=Ly/2 ────────────────────────────────────────────────
print("\n" + "═" * 60)
print("3. SYMETRIE kolem y = Ly/2  (platí pro centrovanou krychli a Re=20)")

# Pro každou buňku (i,j,k) najdi zrcadlovou (i, ny-1-j, k)
# Porovnej u_x(i,j,k) s u_x(i, ny-1-j, k) – musí být stejné
# Porovnej u_y(i,j,k) s -u_y(i, ny-1-j, k) – musí být opačné
idx_3d = np.zeros((nx,ny,nz), dtype=int)
for n in range(len(i_arr)):
    idx_3d[i_arr[n], j_arr[n], k_arr[n]] = n

errs_ux = []; errs_uy = []
for ki in range(nz):
    for ji in range(ny//2):
        ji_mirror = ny-1-ji
        for ii in range(nx):
            c  = idx_3d[ii, ji, ki]
            cm = idx_3d[ii, ji_mirror, ki]
            if solid[c] or solid[cm]: continue
            ref = max(abs(ux[c]), abs(ux[cm]), 1e-10)
            errs_ux.append(abs(ux[c] - ux[cm]) / ref)
            ref2 = max(abs(uy[c]), abs(uy[cm]), 1e-10)
            errs_uy.append(abs(uy[c] + uy[cm]) / ref2)   # uy musí být opačné

e_ux = np.max(errs_ux) if errs_ux else 0
e_uy = np.max(errs_uy) if errs_uy else 0
tol_sym = 5e-3   # pro hrubou síť (20 buněk na průřez) je 5e-3 realistické
ok3a = e_ux < tol_sym; ok3b = e_uy < tol_sym; ok_all &= ok3a & ok3b
print(f"  Max rel. chyba u_x(y) vs u_x(Ly-y) : {e_ux:.3e}{PASS if ok3a else FAIL}")
print(f"  Max rel. chyba u_y(y) vs -u_y(Ly-y): {e_uy:.3e}{PASS if ok3b else FAIL}")

# ─── 4. SYMETRIE KOLEM z=Lz/2 ────────────────────────────────────────────────
print("\n" + "═" * 60)
print("4. SYMETRIE kolem z = Lz/2")

errs_ux2 = []; errs_uz = []
for ki in range(nz//2):
    ki_mirror = nz-1-ki
    for ji in range(ny):
        for ii in range(nx):
            c  = idx_3d[ii, ji, ki]
            cm = idx_3d[ii, ji, ki_mirror]
            if solid[c] or solid[cm]: continue
            ref = max(abs(ux[c]), abs(ux[cm]), 1e-10)
            errs_ux2.append(abs(ux[c] - ux[cm]) / ref)
            ref2 = max(abs(uz[c]), abs(uz[cm]), 1e-10)
            errs_uz.append(abs(uz[c] + uz[cm]) / ref2)

e_ux2 = np.max(errs_ux2) if errs_ux2 else 0
e_uz  = np.max(errs_uz)  if errs_uz  else 0
ok4a = e_ux2 < tol_sym; ok4b = e_uz < tol_sym; ok_all &= ok4a & ok4b
print(f"  Max rel. chyba u_x(z) vs u_x(Lz-z) : {e_ux2:.3e}{PASS if ok4a else FAIL}")
print(f"  Max rel. chyba u_z(z) vs -u_z(Lz-z): {e_uz:.3e}{PASS if ok4b else FAIL}")

# ─── 5. TLAKOVÁ BILANCE (impulzový teorém) ────────────────────────────────────
print("\n" + "═" * 60)
print("5. IMPULZOVÁ BILANCE v ose x")

# Celková síla na kanál v ose x:
#   F_inlet (tlak) + F_outlet (tlak) + F_walls (viskozita, ≈0 pro plochá dna)
#   = změna hybnostního toku (= 0 pro ustálené proudění s Neumann na outlet)
#
# Jednoduchý test: tlak klesá od vstupu k výstupu (dp < 0)
# Inlet: průměrný tlak
mask_in   = (i_arr == 0) & (solid == 0)
p_in_avg  = p[mask_in].mean()
p_out_avg = p[mask_out].mean()   # ≈ 0 (Dirichlet BC)
dp = p_out_avg - p_in_avg
print(f"  Průměrný tlak inlet  : {p_in_avg:.4f} Pa")
print(f"  Průměrný tlak outlet : {p_out_avg:.4f} Pa  (BC = 0)")
print(f"  Tlakový pokles Δp    : {dp:.4f} Pa")
ok5 = dp < 0
ok_all &= ok5
print(f"  Δp < 0 (tlak klesá po proudu){PASS if ok5 else FAIL}")

# ─── 6. MAXIMÁLNÍ RYCHLOST ───────────────────────────────────────────────────
print("\n" + "═" * 60)
print("6. RYCHLOSTNÍ POLE")
spd = np.sqrt(ux**2 + uy**2 + uz**2)
print(f"  Max |U| celkové         : {spd.max():.4f} m/s")
print(f"  Max |U| fluid           : {spd[solid==0].max():.4f} m/s")
print(f"  Max |U| solid (=0?)     : {spd[solid==1].max():.6f} m/s")
ok6 = spd[solid==1].max() < 1e-10
ok_all &= ok6
print(f"  U_solid ≈ 0{PASS if ok6 else FAIL}")
print(f"  Min u_x fluid           : {ux[solid==0].min():.4f} m/s")
print(f"  Max u_x fluid           : {ux[solid==0].max():.4f} m/s")

# ─── SOUHRN ──────────────────────────────────────────────────────────────────
print("\n" + "═" * 60)
if ok_all:
    print("  VÝSLEDEK: všechny kontroly PROŠLY ✓")
else:
    print("  VÝSLEDEK: některé kontroly SELHALY ✗")
print("═" * 60)
