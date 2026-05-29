"""
Grafy konvergence — parsuje log solveru a vykreslí historii reziduí.
Spuštění: python plot_convergence.py [log_soubor]
          (výchozí: convergence.log)
"""
import numpy as np
import matplotlib.pyplot as plt
import re, sys, os

logfile = sys.argv[1] if len(sys.argv) > 1 else "convergence.log"
if not os.path.exists(logfile):
    print(f"Log '{logfile}' nenalezen.")
    print("Tip: zkopíruj výstup běhu, např.  cp <task_output> convergence.log")
    sys.exit(1)

# Parsuj řádky typu:
#   iter   250   cont=1.234e-05   dU/U=5.6e-04   dP/P=7.8e-04
# i starší formát:
#   iter   250   continuity res = 1.234e-05
pat_new = re.compile(
    r"iter\s+(\d+)\s+cont=([0-9.eE+-]+)\s+dU/U=([0-9.eE+-]+)\s+dP/P=([0-9.eE+-]+)")
pat_old = re.compile(r"iter\s+(\d+)\s+continuity res\s*=\s*([0-9.eE+-]+)")

iters, cont, dU, dP = [], [], [], []
snap_iters = []
with open(logfile) as f:
    for line in f:
        m = pat_new.search(line)
        if m:
            iters.append(int(m.group(1)))
            cont.append(float(m.group(2)))
            dU.append(float(m.group(3)))
            dP.append(float(m.group(4)))
            continue
        m = pat_old.search(line)
        if m:
            iters.append(int(m.group(1)))
            cont.append(float(m.group(2)))
            dU.append(np.nan); dP.append(np.nan)
            continue
        m = re.search(r"snapshot ulozen: fields_(\d+)\.dat", line)
        if m:
            snap_iters.append(int(m.group(1)))

if not iters:
    print("Žádné iterace nenalezeny v logu.")
    sys.exit(1)

iters = np.array(iters); cont = np.array(cont)
dU = np.array(dU); dP = np.array(dP)
print(f"Načteno {len(iters)} záznamů (iter {iters.min()}–{iters.max()})")

has_fields = not np.all(np.isnan(dU))

# ── Export dat do textového souboru ──────────────────────────────────────────
datfile = "convergence.dat"
with open(datfile, "w") as f:
    f.write("# Historie konvergence SIMPLE (3D kanál + krychle)\n")
    f.write("# iter  cont_rel_L1   dU/U_rel_L2   dP/P_rel_L2\n")
    for n in range(len(iters)):
        f.write(f"{iters[n]:6d}  {cont[n]:.6e}  "
                f"{dU[n]:.6e}  {dP[n]:.6e}\n")
print(f"Uloženo: {datfile}  ({len(iters)} řádků)")

# ── Graf ──────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2 if has_fields else 1,
                         figsize=(15 if has_fields else 8, 5.5), squeeze=False)
axes = axes[0]

# Panel 1: kontinuitní reziduum
ax = axes[0]
ax.semilogy(iters, cont, 'o-', ms=3, lw=1.2, color='tab:blue',
            label="kontinuitní reziduum (rel. L1)")
ax.axhline(2e-5, color='red', ls='--', lw=1, label="tol = 2e-5")
# Vyznač minimum
imin = np.argmin(cont)
ax.plot(iters[imin], cont[imin], 'v', color='darkgreen', ms=10,
        label=f"min = {cont[imin]:.2e} @ iter {iters[imin]}")
# Snapshoty
for si in snap_iters:
    ax.axvline(si, color='gray', ls=':', lw=0.8)
    ax.text(si, ax.get_ylim()[1], f"{si}", rotation=90,
            va='top', ha='right', fontsize=7, color='gray')
ax.set_xlabel("iterace"); ax.set_ylabel("relativní reziduum")
ax.set_title("Konvergence — kontinuita")
ax.grid(True, which='both', alpha=0.3); ax.legend(fontsize=9)

# Panel 2: změna polí mezi iteracemi
if has_fields:
    ax = axes[1]
    ax.semilogy(iters, dU, 's-', ms=3, lw=1.2, color='tab:orange',
                label="dU/U  (změna rychlosti)")
    ax.semilogy(iters, dP, '^-', ms=3, lw=1.2, color='tab:green',
                label="dP/P  (změna tlaku)")
    for si in snap_iters:
        ax.axvline(si, color='gray', ls=':', lw=0.8)
        ax.text(si, ax.get_ylim()[1], f"{si}", rotation=90,
                va='top', ha='right', fontsize=7, color='gray')
    ax.set_xlabel("iterace"); ax.set_ylabel("relativní L2 změna / iter")
    ax.set_title("Konvergence — změna polí mezi iteracemi")
    ax.grid(True, which='both', alpha=0.3); ax.legend(fontsize=9)

plt.suptitle("Historie konvergence SIMPLE  (3D kanál + krychle, 320×80×80)",
             fontsize=12)
plt.tight_layout()
out = "convergence.png"
plt.savefig(out, dpi=150, bbox_inches='tight')
print(f"Uloženo: {out}")
plt.show()
