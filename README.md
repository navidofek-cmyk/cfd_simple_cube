# cfd_simple_cube

3D řešič nestlačitelného laminárního proudění kolem krychle v hranatém tunelu.
Metoda konečných objemů (FVM), algoritmus SIMPLE, C++20, paralelizace OpenMP.

![Re=20](https://img.shields.io/badge/Re-20-blue) ![C++20](https://img.shields.io/badge/C%2B%2B-20-orange) ![OpenMP](https://img.shields.io/badge/OpenMP-paralelní-green)

## Co to počítá

Ustálené obtékání krychle (0.2³ m) ve středu tunelu 4×1×1 m při Re = 20.
Výsledek: zpomalení před krychlí, urychlení po bocích, symetrický úplav.

## Rychlý start

```bash
make
OMP_NUM_THREADS=6 ./channelCube      # spustí výpočet

python verify.py                      # ověří fyzikální správnost
python plot_fields.py                 # vizualizace proudění
python plot_mesh.py                   # vizualizace sítě
python plot_convergence.py            # grafy konvergence
```

## Struktura

```
include/       hlavičky (Mesh, Field, fvMatrix, CgSolver, Simple, Vec3, io)
src/           Simple.cpp (jádro), io.cpp, main.cpp
docs/          kompletní dokumentace (FastAPI) — teorie + C++ + řešič
*.py           vizualizace a validace
Makefile       build (g++ -O2 -std=c++20 -fopenmp)
```

## Dokumentace

Kompletní dokumentace (teorie od Navier-Stokes přes SIMPLE po CG, výklad C++ kódu
včetně „proč je tam co je") je ve složce `docs/` jako FastAPI aplikace:

```bash
cd docs
uv run uvicorn app.main:app --port 8084
# → http://localhost:8084
```

## Numerika

| Komponenta | Metoda |
|------------|--------|
| Diskretizace | Metoda konečných objemů, kolokovaná mřížka |
| Konvekce | Upwind 1. řádu |
| Tlak-rychlost | SIMPLE + Rhie–Chow interpolace |
| Tlaková rovnice | Conjugate Gradient |
| Hybnostní rovnice | Jacobi (paralelní) |
| Paralelizace | OpenMP (4,2× na 6 jádrech) |

## Výkon

| Síť | Buněk | Čas (6 jader) |
|-----|-------|---------------|
| 80×20×20 | 32k | ~3 s |
| 320×80×80 | 2 mil. | desítky minut |
