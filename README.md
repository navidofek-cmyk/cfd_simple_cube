# cfd_simple_cube

**3D řešič nestlačitelného laminárního proudění** kolem krychle v hranatém
tunelu. Metoda konečných objemů (FVM), algoritmus **SIMPLE**, napsáno v
**C++20**, paralelizováno přes **OpenMP**. Součástí je kompletní dokumentace
(FastAPI + statická verze) a sada Python skriptů pro vizualizaci a validaci.

![Re=20](https://img.shields.io/badge/Re-20-blue)
![C++20](https://img.shields.io/badge/C%2B%2B-20-orange)
![OpenMP](https://img.shields.io/badge/OpenMP-paraleln%C3%AD-green)
![FVM](https://img.shields.io/badge/metoda-kone%C4%8Dn%C3%A9%20objemy-purple)

---

## Obsah

1. [Co řešič počítá](#co-řešič-počítá)
2. [Fyzika a numerika](#fyzika-a-numerika)
3. [Rychlý start](#rychlý-start)
4. [Spouštěcí skripty](#spouštěcí-skripty)
5. [Struktura projektu](#struktura-projektu)
6. [Dokumentace](#dokumentace)
7. [Vizualizace a validace](#vizualizace-a-validace)
8. [Parametry a velikost sítě](#parametry-a-velikost-sítě)
9. [Výkon a paralelizace](#výkon-a-paralelizace)
10. [Výstupní soubory](#výstupní-soubory)

---

## Co řešič počítá

Ustálené laminární obtékání **krychle** (0,2 × 0,2 × 0,2 m) umístěné uprostřed
tunelu se čtvercovým průřezem (4 × 1 × 1 m) při Reynoldsově čísle **Re = 20**.

```
                  WALL (no-slip)
        ┌───────────────────────────────┐
INLET   │             ┌──┐               │  OUTLET
U=1 m/s →│            │██│ krychle        │  p = 0
      → │   proud     └──┘   úplav        │  →
        └───────────────────────────────┘
                  WALL (no-slip)
```

Výsledek ukazuje typické rysy obtékání: **zpomalení** před krychlí (stagnační
bod), **urychlení** po bocích (zúžení průřezu), a symetrický **úplav** s
recirkulační zónou za krychlí. Při Re = 20 je proudění ustálené a symetrické
(žádné odtrhávání vírů).

---

## Fyzika a numerika

Řešič diskretizuje nestlačitelné Navierovy–Stokesovy rovnice:

- **kontinuita:** ∇·u = 0
- **hybnost:** ∇·(u u) = −(1/ρ)∇p + ν∇²u

| Komponenta | Metoda |
|------------|--------|
| Diskretizace | Metoda konečných objemů, kolokovaná kartézská mřížka |
| Překážka | Maska pevných buněk (immersed boundary) |
| Konvekce | Upwind schéma 1. řádu |
| Vazba tlak–rychlost | SIMPLE + Rhie–Chow interpolace |
| Tlaková rovnice | Conjugate Gradient (symetrická SPD soustava) |
| Hybnostní rovnice | Jacobi (paralelizovatelný) |
| Paralelizace | OpenMP (datová, 4,2× na 6 jádrech) |

Podrobné odvození všeho viz [dokumentace](#dokumentace).

---

## Rychlý start

```bash
# 1. přelož a spusť výpočet
make
OMP_NUM_THREADS=6 ./channelCube

# 2. ověř a vizualizuj
python3 verify.py             # 6 nezávislých fyzikálních kontrol
python3 plot_fields.py        # rychlost, tlak, proudnice
python3 plot_mesh.py          # síť + geometrie
python3 plot_convergence.py   # grafy konvergence
```

Potřebuješ `g++` (C++20, OpenMP) a Python s `numpy` + `matplotlib`.

---

## Spouštěcí skripty

Pro pohodlí jsou přiloženy bash skripty (fungují odkudkoli):

```bash
./run.sh            # přeloží + spustí řešič
./run.sh plots      # navíc spustí Python vizualizace a validaci
./run.sh docs       # spustí dokumentační web (FastAPI, port 8084)
```

---

## Struktura projektu

```
cfd_simple_cube/
├── include/                hlavičkové soubory
│   ├── Vec3.hpp            3D vektor + přetížené operátory
│   ├── Field.hpp           šablonové buňkové pole Field<T>
│   ├── Mesh.hpp            geometrie sítě + maska pevných buněk
│   ├── fvMatrix.hpp        7bodový stencil (šablona), Jacobi řešič
│   ├── CgSolver.hpp        Conjugate Gradient + workspace
│   ├── Simple.hpp          deklarace ChannelSolver
│   └── io.hpp              deklarace vstup/výstup
├── src/
│   ├── Simple.cpp          JÁDRO — SIMPLE smyčka, sestavení matic
│   ├── io.cpp              zápis VTK / dat / matic
│   └── main.cpp            vstupní bod, parametry úlohy
├── docs/                   DOKUMENTACE (FastAPI)
│   ├── app/                FastAPI aplikace (main, render, nav)
│   ├── content/            obsah v Markdownu (teorie, C++, řešič)
│   ├── templates/          Jinja2 šablony
│   ├── static/             CSS, JS
│   ├── site/               STATICKÝ web (generovaný, pro GitHub Pages)
│   ├── build_static.py     generátor statického webu
│   └── run.sh              spuštění dokumentace
├── plot_fields.py          vizualizace proudění
├── plot_mesh.py            vizualizace sítě
├── plot_convergence.py     grafy konvergence
├── plot_matrix.py          struktura řídkých matic
├── verify.py               validace (bilance, symetrie, no-slip…)
├── run.sh                  hlavní spouštěcí skript
├── Makefile                build (g++ -O2 -std=c++20 -fopenmp)
└── .github/workflows/      automatické nasazení docs na GitHub Pages
```

---

## Dokumentace

Kompletní dokumentace má tři části:

- **Teorie** (9 stránek) — od Navier–Stokes přes FVM diskretizaci, SIMPLE,
  Rhie–Chow interpolaci, řídké matice až po Conjugate Gradient a paralelizaci.
- **C++ implementace** (10 stránek) — výklad každého souboru řádek po řádku
  **plus „proč je tam co je"**: proč OOP, proč šablony, proč přetěžování
  operátorů, RAII, lambda, OpenMP pragmy…
- **Řešič** (3 stránky) — SIMPLE smyčka krok po kroku, spuštění, validace.

### Live verze (FastAPI)

```bash
cd docs
uv run uvicorn app.main:app --host 0.0.0.0 --port 8084
# → http://localhost:8084   nebo   http://<IP>:8084
```

Spravováno přes [`uv`](https://docs.astral.sh/uv/). Obsah je v Markdownu,
rovnice přes MathJax, zvýraznění kódu přes Pygments.

### Statická verze (GitHub Pages)

Statický web (čisté HTML, žádný server) **včetně galerie výsledků** se generuje:

```bash
cd docs
uv run python build_static.py     # → docs/site/
```

Po pushnutí na `main` ho GitHub Action automaticky nasadí na GitHub Pages
(`.github/workflows/pages.yml`). Stačí v nastavení repozitáře zapnout
**Settings → Pages → Source: GitHub Actions**.

Lokálně otevřeš třeba:

```bash
cd docs/site && python3 -m http.server 8099   # → http://localhost:8099
```

---

## Vizualizace a validace

### `verify.py` — 6 fyzikálních kontrol

| # | Kontrola | Práh |
|---|----------|------|
| 1 | Globální bilance hmoty | < 10⁻⁶ |
| 2 | Max lokální \|div u\| | < 10⁻⁴ |
| 3 | Symetrie kolem y = ½ | < 5·10⁻³ |
| 4 | Symetrie kolem z = ½ | < 5·10⁻³ |
| 5 | Tlakový pokles (Δp < 0) | — |
| 6 | No-slip uvnitř krychle (\|U\| = 0) | < 10⁻¹⁰ |

### Vizualizační skripty

| Skript | Vstup | Výstup |
|--------|-------|--------|
| `plot_fields.py` | `fields.dat` | rychlost, tlak, proudnice (řezy XY/XZ) |
| `plot_mesh.py` | `mesh.dat` | 3D schéma tunelu + detail sítě |
| `plot_convergence.py` | `convergence.dat` | historie reziduí + změny polí |
| `plot_matrix.py` | `pmatrix.dat`, `umatrix.dat` | sparsity pattern, koeficienty |

Pro plné 3D zkoumání lze `channel_cube.vtk` otevřít v **ParaView**.

---

## Parametry a velikost sítě

Vše se nastavuje v `src/main.cpp`:

| Parametr | Význam |
|----------|--------|
| `NX, NY, NZ` | rozlišení sítě |
| `LX, LY, LZ` | rozměry domény (m) |
| `CX0..CZ1` | hranice krychle |
| `RE` | Reynoldsovo číslo |
| `U_IN` | vstupní rychlost |
| `αU, αP` | podrelaxace (0,7 / 0,3) |

| Síť | Buněk | Buněk/strana krychle | Čas (6 jader) |
|-----|-------|----------------------|---------------|
| 80×20×20 | 32 tis. | 4 | ~3 s |
| 160×40×40 | 256 tis. | 8 | minuty |
| 320×80×80 | 2 mil. | 16 | desítky minut |

**Reynoldsovo číslo:** Re = 20 dává ustálené symetrické proudění. Při Re ≈ 40–50
se úplav protahuje; nad ~150 by bylo potřeba časové schéma (von Kármánovy víry).

---

## Výkon a paralelizace

Paralelizace přes OpenMP dala **4,2× zrychlení** na 6 jádrech (32k buněk:
12,9 s → 3,1 s). Dvě klíčová rozhodnutí:

1. **Buňkový přístup** při sestavení matic — každá buňka zapisuje jen do sebe,
   sousedy jen čte → žádný *race condition* bez nutnosti atomik.
2. **Jacobi místo Gauss–Seidel** pro hybnostní rovnici — plně paralelizovatelné
   (čte hodnoty z minulého sweepu).

Úzkým hrdlem je paměťová propustnost (CFD je *memory-bound*), proto se běh pouští
na fyzická jádra (`OMP_NUM_THREADS=6`), ne na logická (hyperthreading).

---

## Výstupní soubory

| Soubor | Obsah |
|--------|-------|
| `channel_cube.vtk` | 3D pole U, p, solid, reziduum (ParaView) |
| `fields.dat` | textová data všech buněk (Python) |
| `mesh.dat` | geometrie sítě |
| `convergence.dat` | historie reziduí (iter, kontinuita, ΔU, Δp) |
| `fields_NNNN.dat` | snapshoty polí v zadaných iteracích |
| `pmatrix.dat`, `umatrix.dat` | koeficienty matic (jen malé sítě) |

> Výstupní data a build artefakty jsou v `.gitignore` — regenerují se spuštěním.

---

## Licence / autorství

Vzdělávací projekt (učení C++ a CFD). Implementace metody konečných objemů a
algoritmu SIMPLE od základu, inspirovaná konvencemi OpenFOAM.
