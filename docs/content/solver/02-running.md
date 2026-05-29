# Spuštění a parametry

Jak řešič přeložit, spustit a co lze nastavit.

## Překlad

```bash
cd cfd_simple_cube
make
```

`Makefile` přeloží s `-O2 -std=c++20 -Wall -fopenmp`. Automaticky najde všechny
`src/*.cpp`, přeloží na objektové soubory a slinkuje do `channelCube`.

```makefile
CXXFLAGS ?= -O2 -std=c++20 -Wall -Iinclude -fopenmp
SRC := $(wildcard src/*.cpp)
OBJ := $(SRC:.cpp=.o)
channelCube: $(OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@
```

| Flag | Význam |
|------|--------|
| `-O2` | optimalizace (klíčové pro výkon) |
| `-std=c++20` | jazykový standard |
| `-fopenmp` | zapne OpenMP paralelizaci |
| `-Iinclude` | kde hledat hlavičky |

## Spuštění

```bash
OMP_NUM_THREADS=6 ./channelCube
```

`OMP_NUM_THREADS=6` omezí na 6 fyzických jader (ne 12 logických — viz
[Paralelizace](/page/theory/09-parallel)).

## Parametry v main.cpp

| Parametr | Výchozí | Význam |
|----------|---------|--------|
| `NX, NY, NZ` | 320, 80, 80 | rozlišení sítě |
| `LX, LY, LZ` | 4, 1, 1 m | rozměry domény |
| `CX0..CZ1` | krychle | hranice překážky |
| `RE` | 20 | Reynoldsovo číslo |
| `U_IN` | 1 m/s | vstupní rychlost |
| `αU, αP` | 0.7, 0.3 | podrelaxace rychlost/tlak |
| `maxIter` | dle běhu | strop iterací |
| `tol` | dle běhu | tolerance konvergence |

### Volba velikosti sítě

| Síť | Buněk | Buněk/strana krychle | Čas (6 jader) |
|-----|-------|----------------------|---------------|
| 80×20×20 | 32k | 4 | ~3 s |
| 160×40×40 | 256k | 8 | ~minuty |
| 320×80×80 | 2 mil. | 16 | desítky minut |

Jemnější síť = lepší rozlišení úplavu a hraničních vrstev, ale výrazně delší
výpočet (počet buněk i počet SIMPLE iterací roste).

### Volitelné nastavení

```cpp
sim.snapIters = {250, 300, 400};           // ulož pole v těchto iteracích
sim.convergencePath = "convergence.dat";   // živý zápis reziduí
sim.matrixExportPath = "pmatrix.dat";      // export tlakové matice (jen malé sítě)
```

## Výstupní soubory

| Soubor | Obsah | Pro |
|--------|-------|-----|
| `channel_cube.vtk` | 3D pole U, p, solid | ParaView |
| `fields.dat` | textová data všech buněk | matplotlib |
| `mesh.dat` | geometrie sítě | `plot_mesh.py` |
| `convergence.dat` | historie reziduí | `plot_convergence.py` |
| `fields_NNNN.dat` | snapshoty polí | porovnání iterací |

## Reynoldsovo číslo — co zkusit

```cpp
const double RE = 20.0;   // → změň na 40, 100, ...
```

- **Re = 20** (výchozí) — ustálené, symetrické, ~110+ iterací
- **Re = 40–50** — protáhlý úplav, na hraně stability
- **Re > 150** — vyžadovalo by časovou smyčku (von Kármánovy víry)

<div class="admonition warning">
<p class="admonition-title">Stagnace na jemné síti</p>
<p>Na 2 mil. buňkách se relativní reziduum ustálí kolem 10⁻⁵ a dál neklesá
(typické chování SIMPLE bez předpodmínění). Pro praktické účely je pole při
10⁻⁵ plně konvergované — proto se používá tolerance ~2·10⁻⁵. Viz
<a href="/page/solver/03-validation">Validace</a>.</p>
</div>
