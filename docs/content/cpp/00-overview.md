# Architektura projektu

Řešič je rozdělen do hlavičkových souborů (`include/`) a zdrojových (`src/`).
Návrh kopíruje styl OpenFOAM — operátory `fvm::`/`fvc::` sestavují matice,
`fvMatrix` je drží, solver je iteruje.

## Struktura souborů

```
cfd_simple_cube/
├── include/
│   ├── Vec3.hpp       hodnotový typ 3D vektoru + operátory
│   ├── Field.hpp      šablonové pole (cell-centered)
│   ├── Mesh.hpp       geometrie + maska pevných buněk
│   ├── fvMatrix.hpp   7bodový stencil, šablona <T>
│   ├── CgSolver.hpp   Conjugate Gradient + workspace
│   ├── Simple.hpp     deklarace ChannelSolver
│   └── io.hpp         deklarace vstup/výstup
├── src/
│   ├── Simple.cpp     JÁDRO — SIMPLE smyčka, sestavení matic
│   ├── io.cpp         zápis VTK / dat / matic
│   └── main.cpp       vstupní bod, parametry
└── Makefile
```

## Datový tok

```
        ┌──────────┐
        │  Vec3    │  základní typ (x,y,z)
        └────┬─────┘
             │ používá
        ┌────▼─────┐
        │ Field<T> │  Field<Vec3>=rychlost, Field<double>=tlak
        └────┬─────┘
             │
   ┌─────────┼──────────┐
   │         │          │
┌──▼───┐ ┌───▼────┐ ┌───▼─────┐
│ Mesh │ │fvMatrix│ │CgSolver │
└──┬───┘ └───┬────┘ └───┬─────┘
   └─────────┼──────────┘
        ┌────▼─────────┐
        │ ChannelSolver│  drží stav + solve()
        └────┬─────────┘
             │
        ┌────▼──┐
        │  io   │  → .vtk, .dat
        └───────┘
```

## Návrhové principy

<div class="admonition note">
<p class="admonition-title">Klíčové vzory</p>
<p><strong>Šablony nad typem pole</strong> — <code>fvMatrix&lt;T&gt;</code> funguje
pro skalár (tlak) i vektor (rychlost) díky tomu, že koeficienty jsou vždy
skalární a jen pravá strana je typu <code>T</code>.</p>
<p><strong>Maska místo nepravidelné sítě</strong> — krychle není vyříznutá
z mřížky, ale označená bitovou maskou <code>Mesh::solid</code>. Síť zůstává
strukturovaná (rychlá indexace).</p>
<p><strong>Předalokovaný workspace</strong> — CG nealokuje paměť při každém
volání, drží ji v <code>CgWorkspace</code>.</p>
</div>

## Závislosti mezi hlavičkami

```
Vec3.hpp   ← Field.hpp ← fvMatrix.hpp ← CgSolver.hpp
                ↑              ↑
            Mesh.hpp ──────────┘
                ↑
           Simple.hpp, io.hpp
```

Žádné cyklické závislosti. `Vec3.hpp` je listový (nic nezávisí pod ním),
`Simple.hpp` je vrchol.

## Build

```bash
make            # přeloží s -O2 -std=c++20 -fopenmp
OMP_NUM_THREADS=6 ./channelCube
```

`Makefile` automaticky najde všechny `src/*.cpp`, přeloží na `.o` a slinkuje.
Detaily v [main.cpp](/page/cpp/08-main) a [Spuštění](/page/solver/02-running).

V dalších kapitolách projdeme každý soubor zvlášť, počínaje nejjednodušším —
[Vec3.hpp](/page/cpp/01-vec3).
