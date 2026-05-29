# main.cpp — vstupní bod

Nastaví geometrii, parametry a spustí výpočet. Krátký, ale ukazuje, **proč je
konfigurace oddělená od řešiče**.

## Celý soubor (zkráceně)

```cpp
#include "Simple.hpp"
#include "io.hpp"
#include <cstdio>

int main(){
    setvbuf(stdout, nullptr, _IOLBF, 0);   // řádkové bufferování → živý progress

    const int    NX=320, NY=80, NZ=80;
    const double LX=4.0, LY=1.0, LZ=1.0;
    const double CX0=0.9, CX1=1.1;          // krychle x
    const double CY0=0.4, CY1=0.6;          // krychle y
    const double CZ0=0.4, CZ1=0.6;          // krychle z
    const double RE=20.0, U_IN=1.0;

    ChannelSolver sim(NX,NY,NZ, LX,LY,LZ,
                      CX0,CX1, CY0,CY1, CZ0,CZ1,
                      RE, U_IN, 0.7, 0.3, 410, 1e-9);
    sim.snapIters = {250, 300, 400};
    sim.convergencePath = "convergence.dat";

    sim.solve();

    io::writeVTK   ("channel_cube.vtk", sim.m, sim.U, sim.p, sim.resField);
    io::writeFields("fields.dat",       sim.m, sim.U, sim.p, sim.resField);
    io::writeMesh  ("mesh.dat",         sim.m);
    return 0;
}
```

## Proč je konfigurace v main, ne v solveru

Solver (`ChannelSolver`) nezná konkrétní čísla — bere je jako parametry
konstruktoru. `main` je jediné místo, kde se rozhoduje *co* se počítá (velikost
sítě, Re, geometrie krychle). Výhoda: změna scénáře = změna `main`, ne zásah do
fyziky. Solver je **znovupoužitelný**.

## Jazykové prvky

### `const` lokální konstanty

```cpp
const int NX=320, ...;
const double RE=20.0, ...;
```

Vše, co se za běhu nemění, je `const` — překladač to může zoptimalizovat a brání
to nechtěné změně. Dokumentuje záměr („tohle je vstupní parametr").

### Konstruktor s výchozími argumenty

```cpp
ChannelSolver sim(NX,NY,NZ, LX,LY,LZ, CX0,CX1,CY0,CY1,CZ0,CZ1,
                  RE, U_IN, 0.7, 0.3, 410, 1e-9);
//                          αU   αP  maxIter tol
```

Poslední čtyři argumenty (`αU=0.7, αP=0.3, maxIter, tol`) mají v deklaraci
**výchozí hodnoty**, takže pro běžný běh se nemusí psát. Tady jsou uvedené,
protože diagnostický běh chtěl jinou toleranci a počet iterací.

### Konfigurace přes veřejné členy

```cpp
sim.snapIters = {250, 300, 400};        // seznam iterací pro snapshot
sim.convergencePath = "convergence.dat"; // kam zapisovat rezidua
```

Volitelné chování se nastaví přiřazením do veřejných členů *po* konstrukci.
Alternativa k přetížení konstruktoru o další parametry — čistší, když je
volitelných nastavení víc. `{250,300,400}` je **initializer list** přiřazený do
`std::vector<int>`.

### `setvbuf` — řádkové bufferování

```cpp
setvbuf(stdout, nullptr, _IOLBF, 0);
```

Přepne `stdout` na **řádkové bufferování** — každý `\n` se hned vypíše. Bez toho
se výstup při přesměrování do souboru drží v bloku a progress není vidět živě.
Praktická drobnost s velkým dopadem na sledování dlouhého výpočtu.

## Pořadí: solve → write

```cpp
sim.solve();                  // 1. spočítá (může trvat minuty)
io::writeVTK(...);            // 2. teprve pak zapíše výsledky
```

Výsledky se zapisují až po návratu ze `solve()`. To znamená: když výpočet
nedoběhne do konvergence (dosáhne `maxIter`), zapíše se i tak poslední stav —
`solve()` se vrátí v obou případech.

<div class="admonition tip">
<p class="admonition-title">Velikost sítě = jediné číslo</p>
<p>Změna <code>NX, NY, NZ</code> přepíná mezi hrubou (rychlou) a jemnou
(přesnou) sítí. <code>80×20×20</code> doběhne za sekundy, <code>320×80×80</code>
= 2 mil. buněk za desítky minut. Nic jiného se měnit nemusí.</p>
</div>
