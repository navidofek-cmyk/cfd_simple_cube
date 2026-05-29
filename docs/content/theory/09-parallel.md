# Paralelizace (OpenMP)

Na 2 milionech buněk je sériový výpočet neúnosný (přes hodinu). Paralelizace přes
**OpenMP** dala **4,2× zrychlení** na 6 jádrech.

## Proč OpenMP a ne async/vlákna

Jde o **datovou paralelizaci** — stejná operace na milionech nezávislých buněk.
To je doména OpenMP (`#pragma omp parallel for`), ne `std::async`/`std::thread`
(ty se hodí na *úlohovou* paralelizaci, např. počítat víc různých Re najednou).

```cpp
#pragma omp parallel for collapse(3)
for(int k=0;k<m.nz;++k)
for(int j=0;j<m.ny;++j)
for(int i=0;i<m.nx;++i){ ... }   // buňky se rozdělí mezi vlákna
```

`collapse(3)` sloučí tři vnořené smyčky do jedné velké iterační domény, kterou
OpenMP rovnoměrně rozdělí mezi vlákna.

## Hlavní překážka: race condition

Naivní paralelizace selže na dvou místech.

### 1. Sestavení konvekce (race)

Původní stěnový přístup: každá stěna zapisovala do **dvou** buněk:

```cpp
// ŠPATNĚ paralelně: dvě vlákna můžou psát do téže buňky
if(f>=0){ aP[c]+=f; aW[e]+=-f; }   // píše do c i do e
```

**Řešení — buňkový přístup**: každá buňka spočítá všechny své příspěvky sama,
sousedy jen *čte*:

```cpp
// západní stěna: tok si přečtu ze sousedovy phiE
double f = -phiE[w];               // jen ČTENÍ
if(f>=0) aP[c]+=f; else aW[c]+=f;  // zápis JEN do c
```

Žádné zápisy do sousedů → žádný race → bezpečné. Detaily v
[Konvekce](/page/theory/03-convection) a [Simple.cpp](/page/cpp/06-simple).

### 2. Gauss–Seidel je sériový

Gauss–Seidel čte *právě aktualizované* hodnoty sousedů — vnitřně sekvenční.
Nahrazen **Jacobi** metodou, která čte hodnoty z minulého sweepu (z odděleného
pole `tmp`):

```cpp
// Jacobi: čte stará psi, píše do tmp → plně paralelní
#pragma omp parallel for collapse(3) reduction(+:res)
for(...) { tmp[c] = (1.0/aP[c]) * (...); }
std::swap(psi.d, tmp);
```

Cena: Jacobi potřebuje o pár sweepů víc než Gauss–Seidel, ale paralelně je
celkově mnohem rychlejší.

## Redukce pro skalární součty

Sčítání přes buňky (reziduum, skalární součiny v CG) potřebuje **reduction**, aby
se sčítání z různých vláken bezpečně sloučilo:

```cpp
double cont=0.0;
#pragma omp parallel for reduction(+:cont)
for(int c=0;c<m.N();++c) cont += std::fabs(resField[c]);
```

OpenMP dá každému vláknu privátní kopii `cont`, na konci je sečte.

## Naměřené zrychlení

| Síť | 1 vlákno | 6 vláken | Zrychlení |
|-----|----------|----------|-----------|
| 80×20×20 (32k) | 12,9 s | 3,1 s | **4,2×** |

## Proč ne ideálních 6×

Úzké hrdlo je **paměťová propustnost**, ne výpočet. CFD stencil je
*memory-bound*: na každou aritmetickou operaci připadá několik čtení z paměti.
Vlákna si konkurují o sdílenou cache a paměťovou sběrnici.

Proto se běh pouští na **6 fyzických jader** (`OMP_NUM_THREADS=6`), ne na 12
logických (hyperthreading) — druhé vlákno na jádře by jen soupeřilo o tutéž
cache bez přínosu.

<div class="admonition note">
<p class="admonition-title">Shrnutí dvou klíčových rozhodnutí</p>
<p>1. <strong>Buňkový přístup</strong> u sestavení matice — eliminuje race bez
nutnosti atomických operací (pomalé) nebo barvení mřížky (složité).</p>
<p>2. <strong>Jacobi místo Gauss–Seidel</strong> — obětuje konvergenční rychlost
za plnou paralelizovatelnost.</p>
</div>
