# 3D CFD řešič — kanál s krychlí

Dokumentace k vlastnímu řešiči nestlačitelného laminárního proudění metodou
konečných objemů (FVM) s algoritmem **SIMPLE**, napsanému v C++20 a
paralelizovanému přes OpenMP.

## Co řešič počítá

Ustálené laminární obtékání **krychle** umístěné uprostřed hranatého tunelu
(kanálu se čtvercovým průřezem):

```
        ┌─────────────────────────────────────┐
 U_in → │            ┌──┐                       │ → outlet (p=0)
   →    │   proud    │██│   úplav               │   →
   →    │            └──┘                       │
        └─────────────────────────────────────┘
              stěny (no-slip) nahoře/dole/bok
```

| Veličina | Hodnota |
|----------|---------|
| Doména | 4 × 1 × 1 m |
| Krychle | 0.2 × 0.2 × 0.2 m, střed (1.0, 0.5, 0.5) |
| Reynoldsovo číslo | Re = 20 (na základě strany krychle) |
| Vstupní rychlost | U_in = 1 m/s |
| Síť | až 320 × 80 × 80 = 2 mil. buněk |

## Jak je dokumentace členěná

<div class="admonition note">
<p class="admonition-title">Tři části</p>
<p><strong>Teorie</strong> — matematika od Navierových–Stokesových rovnic přes
FVM diskretizaci, algoritmus SIMPLE, Rhie–Chow interpolaci až po řešení
řídkých soustav metodou sdružených gradientů.</p>
<p><strong>C++ implementace</strong> — výklad každého souboru řádek po řádku,
včetně detailního vysvětlení použitých jazykových prvků (šablony, struktury,
přetěžování operátorů, RAII, lambda výrazy, OpenMP pragmy).</p>
<p><strong>Řešič</strong> — jak do sebe kusy zapadají, spuštění, parametry,
validace.</p>
</div>

## Architektura jedním pohledem

```
main.cpp
   │  vytvoří
   ▼
ChannelSolver ──── Mesh           (geometrie + maska pevných buněk)
   │           ├── Field<Vec3> U  (rychlost)
   │           ├── Field<double> p (tlak)
   │           └── CgWorkspace
   │
   │  solve() = SIMPLE smyčka:
   │     1. buildUEqn  → fvMatrix<Vec3>   (hybnost)
   │     2. solveWithExtra (Jacobi)        → U*
   │     3. buildPEqn  → fvMatrix<double> (tlak)
   │     4. cgSolve    (Conjugate Gradient) → p
   │     5. Rhie–Chow korekce toků
   │     6. korekce rychlosti
   │
   ▼
io::writeVTK / writeFields → vizualizace (ParaView, matplotlib)
```

## Technologie

- **C++20** — šablony, `std::vector`, RAII, agregátní inicializace
- **OpenMP** — datová paralelizace (4.2× zrychlení na 6 jádrech)
- **Python** — vizualizace (matplotlib) a tato dokumentace (FastAPI)
- **uv** — správa Python prostředí

Začni částí [Teorie → Navierovy–Stokesovy rovnice](/page/theory/01-navier-stokes).
