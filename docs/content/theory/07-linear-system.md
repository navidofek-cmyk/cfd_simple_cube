# Řídké matice a stencil

Po diskretizaci je každá rovnice (hybnost i tlak) lineární soustava
$\mathbf{A}\boldsymbol{\psi} = \mathbf{b}$, kde $\mathbf{A}$ je obrovská, ale
**řídká**.

## Struktura matice

Pro síť s $N = n_x n_y n_z$ buňkami je $\mathbf{A}$ rozměru $N\times N$. Pro
2 mil. buněk to je matice $2\,000\,000 \times 2\,000\,000$ — uložit ji hustě by
znamenalo $4\cdot10^{12}$ čísel (32 TB). Nemožné.

Naštěstí je **řídká**: 7bodový stencil znamená **maximálně 7 nenul na řádek**
(diagonála + 6 sousedů). Skutečná paměť je tedy $\sim 7N$ čísel.

## Uložení po diagonálách

Místo obecného řídkého formátu (CSR) využíváme strukturu mřížky: matice má
**7 diagonál** na pevných pozicích. Stačí 7 polí délky $N$:

| Pole | Diagonála | Soused |
|------|-----------|--------|
| `aP` | hlavní (0) | buňka sama |
| `aE` | +1 | i+1 (východ) |
| `aW` | −1 | i−1 (západ) |
| `aN` | +$n_x$ | j+1 (sever) |
| `aS` | −$n_x$ | j−1 (jih) |
| `aF` | +$n_x n_y$ | k+1 (předek) |
| `aB` | −$n_x n_y$ | k−1 (zadek) |

Toto je přesně to, co drží struktura [`fvMatrix`](/page/cpp/04-fvmatrix). Násobení
matice vektorem se pak dělá bez explicitní matice — jen procházením sousedů.

## Vzor řídkosti (sparsity pattern)

Pro 2D řez vypadá matice takto (modré = nenuly):

```
          sloupec →
řádek ┌─╲──────╲──────────┐   ╲ = hlavní diagonála (aP)
  ↓   │  ╲──────╲          │   sousední ╲ = aE, aW (±1)
      │   ╲──────╲─────────│   vzdálené ╲ = aN, aS (±nx)
      │    ╲──────╲        │
      └───────────────────┘   červené řádky = pevné buňky (identita)
```

Tři „pásy" odpovídají diagonálám. Pevné buňky (krychle) mají jen diagonálu
(identita). Tlaková matice je **symetrická** (Laplaceův operátor), hybnostní
**nesymetrická** (upwind konvekce posílá různé příspěvky do aE vs aW).

<div class="admonition note">
<p class="admonition-title">Vizualizace</p>
<p>Skript <code>plot_matrix.py</code> vykreslí sparsity pattern, prostorové
rozložení koeficientů a hustý blok kolem krychle. Export koeficientů zajišťují
<code>io::writePMatrix</code> a <code>io::writeUMatrix</code>.</p>
</div>

## Diagonální dominance

Pro stabilitu iteračních solverů potřebujeme **diagonální dominanci**:

$$|a_P| \geq \sum_{nb} |a_{nb}|$$

- **Tlaková matice**: $|a_P| = \sum|a_{nb}|$ uvnitř domény (slabá dominance);
  na hranici s Dirichlet (outlet) se diagonála zvětší → ostrá dominance.
- **Hybnostní matice**: upwind konvekce + podrelaxace dělají diagonálu výrazně
  dominantní (poměr 1,5–3×).

Bez dominance by Jacobi/Gauss–Seidel divergovaly a CG by ztratil
pozitivní definitnost.

## Dvě různé soustavy

| | Hybnost (`UEqn`) | Tlak (`pEqn`) |
|--|------------------|----------------|
| Typ | `fvMatrix<Vec3>` | `fvMatrix<double>` |
| Symetrie | nesymetrická | symetrická |
| Řešič | Jacobi (40 sweepů) | Conjugate Gradient |
| Definitnost | diag. dominantní | neg. semidefinitní → CG na $-\mathbf{A}$ |

Tlaková rovnice je nejdražší část (řeší se přesně každou iteraci) → používá
efektivní [CG solver](/page/theory/08-cg).
