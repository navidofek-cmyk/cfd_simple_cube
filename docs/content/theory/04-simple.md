# Algoritmus SIMPLE

**SIMPLE** = *Semi-Implicit Method for Pressure-Linked Equations* (Patankar &
Spalding, 1972). Řeší zásadní problém nestlačitelných rovnic: **pro tlak
neexistuje vlastní rovnice**.

## Problém spojení tlaku a rychlosti

Máme:
- 3 hybnostní rovnice (pro $u, v, w$) — obsahují tlakový gradient $\nabla p$
- 1 rovnici kontinuity $\nabla\cdot\mathbf{u} = 0$ — tlak vůbec neobsahuje

Tlak tedy nemá svou evoluční rovnici. Jeho role je nepřímá: je to ta veličina,
která *donutí* rychlostní pole splnit kontinuitu. SIMPLE tuto vazbu řeší
iteračně.

## Princip: prediktor–korektor

Každá iterace SIMPLE má dva kroky:

1. **Prediktor** — vyřeš hybnost se *starým* tlakem → rychlost $\mathbf{u}^*$,
   která obecně **nesplňuje** kontinuitu.
2. **Korektor** — najdi tlakovou korekci $p'$, která rychlost upraví tak, aby
   kontinuitu splnila.

### Odvození tlakové rovnice

Diskrétní hybnostní rovnice pro buňku P se dá zapsat:

$$a_P \mathbf{u}_P = \underbrace{\sum_{nb} a_{nb}\mathbf{u}_{nb} + \mathbf{b}}_{\mathbf{H}(\mathbf{u})} - V\nabla p$$

Vydělíme $a_P$ a zavedeme operátor $\mathbf{H}/a_P$:

$$\mathbf{u}_P = \underbrace{\frac{\mathbf{H}(\mathbf{u})}{a_P}}_{\mathbf{HbyA}} - \frac{V}{a_P}\nabla p$$

Označme $r_{AU} = 1/a_P$. Dosadíme do kontinuity $\nabla\cdot\mathbf{u} = 0$:

$$\nabla\cdot\left(\mathbf{HbyA} - r_{AU}\,V\,\nabla p\right) = 0$$

Přeskupením vznikne **tlaková Poissonova rovnice**:

$$\boxed{\nabla\cdot(r_{AU}\,V\,\nabla p) = \nabla\cdot\mathbf{HbyA}}$$

Levá strana je Laplaceův operátor tlaku (řeší se přes [CG](/page/theory/08-cg)),
pravá je divergence rychlosti bez tlakové korekce.

## Kroky jedné iterace (mapováno na kód)

| Krok | Rovnice | V kódu (`Simple.cpp`) |
|------|---------|------------------------|
| 1 | Sestav hybnostní matici | `buildUEqn` |
| 2 | Podrelaxace | `UEqn.relax(αU)` |
| 3 | Vyřeš $\mathbf{u}^*$ | `solveWithExtra` (Jacobi) |
| 4 | Spočti $a_P$, $\mathbf{HbyA}$ | `UEqn.A()`, `UEqn.H()` |
| 5 | Toky z HbyA | `interpFlux` |
| 6 | Sestav tlakovou rovnici | `buildPEqn` |
| 7 | Vyřeš tlak | `cgSolve` |
| 8 | Korekce toků (Rhie–Chow) | viz krok 5 v solve() |
| 9 | Korekce rychlosti $\mathbf{u} = \mathbf{HbyA} - r_{AU}\nabla p$ | krok 6 v solve() |
| 10 | Podrelaxace tlaku | $p = p^{old} + \alpha_P(p - p^{old})$ |

## Podrelaxace

SIMPLE je iterační a bez tlumení by oscilovala. Proto se změny **podrelaxují**:

$$p^{new} = p^{old} + \alpha_P\,(p^{vyřešené} - p^{old})$$

s $\alpha_P = 0{,}3$ pro tlak a $\alpha_U = 0{,}7$ pro rychlost (implicitně, přes
úpravu diagonály v `relax()`). Tyto hodnoty jsou klasická konzervativní volba.

<div class="admonition note">
<p class="admonition-title">Podrelaxace rychlosti implicitně</p>
<p>Pro rychlost se nepodrelaxuje výsledek, ale rovnice: diagonála se vydělí
$\alpha_U$ a do zdroje se přidá $(a_P/\alpha_U - a_P)\,\mathbf{u}^{old}$. To je
stabilnější. Viz <code>fvMatrix::relax</code>.</p>
</div>

## Konvergence

Iterace se opakují, dokud relativní reziduum kontinuity neklesne pod toleranci:

$$\frac{\sum_c |\nabla\cdot\mathbf{u}|_c}{\text{(hodnota v 1. iteraci)}} < \text{tol}$$

Pro Re = 20 na hrubé síti (32k buněk) to je ~110 iterací; na jemné (2 mil.) se
reziduum ustálí kolem $10^{-5}$ — viz [Validace](/page/solver/03-validation).
