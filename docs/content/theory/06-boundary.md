# Okrajové podmínky

Doména má pět typů hranic. Každá vyžaduje jiné zacházení pro rychlost a tlak.

## Přehled

```
                  WALL (no-slip)
        ┌───────────────────────────────┐
INLET   │             ┌──┐              │  OUTLET
U=(U_in,│             │██│ CUBE         │  ∂U/∂n=0
0,0)  → │             └──┘ (no-slip)    │  p=0
∂p/∂n=0 │                               │
        └───────────────────────────────┘
                  WALL (no-slip)
```

| Hranice | Rychlost | Tlak |
|---------|----------|------|
| Inlet (x=0) | Dirichlet $\mathbf{u}=(U_{in},0,0)$ | Neumann $\partial p/\partial n = 0$ |
| Outlet (x=L) | Neumann $\partial \mathbf{u}/\partial n = 0$ | Dirichlet $p = 0$ |
| Stěny (±y, ±z) | Dirichlet $\mathbf{u}=0$ (no-slip) | Neumann |
| Krychle | Dirichlet $\mathbf{u}=0$ (no-slip) | Neumann |

## Dirichletova podmínka (předepsaná hodnota)

Pro stěnu a krychli (no-slip, $\mathbf{u}=0$) a inlet ($\mathbf{u}=U_{in}$) je
hodnota na hranici známá. Difúzní člen na hraniční stěně používá vzdálenost
**poloviny buňky** (střed buňky → stěna):

$$\nu\frac{A}{\Delta x/2}(\phi_{stěna} - \phi_P)$$

Koeficient je tedy dvojnásobný ($g_b = \nu A/(0{,}5\,\Delta x)$). V kódu:

```cpp
eq.aP[c] += gbx;                       // přidá do diagonály
eq.src[c] += gbx * Uin;                // známá hodnota → do zdroje
```

Pro no-slip ($\mathbf{u}=0$) zdrojový člen vypadne (násobí se nulou), zůstane jen
příspěvek do diagonály.

### Inlet — i konvektivní tok

Na vstupu je znám i tok ($U_{in} \cdot A_x$), takže do zdroje přibude i
konvektivní příspěvek:

```cpp
eq.src[c] += gbx * Uin;          // difúze ze vstupu
eq.src[c] += U_in * m.Ax * Uin;  // konvekce ze vstupu
```

## Neumannova podmínka (předepsaný gradient)

**Outlet pro rychlost** ($\partial\mathbf{u}/\partial n = 0$): rychlost na výstupu
se nemění → stěna nepřispívá k difúzi, jen ke konvekci (odtok). Tekutina volně
odtéká.

**Stěny/inlet pro tlak** ($\partial p/\partial n = 0$): tlak nemá normálový
gradient → příslušná stěna do tlakové matice **nepřispívá vůbec**. V kódu se ten
směr prostě přeskočí.

## Outlet pro tlak — Dirichlet p = 0

Tlak na výstupu je ukotven na nulu (referenční hodnota). To zároveň **odstraňuje
singularitu** tlakové rovnice: čistě Neumannův problém (tlak všude jen
gradientem) má řešení určené až na konstantu. Ukotvení $p=0$ na výstupu udělá
soustavu jednoznačnou.

```cpp
// outlet: Dirichlet p_stěna = 0, vzdálenost dx/2
double g = rAU[c] * m.Ax/(0.5*m.dx);
L.aP[c] -= g;   // src += g*0 = 0
```

## Krychle jako vnitřní hranice

Buňky uvnitř krychle nejsou součástí výpočtu. V matici dostanou **identitu**:

```cpp
if (m.solid[c]) { eq.aP[c] = 1.0; continue; }   // → U = 0, p = 0
```

Tekutinová buňka sousedící s pevnou buňkou s ní zachází jako s no-slip
Dirichletovou stěnou (difúze přes půl buňky, žádný tok). Tím vznikne přirozeně
povrch krychle s nulovou rychlostí.

<div class="admonition note">
<p class="admonition-title">Rohy a hrany krychle</p>
<p>Tekutinová buňka v rohu sousedí s pevnými buňkami ze 2–3 stran a dostane
příspěvek $g_b$ z každé z nich. Je to fyzikálně správné — buňka „cítí" stěnu
z více směrů a její rychlost je silněji tlačena k nule.</p>
</div>
