# Metoda konečných objemů (FVM)

FVM převádí spojité PDE na algebraické rovnice tím, že doménu rozdělí na malé
**kontrolní objemy** (buňky) a v každém integruje zákony zachování.

## Diskretizace domény

Doména 4×1×1 m se rozdělí na uniformní kartézskou mřížku $n_x \times n_y \times n_z$
buněk. Každá buňka je kvádr o rozměrech $\Delta x \times \Delta y \times \Delta z$.
V řešiči je síť **kolokovaná** — všechny veličiny (rychlost i tlak) jsou uloženy
ve **středech buněk**.

```
        N (j+1)
        │
 W ───┌─┼─┐─── E          každá buňka má 6 sousedů ve 3D:
(i-1) │ P │ (i+1)         W/E (osa x), S/N (osa y), B/F (osa z)
      └─┼─┘
        │
        S (j-1)        (+ B = k-1 vzadu, F = k+1 vepředu)
```

Indexace v kódu (`Mesh::id`):

$$\text{id}(i,j,k) = i + j\,n_x + k\,n_x n_y$$

což je lineární index do 1D pole — řádkové uspořádání (x se mění nejrychleji).

## Integrace přes kontrolní objem

Klíčový krok: rovnici zintegrujeme přes objem buňky $V$ a použijeme
**Gaussovu–Ostrogradského větu** (objemový integrál divergence = plošný integrál
přes stěny):

$$\int_V \nabla \cdot \mathbf{F} \, dV = \oint_{\partial V} \mathbf{F} \cdot \mathbf{n} \, dS
\approx \sum_{f} \mathbf{F}_f \cdot \mathbf{n}_f \, A_f$$

Tedy **objemový integrál se změní na součet toků přes 6 stěn**. To je podstata
FVM: pracujeme s toky přes stěny, ne s derivacemi v bodech.

### Difúzní člen

Difúzi $\nu\nabla^2\phi$ zintegrujeme:

$$\int_V \nu\nabla^2\phi\,dV = \sum_f \nu (\nabla\phi)_f \cdot \mathbf{n}_f A_f$$

Gradient na stěně mezi buňkami $P$ a $E$ aproximujeme centrální diferencí:

$$(\nabla\phi)_f \cdot \mathbf{n}_f \approx \frac{\phi_E - \phi_P}{\Delta x}$$

Příspěvek stěny do rovnice buňky $P$ je tedy:

$$\nu \frac{A_x}{\Delta x}(\phi_E - \phi_P)$$

V kódu je $A_x = \Delta y\,\Delta z$ (plocha stěny kolmé na x) a koeficient
$g_x = \nu A_x/\Delta x$. Tento člen přispívá:

- $+g_x$ k diagonále buňky P (člen $-\phi_P$ se záporným znaménkem → kladná diagonála)
- $-g_x$ k mimodiagonále (soused E)

### Plochy stěn ve 3D

| Stěna | Normála | Plocha |
|-------|---------|--------|
| W, E | osa x | $A_x = \Delta y \, \Delta z$ |
| S, N | osa y | $A_y = \Delta x \, \Delta z$ |
| B, F | osa z | $A_z = \Delta x \, \Delta y$ |

Objem buňky: $V = \Delta x\,\Delta y\,\Delta z$. Tyto hodnoty drží struktura
[`Mesh`](/page/cpp/03-mesh).

## Výsledný tvar: 7bodový stencil

Po diskretizaci difúze i konvekce dostane každá buňka rovnici tvaru:

$$a_P \phi_P = a_E\phi_E + a_W\phi_W + a_N\phi_N + a_S\phi_S + a_F\phi_F + a_B\phi_B + b$$

To je **7bodový stencil** (buňka P + 6 sousedů). Koeficienty $a_E,\dots,a_B$ a
zdroj $b$ jsou přesně to, co ukládá šablonová struktura
[`fvMatrix`](/page/cpp/04-fvmatrix). Tato lineární soustava se pak řeší
iteračně.

## Pevná překážka (krychle)

Buňky uvnitř krychle jsou označeny jako **solid** (maska `Mesh::solid`).
V matici dostanou **identitu** ($a_P = 1$, ostatní = 0, $b = 0$), takže jejich
hodnota je nulová (no-slip). Sousední tekutinové buňky vidí stěnu krychle jako
Dirichletovu okrajovou podmínku — viz [Okrajové podmínky](/page/theory/06-boundary).
