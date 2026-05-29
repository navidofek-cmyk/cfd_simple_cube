# Rhie–Chow interpolace

Rhie–Chow (1983) je technika, bez které **kolokovaná** mřížka (rychlost i tlak
ve středech buněk) selže kvůli oscilacím tlaku.

## Problém šachovnice (checkerboard)

Na kolokované mřížce se tlakový gradient ve středu buňky počítá centrální
diferencí přes **obě** sousední buňky:

$$\left(\frac{\partial p}{\partial x}\right)_P \approx \frac{p_E - p_W}{2\Delta x}$$

Všimni si: hodnota $p_P$ ve vzorci **chybí**. Důsledek — pole tlaku ve tvaru
šachovnice:

```
p:  10  0  10  0  10  0      ← střídání vysoká/nízká
```

dává *nulový* gradient ve všech buňkách (protože $p_E = p_W$ vždy). Řešič takové
nefyzikální pole „nevidí" a nechá ho být. Tlak a rychlost se rozpojí.

## Řešení: tok počítat přes sousední buňky

Rhie–Chow myšlenka: **objemové toky přes stěny** $\phi_f$ nepočítej z
interpolované rychlosti, ale rovnou zahrň do nich tlakový rozdíl mezi
**bezprostředně sousedícími** buňkami:

$$\phi_f = \underbrace{\tfrac{1}{2}(\mathbf{HbyA}_P + \mathbf{HbyA}_E)\cdot\mathbf{n}\,A_f}_{\text{interpolace HbyA}}
- \underbrace{\overline{r_{AU}}_f \frac{A_f}{\Delta x}(p_E - p_P)}_{\text{tlaková korekce přes }P\text{-}E}$$

Klíč je druhý člen: tlakový gradient na stěně používá $p_E - p_P$ (přímí
sousedé), takže **šachovnicový tlak už dává nenulový gradient** a řešič ho
potlačí.

$\overline{r_{AU}}_f = \tfrac{1}{2}(r_{AU,P} + r_{AU,E})$ je interpolace
koeficientu $1/a_P$ na stěnu.

## V kódu

Dvě místa v `Simple.cpp`:

**1. Interpolace HbyA na stěny** (`interpFlux`) — první člen, tok bez tlaku:

```cpp
pE[c] = 0.5*(HbyA[c].x + HbyA[e].x) * m.Ax;   // tok přes východní stěnu
```

**2. Korekce toku po vyřešení tlaku** (krok 5 v `solve()`) — odečte tlakový člen:

```cpp
phiE[c] = pE[c] - rf(c,e) * m.Ax/m.dx * (p[e] - p[c]);
//        ^HbyA tok          ^r_AU na stěně   ^p_E - p_P
```

kde `rf(c,e)` je $\overline{r_{AU}}_f$. Tento `phiE` (Rhie–Chow konzistentní tok)
se pak použije v divergenci pro reziduum i v příští hybnostní rovnici.

<div class="admonition tip">
<p class="admonition-title">Proč to funguje</p>
<p>Stejný diskrétní tlakový operátor $\frac{A_f}{\Delta x}(p_E - p_P)$ se
objevuje jak v <strong>sestavení tlakové matice</strong> (<code>buildPEqn</code>),
tak v <strong>korekci toku</strong>. Tato konzistence zaručuje, že po vyřešení
tlaku jsou toky přesně divergence-free (kontinuita splněna do strojové přesnosti
— viz reziduum $\sim10^{-16}$ ve validaci).</p>
</div>

## Souvislost s gradientem tlaku ve středu buňky

Korekce rychlosti ve středu buňky (krok 6) používá *jiný*, centrální gradient
(`gradP`):

$$\mathbf{u}_P = \mathbf{HbyA}_P - r_{AU,P}\,(\nabla p)_P$$

Rozdíl mezi tímto centrálním gradientem (pro rychlost ve středu) a stěnovým
gradientem (pro toky) je přesně tou Rhie–Chow stabilizací, která tlumí
šachovnici.
