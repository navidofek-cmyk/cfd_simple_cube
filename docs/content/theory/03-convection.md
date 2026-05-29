# Konvekce a upwind schéma

Konvektivní člen $\nabla\cdot(\mathbf{u}\,\phi)$ je nelineární a numericky
citlivý. Způsob jeho diskretizace rozhoduje o stabilitě celého řešiče.

## Tok přes stěnu

Po integraci přes buňku (Gaussova věta) dostaneme součet konvektivních toků:

$$\int_V \nabla\cdot(\mathbf{u}\phi)\,dV = \sum_f \phi_f \, (\mathbf{u}_f\cdot\mathbf{n}_f)\,A_f
= \sum_f \phi_f \, F_f$$

kde $F_f = (\mathbf{u}_f\cdot\mathbf{n}_f)\,A_f$ je **objemový tok** přes stěnu
(kolik tekutiny protéká stěnou za jednotku času). V kódu jsou tyto toky uloženy
jako `phiE`, `phiN`, `phiF` (tok přes východní / severní / přední stěnu vlastníka).

Otázka je: jakou hodnotu $\phi_f$ dosadit na stěnu, když máme hodnoty jen ve
středech buněk?

## Centrální schéma — a proč nestačí

Nejpřesnější (2. řád) by bylo lineární interpolovat:

$$\phi_f = \tfrac{1}{2}(\phi_P + \phi_E)$$

Pro konvekci ale centrální schéma vede k **nefyzikálním oscilacím** při vyšším
toku (záporné koeficienty v matici, ztráta diagonální dominance). Řešení by
mohlo kmitat a divergovat.

## Upwind schéma (1. řád)

Robustní volba: hodnotu na stěně vezmeme z buňky **proti proudu** (odkud
tekutina teče):

$$\phi_f = \begin{cases}
\phi_P & \text{když } F_f \geq 0 \quad (\text{teče z } P \text{ ven}) \\
\phi_E & \text{když } F_f < 0 \quad (\text{teče dovnitř } P)
\end{cases}$$

To zaručuje, že příspěvky do matice mají správné znaménko a diagonála zůstává
dominantní → stabilní řešení.

### Příspěvky do matice

Pro východní stěnu buňky P (s tokem $F = \phi_E^{\text{flux}}$):

```
if (F >= 0) { aP += F;  aW[soused] += -F; }   // odtok z P
else        { aE += F;  aP[soused] += -F; }   // přítok do P
```

Tento vzorec se v kódu objevuje v
[`buildUEqn`](/page/cpp/06-simple) pro každou ze 6 stěn.

<div class="admonition warning">
<p class="admonition-title">Cena za stabilitu</p>
<p>Upwind 1. řádu zavádí <strong>numerickou difúzi</strong> — uměle vyhlazuje
ostré gradienty, jako by byla viskozita o něco vyšší. Pro Re = 20 na jemné síti
je to zanedbatelné; pro vysoké Re by bylo potřeba schéma vyššího řádu (např.
linear-upwind, QUICK).</p>
</div>

## Buňkový vs. stěnový přístup (důležité pro paralelizaci)

Konvekci lze sestavit dvěma způsoby:

**Stěnový** (původní): projdu stěny, každá zapíše do **dvou** buněk (vlastník +
soused). Jednoduché, ale při paralelizaci vzniká *race condition* — dvě vlákna
píšou do stejné buňky.

**Buňkový** (paralelní): každá buňka projde svých 6 stěn a zapíše **jen do
sebe**. Tok přes západní stěnu si přečte ze sousedovy hodnoty `phiE[w]`. Žádné
zápisy do sousedů → bezpečné pro OpenMP.

Tento přepis byl klíčový pro paralelizaci — viz
[Paralelizace](/page/theory/09-parallel) a
[Simple.cpp](/page/cpp/06-simple).
