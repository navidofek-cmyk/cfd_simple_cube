# Navierovy–Stokesovy rovnice

Řešič počítá **nestlačitelné laminární proudění**. Výchozím bodem jsou
Navierovy–Stokesovy rovnice — vyjádření zákona zachování hmoty a hybnosti pro
tekutinu.

## Zákon zachování hmoty (kontinuita)

Pro nestlačitelnou tekutinu (konstantní hustota $\rho$) říká, že kolik tekutiny
do libovolného objemu vteče, tolik musí i vytéct:

$$\nabla \cdot \mathbf{u} = 0$$

kde $\mathbf{u} = (u, v, w)$ je vektor rychlosti. Rozepsáno do složek:

$$\frac{\partial u}{\partial x} + \frac{\partial v}{\partial y} + \frac{\partial w}{\partial z} = 0$$

Tato rovnice je **klíčová obtíž** celého řešiče: neobsahuje tlak, ale přitom
musí být splněna. Algoritmus SIMPLE existuje právě proto, aby spojil tlak
s touto podmínkou (viz [SIMPLE](/page/theory/04-simple)).

## Zákon zachování hybnosti

Pro každou složku rychlosti (ustálený stav, $\partial/\partial t = 0$):

$$\underbrace{\nabla \cdot (\mathbf{u}\,\mathbf{u})}_{\text{konvekce}}
= \underbrace{-\frac{1}{\rho}\nabla p}_{\text{tlak}}
+ \underbrace{\nu\,\nabla^2 \mathbf{u}}_{\text{difúze}}$$

Jednotlivé členy:

- **Konvekce** $\nabla\cdot(\mathbf{u}\mathbf{u})$ — přenos hybnosti samotným
  prouděním. Nelineární (rychlost krát rychlost) → zdroj největších numerických
  potíží.
- **Tlakový gradient** $-\frac{1}{\rho}\nabla p$ — síla tlačící tekutinu z míst
  vysokého do nízkého tlaku.
- **Difúze** $\nu\nabla^2\mathbf{u}$ — vyhlazování rychlostního pole vlivem
  vazkosti $\nu$ (kinematická viskozita).

<div class="admonition note">
<p class="admonition-title">Konvence v kódu</p>
<p>V řešiči je hybnostní rovnice sestavena jako
<code>UEqn = div(φ, U) − ν·laplacian(U)</code> a postavena tak, aby se rovnala
zápornému tlakovému gradientu. Viz
<a href="/page/cpp/06-simple">Simple.cpp → buildUEqn</a>.</p>
</div>

## Reynoldsovo číslo

Bezrozměrné číslo charakterizující poměr setrvačných a vazkých sil:

$$\mathrm{Re} = \frac{U \, D}{\nu}$$

kde $U$ je charakteristická rychlost (zde $U_{in} = 1$ m/s) a $D$ délkový rozměr
(zde strana krychle $D = 0{,}2$ m). Z toho se v řešiči dopočítává viskozita:

$$\nu = \frac{U_{in} \, D}{\mathrm{Re}} = \frac{1 \cdot 0{,}2}{20} = 0{,}01 \;\text{m}^2/\text{s}$$

**Proč Re = 20?** Pro obtékání krychle/válce je proudění:

| Re | Charakter |
|----|-----------|
| < ~1 | plíživé (Stokesovo), symetrické |
| ~20 | **ustálené, laminární, symetrický úplav** ← náš případ |
| ~40–50 | začínají oscilace, odtrhávání vírů |
| > ~150 | von Kármánova vírová cesta (nestacionární) |

Při Re = 20 existuje **jedno ustálené řešení**, které SIMPLE najde iterací. To je
ideální pro řešič bez časové smyčky.

## Od spojitých rovnic k diskrétním

Tyto parciální diferenciální rovnice nelze řešit analyticky pro obecnou
geometrii. Proto je převedeme na soustavu algebraických rovnic metodou
**konečných objemů** — viz [další stránka](/page/theory/02-fvm).
