# Validace a vizualizace

Jak ověřit, že výsledky jsou fyzikálně správné, a jak je zobrazit.

## Ověřovací skript

`verify.py` provádí šest nezávislých kontrol na `fields.dat`:

| # | Kontrola | Co ověřuje | Práh |
|---|----------|------------|------|
| 1 | Bilance hmoty | $\sum \nabla\cdot\mathbf{u} \approx 0$ globálně | < 10⁻⁶ |
| 2 | Kontinuitní reziduum | max lokální $\|\nabla\cdot\mathbf{u}\|$ | < 10⁻⁴ |
| 3 | Symetrie kolem y=½ | $u_x(y) = u_x(L_y{-}y)$, $u_y$ antisym. | < 5·10⁻³ |
| 4 | Symetrie kolem z=½ | totéž v ose z | < 5·10⁻³ |
| 5 | Tlakový pokles | $\Delta p < 0$ (klesá po proudu) | — |
| 6 | No-slip | $\|U\|=0$ uvnitř krychle | < 10⁻¹⁰ |

### Proč právě tyto kontroly

- **Bilance + reziduum** — přímo testují, že SIMPLE splnil kontinuitu (jádro
  algoritmu). Hodnota ~10⁻¹⁶ = strojová přesnost.
- **Symetrie** — pro Re = 20 a centrovanou krychli *musí* být řešení symetrické
  ve y i z. Asymetrie by signalizovala chybu. Je to nezávislý test, který
  nevychází z žádné rovnice, jen z geometrie.
- **Tlakový pokles** — fyzikální zdravý rozum (tekutina teče z vyššího do
  nižšího tlaku).
- **No-slip** — ověří, že maska pevných buněk funguje (rychlost přesně 0).

## Grafy konvergence

`plot_convergence.py` čte `convergence.dat` a vykreslí historii reziduí:

- **Panel 1** — kontinuitní reziduum (log osa) s vyznačeným minimem
- **Panel 2** — `dU/U` a `dP/P` (změna polí mezi iteracemi)

Skript zároveň exportuje očištěná data zpět zpět do `convergence.dat`
(sloupce `iter cont dU/U dP/P`).

### Reziduum vs. změna polí — dvě různé věci

```
cont        = Σ|div U|        ... jak je porušena bilance hmoty
dU/U, dP/P  = ‖ΔU‖/‖U‖        ... jak moc se řešení ještě hýbe
```

Důležité při diagnostice stagnace: na jemné síti se kontinuitní reziduum ustálí
kolem 10⁻⁵, ale `dU/U` a `dP/P` můžou být mnohem nižší — pole je prakticky
zamrzlé, jen se kontinuita mírně přeskupuje přes Rhie–Chow korekci. Snapshoty
(`fields_0250/0300/0400.dat`) to potvrdí přímým porovnáním polí.

## Vizualizace polí

`plot_fields.py` (čte `fields.dat`) vykreslí řezy XY a XZ:

| Veličina | Co ukazuje |
|----------|------------|
| Rychlost $\|U\|$ | zpomalení před krychlí, urychlení po bocích, úplav za ní |
| Tlak $p$ | vysoký před překážkou, nízký za ní (tlakový odpor) |
| Proudnice | obtékání + symetrická recirkulační zóna |

`plot_mesh.py` (čte `mesh.dat`) zobrazí 3D schéma tunelu, polohu řezů a detail
sítě kolem krychle.

## Vizualizace matic

`plot_matrix.py` (čte `pmatrix.dat`, `umatrix.dat`) ukáže:

- **Sparsity pattern** — tři diagonální pásy, identity řádky pevných buněk
- **Prostorové rozložení koeficientů** — aP, aE, aN po doméně
- **Diagonální dominance** — |aP|/Σ|off| (stabilita)
- **Hustý blok** kolem krychle — 7bodový stencil zblízka

## ParaView

Pro plné 3D zkoumání otevři `channel_cube.vtk` v ParaView — řezy v libovolné
rovině, izoplochy rychlosti, sledování částic (streamlines), filtrování krychle
přes pole `solid`.

<div class="admonition tip">
<p class="admonition-title">Typický pracovní postup</p>
<p>1. <code>make &amp;&amp; OMP_NUM_THREADS=6 ./channelCube</code><br>
2. <code>python verify.py</code> — všechny kontroly musí projít<br>
3. <code>python plot_convergence.py</code> — zkontroluj konvergenci<br>
4. <code>python plot_fields.py</code> — podívej se na proudění<br>
5. ParaView <code>channel_cube.vtk</code> — 3D detail</p>
</div>
