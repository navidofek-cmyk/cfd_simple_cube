# Použité prvky C++ — a proč

Tato stránka odpovídá na otázku **„proč je tam co je"**. Každý jazykový prvek
byl zvolen z konkrétního důvodu — ne proto, že „se to tak dělá", ale proto, že
řeší konkrétní problém řešiče.

## Proč vůbec OOP (objekty, struktury)?

Mohli bychom napsat čisté C — globální pole a funkce. Proč třídy?

**Důvod: stav, který patří k sobě.** Řešič má desítky polí (rychlost, tlak,
toky, koeficienty, workspace) a parametry (Re, viskozita, podrelaxace). Bez
objektů by se to předávalo jako dvacet argumentů do každé funkce, nebo by to
byly globální proměnné (nebezpečné, neznovupoužitelné).

`ChannelSolver` **sdružuje stav a operace nad ním**:

```cpp
struct ChannelSolver {
    Mesh m; double Re, nu; volVector U; volScalar p; ...   // stav
    void solve();                                           // operace
};
```

Důsledek: spustit dva výpočty s různým Re = vytvořit dva objekty. Stav každého
je izolovaný. To by s globály nešlo.

<div class="admonition note">
<p class="admonition-title">struct vs class</p>
<p>Používáme <code>struct</code> (členy veřejné) místo <code>class</code>
(soukromé). Proč? Tohle je numerický kód, ne knihovna pro cizí uživatele.
Zapouzdření s gettery/settery by byla režie bez užitku — vizualizace a main
potřebují přímý přístup k polím (<code>sim.U</code>, <code>sim.p</code>).
Jednoduchost &gt; dogma.</p>
</div>

## Proč přetěžování operátorů?

`Vec3` má `operator+`, `operator-`, `operator*`. Proč ne funkce `add(a,b)`?

**Důvod: čitelnost matematického kódu.** Porovnej korekci rychlosti:

```cpp
U[c] = HbyA[c] - rAU[c]*gpn[c];              // s operátory
U[c] = sub(HbyA[c], scale(rAU[c], gpn[c]));  // bez operátorů
```

První řádek je doslova fyzikální vzorec $\mathbf{u} = \mathbf{HbyA} - r_{AU}\nabla p$.
Druhý je hlavolam. V CFD kódu, kde je vzorců stovky, je čitelnost zásadní —
chyba znaménka se v prvním zápisu pozná na první pohled.

**Proč `double * Vec3` a ne `Vec3 * double`?** Protože vzorce píšeme „skalár krát
vektor" ($r_{AU}\,\mathbf{H}$). Operátor je definovaný v tom pořadí, v jakém ho
používáme.

## Proč šablony (templates)?

`fvMatrix<T>` a `Field<T>` jsou šablonované. Proč ne dvě kopie kódu — jedna pro
tlak, jedna pro rychlost?

**Důvod: jeden algoritmus, dva typy dat.** Sestavení matice, Jacobi sweep,
výpočet rezidua — to je **identický kód** pro skalární tlak i vektorovou
rychlost. Liší se jen typ pravé strany (`double` vs `Vec3`). Šablona napíše
algoritmus jednou:

```cpp
template<class T>
struct fvMatrix {
    std::vector<double> aP, aE, ...;   // koeficienty VŽDY skalární
    std::vector<T>      src;           // pravá strana: double NEBO Vec3
};
```

Bez šablon bychom měli `fvMatrixScalar` a `fvMatrixVector` se zduplikovaným
kódem — dvojnásobek chyb, dvojnásobek údržby. Šablona je **bezpečnější než
kopie** a **rychlejší než dědičnost** (žádné virtuální volání, vše se vyřeší
při překladu).

**Co šablonu umožní:** přetížené `mag2(double)` i `mag2(Vec3)` (z
[Vec3.hpp](/page/cpp/01-vec3)). Šablonový kód volá `mag2(np - psi[c])` a
překladač podle `T` vybere správnou verzi.

## Proč reference (`&`) místo kopií a ukazatelů?

Skoro všude se předává `const Mesh& m`, `Field<T>& psi`. Proč?

**Důvod 1 — výkon:** `Mesh` a `Field` jsou velké (miliony prvků). Předání
hodnotou by je kopírovalo. Reference předá jen adresu.

**Důvod 2 — bezpečnost oproti ukazatelům:** reference nemůže být `null` a nejde
„přesměrovat". `const Mesh&` navíc překladačem vynutí, že funkce síť nezmění:

```cpp
static volVector gradP(const Mesh& m, const volScalar& p)  // ani m ani p nezměním
```

Kdybych se v těle pokusil `m` změnit, překlad selže. `const` je **kontrakt
vynucený překladačem**, ne jen komentář.

## Proč `auto`?

```cpp
auto UEqn = buildUEqn(...);   // místo fvMatrix<Vec3> UEqn = ...
```

**Důvod: typ je dlouhý a zřejmý z pravé strany.** `auto` ušetří psaní a hlavně
údržbu — když změním návratový typ `buildUEqn`, nemusím opravovat každé volání.
Nepoužívá se tam, kde by typ nebyl jasný (parametry funkcí jsou vždy explicitní).

## Proč lambda výrazy?

```cpp
auto rnd = [](double v){ return (int)(v+0.5); };          // v Mesh
auto negAv = [&](const auto& v, auto& y){ ... };           // v CgSolver
```

**Důvod: krátká lokální logika bez zaplevelení jmenného prostoru.** `rnd` se
použije třikrát v jedné metodě — nemá smysl z toho dělat globální funkci.
`negAv` v CG zachytí lokální `m`, `M` přes `[&]` a volá se třikrát uvnitř
algoritmu. Bez lambdy bych musel ty proměnné předávat jako argumenty pořád
dokola.

## Proč RAII (std::vector, ofstream)?

**Důvod: žádné úniky paměti, žádné zapomenuté `fclose`.** `std::vector` se sám
uvolní, `std::ofstream` se sám zavře — když objekt zanikne (konec funkce, i při
výjimce). Porovnej:

```cpp
// C styl: musím nezapomenout uvolnit, i při chybě
double* arr = malloc(N*sizeof(double));
... if(chyba) { free(arr); return; }   // snadno zapomenu
free(arr);

// C++ RAII: uvolní se samo
std::vector<double> arr(N);
... if(chyba) return;                  // arr se uvolní automaticky
```

V solveru má `CgWorkspace` pět vektorů — RAII zaručí, že zaniknou se solverem.

## Proč inicializační seznamy konstruktorů?

```cpp
Mesh(int nx_, ...) : nx(nx_), dx(Lx/nx_), V(dx*dy*dz), ... {}
```

**Důvod: členy se inicializují přímo, ne dvoufázově.** Alternativa (přiřazení
v těle) by členy nejdřív default-zkonstruovala a pak přepsala — dvojí práce.
Pro `const` členy a reference je inicializační seznam dokonce **jediná možnost**.
Pozor: členy se inicializují v pořadí *deklarace*, na čemž tady závisí výpočet
`V` z `dx,dy,dz`.

## Proč OpenMP pragmy a ne std::thread?

```cpp
#pragma omp parallel for collapse(3)
for(...) { ... }
```

**Důvod: datová paralelizace milionů nezávislých buněk.** `std::thread` by
znamenal ruční rozdělení rozsahu, ruční join, ruční redukce — desítky řádků
boilerplate na každou smyčku. OpenMP to udělá jedním řádkem `#pragma`. Pro
„udělej tuhle smyčku paralelně" je OpenMP přesně ten nástroj. `std::thread`/
`std::async` by se hodily na *úlohovou* paralelizaci (různé úlohy najednou), ne
na rozsekání jedné smyčky. Viz [Paralelizace](/page/theory/09-parallel).

## Proč `static` funkce v Simple.cpp?

```cpp
static fvMatrix<Vec3> buildUEqn(...) { ... }
```

**Důvod: jsou to interní pomocníci, ne veřejné API.** `static` u funkce na
úrovni souboru = **interní vazba** (viditelná jen v tomto .cpp). Nezanáší
globální jmenný prostor a dává překladači volnost je agresivně inlinovat. Veřejné
je jen `ChannelSolver::solve()`.

## Shrnutí: každý prvek řeší konkrétní problém

| Prvek | Řeší problém |
|-------|--------------|
| OOP / struct | sdružení stavu solveru, izolace instancí |
| přetěžování operátorů | čitelnost fyzikálních vzorců |
| šablony | jeden algoritmus pro skalár i vektor |
| reference + const | výkon + překladačem vynucená bezpečnost |
| auto | méně psaní, snazší údržba |
| lambda | lokální logika bez globálních funkcí |
| RAII | automatická správa paměti a souborů |
| init seznamy | přímá, jednofázová inicializace |
| OpenMP pragma | datová paralelizace jedním řádkem |
| static funkce | interní pomocníci mimo veřejné API |
