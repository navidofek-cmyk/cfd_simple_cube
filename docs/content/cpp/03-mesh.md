# Mesh.hpp — síť a maska

`Mesh` drží geometrii uniformní 3D mřížky a masku pevných buněk (krychle).

## Celý soubor

```cpp
#pragma once
#include <vector>

struct Mesh {
    int nx, ny, nz;
    double dx, dy, dz, V;
    double Ax, Ay, Az;           // plochy stěn (E/W, N/S, F/B)
    std::vector<bool> solid;

    Mesh(int nx_, int ny_, int nz_, double Lx, double Ly, double Lz)
        : nx(nx_), ny(ny_), nz(nz_),
          dx(Lx/nx_), dy(Ly/ny_), dz(Lz/nz_),
          V(dx*dy*dz),
          Ax(dy*dz), Ay(dx*dz), Az(dx*dy),
          solid(nx_*ny_*nz_, false) {}

    int N() const { return nx*ny*nz; }
    int id(int i, int j, int k) const { return i + j*nx + k*nx*ny; }

    void markSolid(double x0, double x1, double y0, double y1,
                   double z0, double z1) {
        auto rnd=[](double v){ return (int)(v+0.5); };
        int i0=rnd(x0/dx), i1=rnd(x1/dx);
        int j0=rnd(y0/dy), j1=rnd(y1/dy);
        int k0=rnd(z0/dz), k1=rnd(z1/dz);
        for (int k=k0; k<k1; ++k)
            for (int j=j0; j<j1; ++j)
                for (int i=i0; i<i1; ++i)
                    solid[id(i,j,k)] = true;
    }
};
```

## Jazykové prvky

### Inicializační seznam s odvozenými hodnotami

```cpp
Mesh(int nx_, ...) : nx(nx_), dx(Lx/nx_), V(dx*dy*dz), Ax(dy*dz), ...
```

Členy se inicializují **v pořadí deklarace** (ne v pořadí v seznamu!). Proto je
důležité, že `dx, dy, dz` jsou deklarovány *před* `V` a `Ax` — když se počítá
`V(dx*dy*dz)`, `dx` už existuje. Toto je častý zdroj chyb; tady je pořadí
správné.

### `std::vector<bool>` — speciální případ

```cpp
std::vector<bool> solid;
```

`vector<bool>` je v STL **bitově komprimovaný** (1 bit na prvek, ne 1 byte) —
2 mil. buněk = 250 kB místo 2 MB. Cena: nelze brát adresu prvku (`&solid[i]`
nefunguje), ale tady to nevadí, jen čteme/zapisujeme hodnoty.

### `const` členské metody

```cpp
int N() const { return nx*ny*nz; }
int id(int i,int j,int k) const { ... }
```

`const` za hlavičkou slibuje, že metoda objekt nemění. Umožňuje volat je na
`const Mesh&` — a `Mesh` se všude předává jako `const Mesh&` (solver síť nemění).

### Linearizace 3D indexu

```cpp
int id(int i, int j, int k) const { return i + j*nx + k*nx*ny; }
```

Převod 3D souřadnic na 1D index do plochého pole. **x se mění nejrychleji**
(řádkové uspořádání) → sousedé v ose x jsou v paměti vedle sebe (`id±1`), což je
cache-friendly pro nejčastější přístupy.

| Soused | Posun indexu |
|--------|--------------|
| východ/západ (±x) | ±1 |
| sever/jih (±y) | ±nx |
| předek/zadek (±z) | ±nx·ny |

### Lambda uvnitř metody

```cpp
auto rnd=[](double v){ return (int)(v+0.5); };
```

**Lambda výraz** — anonymní funkce definovaná lokálně. `rnd` zaokrouhlí na
nejbližší celé číslo. Řeší konkrétní chybu: `(int)(0.6/0.05)` dá kvůli
zaobrouhlovací chybě plovoucí čárky `11` místo `12`. Přidání `+0.5` před
oříznutím to opraví → krychle má správně 16 buněk na stranu, ne 15.

<div class="admonition warning">
<p class="admonition-title">Plovoucí čárka a markSolid</p>
<p>Bez zaokrouhlení dávalo <code>markSolid(0.9, 1.1, ...)</code> krychli
4×3×3 místo 4×4×4 buněk — hranice 0.6/0.05 se uložila jako 11.9999 a oříznutí
udělalo 11. Lambda <code>rnd</code> je oprava této reálné chyby.</p>
</div>

## Role masky `solid`

Krychle není z mřížky vyříznutá — síť zůstává plně strukturovaná (rychlá
indexace `id`). Místo toho jsou buňky uvnitř krychle označené `solid[c]=true` a
solver je v matici nahradí identitou. Tomu se říká **immersed boundary /
masking** přístup. Detaily v [Okrajové podmínky](/page/theory/06-boundary).
