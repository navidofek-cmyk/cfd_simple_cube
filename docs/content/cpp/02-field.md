# Field.hpp — šablony a kontejnery

`Field<T>` je tenká obálka nad `std::vector<T>` reprezentující buňkové pole.
Ukazuje **šablony tříd**, `std::vector` a **type aliasy**.

## Celý soubor

```cpp
#pragma once
#include <vector>
#include "Vec3.hpp"

template<class T>
struct Field {
    std::vector<T> d;
    Field() = default;
    explicit Field(int n, T v = T{}) : d(n, v) {}
    T&       operator[](int i)       { return d[i]; }
    const T& operator[](int i) const { return d[i]; }
    int size() const { return (int)d.size(); }
};

using volScalar = Field<double>;
using volVector = Field<Vec3>;
```

## Jazykové prvky

### Šablona třídy

```cpp
template<class T>
struct Field { std::vector<T> d; ... };
```

`T` je typový parametr. Z jedné definice překladač vygeneruje konkrétní typy
podle použití — `Field<double>` a `Field<Vec3>`. Tomu se říká **generické
programování**: kód napsaný jednou, funguje pro libovolný typ s potřebnými
operacemi.

### `std::vector<T>`

Dynamické pole z STL. Samo spravuje paměť (alokace, dealokace, růst) —
**RAII**: když `Field` zanikne, `vector` automaticky uvolní svou paměť. Žádné
ruční `new`/`delete`.

### `explicit` konstruktor

```cpp
explicit Field(int n, T v = T{}) : d(n, v) {}
```

- `explicit` zabrání nechtěné implicitní konverzi — `Field<double> f = 5;` se
  nepřeloží, musí se `Field<double> f(5);`. Brání chybám.
- `T v = T{}` — **výchozí argument** s **value-initialization** `T{}`
  (`0.0` pro double, `{0,0,0}` pro Vec3).
- `: d(n, v)` — inicializační seznam vytvoří vektor délky `n` vyplněný hodnotou
  `v`.

### Dvě verze `operator[]`

```cpp
T&       operator[](int i)       { return d[i]; }   // pro úpravu
const T& operator[](int i) const { return d[i]; }   // pro čtení z const objektu
```

**Const overloading**: nekonstantní verze vrací referenci pro zápis
(`U[c] = ...`), konstantní verze (volaná na `const Field&`) vrací
`const T&` jen pro čtení. Umožňuje předávat pole jako `const Field<T>&` tam, kde
je nemáme měnit (např. v `interpFlux`).

### Type aliasy

```cpp
using volScalar = Field<double>;   // skalární pole (tlak, reziduum)
using volVector = Field<Vec3>;     // vektorové pole (rychlost)
```

`using` (moderní náhrada za `typedef`) dává čitelná jména. Název `vol*`
navazuje na konvenci OpenFOAM (`volScalarField`, `volVectorField`) — „vol" =
veličina definovaná na objemech (buňkách).

<div class="admonition note">
<p class="admonition-title">Proč obálka, ne přímo std::vector</p>
<p>Field by mohl být jen <code>std::vector</code>, ale tenká obálka dává:
(1) sémantický typ <code>volScalar</code>/<code>volVector</code>,
(2) <code>int size()</code> místo <code>size_t</code> (méně varování při
porovnání se znaménkovými indexy ve smyčkách),
(3) prostor pro budoucí rozšíření (např. okrajové hodnoty) beze změny volajícího
kódu.</p>
</div>
