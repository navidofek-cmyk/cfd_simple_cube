# fvMatrix.hpp — matice a stencil

`fvMatrix<T>` je srdce diskretizace: drží 7bodový stencil jako sedm polí
koeficientů a umí Jacobi řešení. Šablona funguje pro tlak (`double`) i rychlost
(`Vec3`).

## Struktura

```cpp
template<class T>
struct fvMatrix {
    const Mesh& m;
    Field<T>& psi;                                 // řešené pole (reference!)
    std::vector<double> aP, aE, aW, aN, aS, aF, aB; // 7 diagonál
    std::vector<T>      src;                        // pravá strana

    fvMatrix(const Mesh& m_, Field<T>& psi_)
        : m(m_), psi(psi_),
          aP(m_.N(),0.0), aE(m_.N(),0.0), /* ... */
          src(m_.N(), T{}) {}
    ...
};
```

## Jazykové prvky

### Reference jako členy

```cpp
const Mesh& m;
Field<T>& psi;
```

- `m` je **konstantní reference** na síť — matice síť nevlastní ani nekopíruje,
  jen na ni odkazuje (síť žije v solveru).
- `psi` je **nekonstantní reference** na řešené pole — Jacobi sweep zapisuje
  přímo do něj.

<div class="admonition warning">
<p class="admonition-title">Reference členy a životnost</p>
<p>Reference se musí inicializovat v konstruktoru (nelze později „přesměrovat")
a objekt, na který ukazují, musí žít déle než <code>fvMatrix</code>. Tady to
platí — matice je dočasná, síť i pole žijí v solveru.</p>
</div>

### Koeficienty vždy `double`, zdroj typu `T`

```cpp
std::vector<double> aP, aE, ...;   // skalární koeficienty
std::vector<T>      src;           // T = double NEBO Vec3
```

Klíčový trik šablony: **koeficienty matice jsou skalární i pro vektorovou
rovnici**. Hybnostní rovnice pro `u, v, w` má stejné koeficienty (stejná
geometrie, stejné toky), liší se jen pravá strana. Proto `src` je typu `T`, ale
`aP...aB` jsou vždy `double`. Tím jedna šablona pokryje obojí.

## Jacobi řešič

```cpp
void solveWithExtra(const Field<T>& extra, double scale, int sweeps) {
    std::vector<T> tmp(m.N());
    for(int s=0;s<sweeps;++s){
        double res=0.0;
        #pragma omp parallel for collapse(3) reduction(+:res)
        for(int k=0;k<m.nz;++k) for(int j=0;j<m.ny;++j) for(int i=0;i<m.nx;++i){
            int c=m.id(i,j,k);
            T rhs=src[c];
            if(!m.solid[c]) rhs=rhs+(scale*extra[c]);
            if(i<m.nx-1) rhs=rhs-(aE[c]*psi[m.id(i+1,j,k)]);
            /* ... ostatní sousedé ... */
            T np=(1.0/aP[c])*rhs;
            res+=mag2(np-psi[c]);     // ← mag2 přetížené pro double i Vec3
            tmp[c]=np;
        }
        std::swap(psi.d, tmp);        // ← prohození místo kopie
        if(s==0) res0=std::sqrt(res)+1e-30;
        if(std::sqrt(res)/res0<1e-3) break;
    }
}
```

### Jacobi vs. Gauss–Seidel

Jacobi čte hodnoty `psi` z **minulého** sweepu (zapisuje do odděleného `tmp`),
takže iterace jsou **nezávislé → paralelizovatelné**. Gauss–Seidel by četl
právě aktualizované hodnoty (sériový). Viz [Paralelizace](/page/theory/09-parallel).

### `std::swap` místo kopie

```cpp
std::swap(psi.d, tmp);
```

Prohodí vnitřní buffery dvou vektorů v konstantním čase (jen přehodí ukazatele),
místo kopírování milionů prvků. Po swapu `psi` obsahuje nové hodnoty.

### `extra` a `scale` — tlakový gradient bez kopie matice

```cpp
UEqn.solveWithExtra(gp, -m.V, 40);   // řeší A·u = src − V·∇p
```

Místo aby se do `src` napevno přičetl tlakový gradient (a matice se musela
kopírovat mezi iteracemi SIMPLE), předá se gradient `extra` zvlášť a přičte se
za běhu (`rhs += scale*extra[c]`). Šetří kopii velké matice.

## Pomocné metody pro SIMPLE

```cpp
volScalar A() const;   // a_P / V  → koeficient pro HbyA
Field<T>  H() const;   // (src − Σ a_nb·ψ_nb) / V  → operátor H
void relax(double alpha);  // implicitní podrelaxace: aP/=α, src+=(aP_new−aP_old)·ψ
```

Tyto metody implementují členy z odvození SIMPLE (viz
[SIMPLE](/page/theory/04-simple)). `A()` a `H()` se používají k sestavení
`HbyA = H/A` a tlakové rovnice. `relax()` provádí implicitní podrelaxaci
rychlosti zvýšením diagonály.

### Pevné buňky v `relax` a `A`

```cpp
void relax(double alpha){
    for(int c=0;c<m.N();++c){
        if(m.solid[c]) continue;        // pevné buňky se nepodrelaxují
        double aOld=aP[c]; aP[c]/=alpha;
        src[c]+=(aP[c]-aOld)*psi[c];
    }
}
```

Pevné buňky mají `aP=1` (identita) a přeskakují se, aby si identitu udržely.
