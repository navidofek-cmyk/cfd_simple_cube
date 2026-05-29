# CgSolver.hpp — lambda a workspace

Conjugate Gradient řešič tlakové rovnice. Ukazuje **lambda s capture**,
**workspace vzor** a paralelní redukce.

## Workspace

```cpp
struct CgWorkspace {
    std::vector<double> b, Ax, r, dir, Ap;
    void ensure(int N) {
        if(static_cast<int>(b.size())==N) return;
        b.resize(N); Ax.resize(N); r.resize(N); dir.resize(N); Ap.resize(N);
    }
};
```

CG potřebuje 5 pomocných vektorů. `CgWorkspace` je drží a `ensure(N)` je
**jednou** naalokuje (při změně velikosti). Solver ho drží jako member a předává
do každého z ~1000 volání `cgSolve` za běh — žádné opakované alokace.

`static_cast<int>` je explicitní, typově bezpečný převod (`b.size()` vrací
`size_t`) — čistší než C-style `(int)`.

## Funkce cgSolve

```cpp
inline void cgSolve(fvMatrix<double>& M, CgWorkspace& ws,
                    int maxIter=1000, double tol=1e-6) {
    const Mesh& m = M.m;
    const int   N = m.N();
    auto& x = M.psi.d;          // reference na řešené pole tlaku
    ws.ensure(N);

    auto negAv = [&](const std::vector<double>& v, std::vector<double>& y) {
        #pragma omp parallel for collapse(3)
        for(/* všechny buňky */){
            int c=m.id(i,j,k);
            y[c] = -M.aP[c]*v[c];
            if(i<m.nx-1) y[c] -= M.aE[c]*v[m.id(i+1,j,k)];
            /* ... ostatních 5 sousedů ... */
        }
    };
    ...
}
```

## Jazykové prvky

### Lambda s capture `[&]`

```cpp
auto negAv = [&](const std::vector<double>& v, std::vector<double>& y) { ... };
```

`negAv` je lokální funkce počítající $\mathbf{y} = (-\mathbf{A})\mathbf{v}$.

- `[&]` = **capture by reference** — lambda vidí všechny lokální proměnné
  (`m`, `M`) odkazem, bez kopírování. Může je číst přímo.
- `auto` — typ lambdy je anonymní, překladač ho odvodí.
- Volá se pak `negAv(ws.dir, ws.Ap)` — třikrát v algoritmu, bez duplikace kódu.

### Maticový součin bez matice

```cpp
y[c] = -M.aP[c]*v[c];
if(i<m.nx-1) y[c] -= M.aE[c]*v[m.id(i+1,j,k)];
```

Násobení $(-\mathbf{A})\mathbf{v}$ se dělá procházením 7 diagonál a sousedů —
matice nikdy neexistuje jako 2D struktura. To je smysl uložení po diagonálách
(viz [Řídké matice](/page/theory/07-linear-system)). Znaménko `-` realizuje obrat
na pozitivně definitní soustavu.

### Reference na vnitřek pole

```cpp
auto& x = M.psi.d;   // přímý přístup k std::vector uvnitř Field
```

`x` je alias pro vektor hodnot tlaku — CG zapisuje řešení přímo do pole `p`
solveru.

## Paralelní redukce

```cpp
double rr=0.0;
#pragma omp parallel for reduction(+:rr)
for(int c=0;c<N;++c){ ws.r[c]=ws.b[c]-ws.Ax[c]; rr+=ws.r[c]*ws.r[c]; }
```

Skalární součin $\mathbf{r}^T\mathbf{r}$ paralelně — `reduction(+:rr)` dá každému
vláknu privátní `rr` a na konci je sečte. Stejně pro $\mathbf{d}^T\mathbf{Ap}$.

## Tělo iterace

```cpp
for(int k=0;k<maxIter;++k){
    negAv(ws.dir, ws.Ap);                         // Ap = (-A)·dir
    double pAp=0.0;
    #pragma omp parallel for reduction(+:pAp)
    for(int c=0;c<N;++c) pAp+=ws.dir[c]*ws.Ap[c];
    double alpha=rr/(pAp+1e-30);                   // délka kroku

    #pragma omp parallel for
    for(int c=0;c<N;++c){ x[c]+=alpha*ws.dir[c]; ws.r[c]-=alpha*ws.Ap[c]; }

    double rr_new=0.0;
    #pragma omp parallel for reduction(+:rr_new)
    for(int c=0;c<N;++c) rr_new+=ws.r[c]*ws.r[c];
    if(std::sqrt(rr_new)/res0<tol) break;          // konvergence

    double beta=rr_new/(rr+1e-30);
    #pragma omp parallel for
    for(int c=0;c<N;++c) ws.dir[c]=ws.r[c]+beta*ws.dir[c];
    rr=rr_new;
}
```

Přesně mapuje vzorce z [teorie CG](/page/theory/08-cg). Všechny smyčky jsou
paralelní. `+1e-30` ve jmenovatelích brání dělení nulou (když reziduum padne na
nulu přesně).

<div class="admonition note">
<p class="admonition-title">inline funkce v hlavičce</p>
<p><code>cgSolve</code> je <code>inline</code> a celá v hlavičce — může se
vložit do více .cpp bez chyby linkeru a překladač ji může optimalizovat
v místě volání.</p>
</div>
