# Simple.cpp — jádro řešiče

Nejdůležitější soubor. Obsahuje sestavení obou matic a hlavní SIMPLE smyčku.
Implementuje teorii z kapitol [SIMPLE](/page/theory/04-simple) a
[Rhie–Chow](/page/theory/05-rhie-chow).

## Statické pomocné funkce

Soubor obsahuje pět `static` funkcí (viditelné jen v tomto .cpp):

| Funkce | Vrací | Účel |
|--------|-------|------|
| `buildUEqn` | `fvMatrix<Vec3>` | sestaví hybnostní matici |
| `buildPEqn` | `fvMatrix<double>` | sestaví tlakovou matici |
| `interpFlux` | (přes reference) | toky z HbyA na stěny |
| `divFlux` | `volScalar` | divergence toků (reziduum, RHS) |
| `gradP` | `volVector` | gradient tlaku ve středech buněk |

`static` u funkce na úrovni souboru znamená **interní vazbu** — funkce není
viditelná z jiných .cpp, nezanáší globální jmenný prostor.

## buildUEqn — buňkový přístup

Klíčová funkce. Sestaví hybnostní rovnici tak, aby šla paralelizovat:

```cpp
static fvMatrix<Vec3> buildUEqn(const Mesh& m,
                                 const volScalar& phiE, const volScalar& phiN,
                                 const volScalar& phiF,
                                 double nu, double U_in, volVector& U)
{
    fvMatrix<Vec3> eq(m, U);
    const double gx = nu*m.Ax/m.dx, /* ... */;
    const double gbx = nu*m.Ax/(0.5*m.dx);   // hraniční (půl buňky)
    const Vec3 Uin{U_in, 0.0, 0.0};

    #pragma omp parallel for collapse(3)
    for(/* všechny buňky */){
        int c = m.id(i,j,k);
        if (m.solid[c]){ eq.aP[c]=1.0; continue; }   // pevná buňka = identita

        // ── Západ ──
        if (i>0){ int w=m.id(i-1,j,k);
            if(!m.solid[w]){
                eq.aW[c]-=gx; eq.aP[c]+=gx;          // difúze
                double f=-phiE[w];                    // tok přes záp. stěnu
                if(f>=0) eq.aP[c]+=f; else eq.aW[c]+=f; // upwind konvekce
            } else { eq.aP[c]+=gbx; }                 // pevný soused = no-slip
        } else {                                       // i==0: INLET
            eq.aP[c]  += gbx;
            eq.src[c] += gbx*Uin;                     // difúze ze vstupu
            eq.src[c] += U_in*m.Ax*Uin;               // konvekce ze vstupu
        }
        // ── Východ, Jih, Sever, Zadek, Předek analogicky ──
    }
    return eq;
}
```

### Proč buňkový přístup

Každá buňka `c` zapisuje **jen do svých** koeficientů (`eq.aW[c]`, `eq.aP[c]`),
sousedy jen **čte** (`phiE[w]`, `m.solid[w]`). Tok přes západní stěnu si přečte
ze sousedovy hodnoty `phiE[w]` se záporným znaménkem (`-phiE[w]` = odtok na
západ). Žádné zápisy do sousedů → **žádný race condition** → bezpečné pro
`#pragma omp parallel for`. Viz [Konvekce](/page/theory/03-convection).

### Znaménková konvence

```cpp
eq.aW[c] -= gx;   // mimodiagonála: záporná
eq.aP[c] += gx;   // diagonála: kladná
```

Rovnice je `div − ν·laplacian`. Laplacián přispívá `+gx` do diagonály (kladné)
a `−gx` mimo. Tato konvence zaručuje diagonální dominanci. (Špatné znaménko zde
bylo původní chybou, která dávala NaN.)

### `return eq` — návrat velkého objektu

Funkce vrací `fvMatrix<Vec3>` hodnotou. Díky **move semantics** / **NRVO**
(named return value optimization) se nekopíruje — překladač sestaví objekt
rovnou na místě volajícího.

## Hlavní smyčka solve()

```cpp
void ChannelSolver::solve(){
    volVector Uprev = U;
    for(int it=0; it<maxIter; ++it){
        Uprev = U;

        // 1. Hybnost → U*
        auto UEqn = buildUEqn(m, phiE,phiN,phiF, nu,U_in, U);
        UEqn.relax(alphaU);
        auto gp = gradP(m, p);
        UEqn.solveWithExtra(gp, -m.V, 40);          // Jacobi → U*

        // 2. A, HbyA
        volScalar A = UEqn.A();
        volVector Hf = UEqn.H();
        volScalar rAU(m.N(),1.0);
        volVector HbyA(m.N());
        #pragma omp parallel for
        for(int c=0;c<m.N();++c){
            if(m.solid[c]){ rAU[c]=0.0; continue; }
            rAU[c]=1.0/A[c]; HbyA[c]=rAU[c]*Hf[c];
        }

        // 3.-4. toky, tlaková rovnice
        volScalar pE(m.N()),pN(m.N()),pF2(m.N());
        interpFlux(m, HbyA, pE,pN,pF2);
        auto pEqn = buildPEqn(m, rAU, p);
        auto divH = divFlux(m, pE,pN,pF2, U_in);
        #pragma omp parallel for
        for(int c=0;c<m.N();++c) pEqn.src[c]+=divH[c];
        cgSolve(pEqn, cgWork, 1000, 1e-6);          // CG → p

        // 5. Rhie–Chow korekce toků
        // 6. korekce rychlosti U = HbyA − rAU·∇p
        // 7. podrelaxace tlaku
        // 8. reziduum + ΔU/Δp + snapshoty
    }
}
```

### `auto` pro odvození typu

```cpp
auto UEqn = buildUEqn(...);   // typ fvMatrix<Vec3> odvozen
auto gp   = gradP(m, p);       // volVector
```

`auto` nechá překladač odvodit typ z pravé strany — méně psaní, bezpečné.

### `ChannelSolver::solve()` — definice členské metody

`::` je **scope resolution operator** — definuje metodu `solve` patřící třídě
`ChannelSolver` (deklarované v [Simple.hpp](/page/cpp/00-overview)). Tělo je
mimo definici třídy, v .cpp.

## Diagnostika: ΔU/Δp a snapshoty

Smyčka navíc počítá relativní L2 změnu polí mezi iteracemi a umí ukládat
snapshoty:

```cpp
#pragma omp parallel for reduction(+:cont,dU2,dP2,U2,P2)
for(int c=0;c<m.N();++c){
    cont += std::fabs(resField[c]);          // kontinuitní reziduum
    dU2  += mag2(U[c]-Uprev[c]);             // změna rychlosti
    dP2  += (p[c]-pOld[c])*(p[c]-pOld[c]);   // změna tlaku
}
```

Více proměnných v jedné `reduction` klauzuli. Výstup se zapisuje živě do
`convergence.dat`. Viz [Validace](/page/solver/03-validation).
