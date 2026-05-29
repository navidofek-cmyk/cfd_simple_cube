# SIMPLE smyčka krok po kroku

Tato stránka spojuje teorii a kód — projde jednu iteraci `ChannelSolver::solve()`
a ukáže, který řádek dělá který krok algoritmu.

## Schéma jedné iterace

```
        ┌─────────────────────────────────────────────┐
        │  Uprev = U          (záloha pro ΔU)          │
        ├─────────────────────────────────────────────┤
   1.   │  buildUEqn  →  fvMatrix<Vec3>  (hybnost)     │
   2.   │  UEqn.relax(αU)                              │
   3.   │  solveWithExtra(∇p, −V, 40)  →  U*  (Jacobi) │
        ├─────────────────────────────────────────────┤
   4.   │  A = UEqn.A();  HbyA = H/A                   │
   5.   │  interpFlux  →  toky pE,pN,pF z HbyA         │
        ├─────────────────────────────────────────────┤
   6.   │  buildPEqn  →  fvMatrix<double>  (tlak)      │
   7.   │  cgSolve  →  p   (Conjugate Gradient)        │
        ├─────────────────────────────────────────────┤
   8.   │  Rhie–Chow korekce toků phiE,phiN,phiF       │
   9.   │  U = HbyA − rAU·∇p   (korekce rychlosti)     │
  10.   │  p = pOld + αP·(p − pOld)  (podrelaxace)      │
        ├─────────────────────────────────────────────┤
  11.   │  reziduum + ΔU/ΔP → convergence.dat          │
        └─────────────────────────────────────────────┘
                         opakuj dokud cont/cont₀ > tol
```

## Krok za krokem

### 1–3. Prediktor rychlosti

```cpp
auto UEqn = buildUEqn(m, phiE,phiN,phiF, nu,U_in, U);  // sestav matici
UEqn.relax(alphaU);                                     // implicitní podrelaxace
auto gp = gradP(m, p);                                  // gradient STARÉHO tlaku
UEqn.solveWithExtra(gp, -m.V, 40);                      // vyřeš → U*
```

Vyřeší hybnost se starým tlakem. Výsledek `U*` obecně **nesplňuje kontinuitu** —
to napraví tlaková korekce. Matematicky viz [SIMPLE](/page/theory/04-simple),
kód [Simple.cpp](/page/cpp/06-simple).

### 4. Operátory A a HbyA

```cpp
volScalar A = UEqn.A();      // a_P / V
volVector Hf = UEqn.H();     // (src − Σ a_nb·u_nb)/V
for(int c=0;c<m.N();++c){
    rAU[c]=1.0/A[c];         // 1/a_P
    HbyA[c]=rAU[c]*Hf[c];    // HbyA = H/A
}
```

`HbyA` je rychlost vyjádřená bez tlakového členu — vstup do tlakové rovnice.

### 5. Toky z HbyA (první půlka Rhie–Chow)

```cpp
interpFlux(m, HbyA, pE,pN,pF2);   // interpolace HbyA na stěny
```

### 6–7. Tlaková rovnice

```cpp
auto pEqn = buildPEqn(m, rAU, p);          // Laplaceův operátor
auto divH = divFlux(m, pE,pN,pF2, U_in);   // pravá strana = div(HbyA)
for(int c=0;c<m.N();++c) pEqn.src[c]+=divH[c];
cgSolve(pEqn, cgWork, 1000, 1e-6);          // vyřeš tlak metodou CG
```

Nejdražší krok — tlaková Poissonova rovnice se řeší přesně každou iteraci přes
[Conjugate Gradient](/page/theory/08-cg).

### 8. Rhie–Chow korekce toků (druhá půlka)

```cpp
phiE[c] = pE[c] - rf(c,e)*m.Ax/m.dx*(p[e]-p[c]);
```

Toky se opraví o tlakový gradient mezi sousedy → divergence-free. Viz
[Rhie–Chow](/page/theory/05-rhie-chow). Tyto `phiE` se použijí v reziduu i v
příští hybnostní rovnici.

### 9–10. Korekce rychlosti a tlaku

```cpp
U[c]  = HbyA[c] - rAU[c]*gpn[c];          // dosadí vyřešený tlak
p[c]  = pOld[c] + alphaP*(p[c]-pOld[c]);  // podrelaxace tlaku
```

### 11. Reziduum a diagnostika

```cpp
resField = divFlux(m, phiE,phiN,phiF, U_in);   // div(U) po korekci
cont = Σ|resField[c]|;                          // L1 norma
// + ΔU/ΔP, zápis do convergence.dat
if(cont/contInit < tol) break;                  // konvergence?
```

## Proč právě toto pořadí

SIMPLE je **prediktor–korektor**: nejdřív odhad rychlosti (1–3), pak tlak, který
ten odhad opraví do kontinuity (4–9). Pořadí je nutné — tlaková rovnice
potřebuje `HbyA` z hybnostní rovnice, korekce rychlosti potřebuje vyřešený tlak.
Iterace se opakuje, protože hybnost a tlak jsou provázané (nelineární vazba přes
konvekci).
