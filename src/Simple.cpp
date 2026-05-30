#include "Simple.hpp"
#include "fvMatrix.hpp"
#include "CgSolver.hpp"
#include "io.hpp"
#include <cmath>
#include <cstdio>

// ─────────────────────────────────────────────
// Konstruktor
// ─────────────────────────────────────────────
ChannelSolver::ChannelSolver(int nx, int ny, int nz,
                             double Lx, double Ly, double Lz,
                             double cubeX0, double cubeX1,
                             double cubeY0, double cubeY1,
                             double cubeZ0, double cubeZ1,
                             double Re_, double U_in_,
                             double aU, double aP_,
                             int maxIter_, double tol_)
    : m(nx, ny, nz, Lx, Ly, Lz),
      Re(Re_), U_in(U_in_),
      nu(U_in_ * (cubeX1-cubeX0) / Re_),
      alphaU(aU), alphaP(aP_),
      maxIter(maxIter_), tol(tol_),
      U(nx*ny*nz,  Vec3{U_in_,0.0,0.0}),
      p(nx*ny*nz,  0.0),
      phiE(nx*ny*nz, 0.0), phiN(nx*ny*nz, 0.0), phiF(nx*ny*nz, 0.0),
      resField(nx*ny*nz, 0.0)
{
    m.markSolid(cubeX0, cubeX1, cubeY0, cubeY1, cubeZ0, cubeZ1);

    for (int c=0;c<m.N();++c) if (m.solid[c]) U[c] = Vec3{0,0,0};

    // Inicializace toků z počátečního pole
    for(int k=0;k<m.nz;++k)
    for(int j=0;j<m.ny;++j)
    for(int i=0;i<m.nx;++i){
        int c = m.id(i,j,k);
        phiE[c] = (i<m.nx-1 && !m.solid[c] && !m.solid[m.id(i+1,j,k)])
                  ? 0.5*(U[c].x+U[m.id(i+1,j,k)].x)*m.Ax : 0.0;
        phiN[c] = (j<m.ny-1 && !m.solid[c] && !m.solid[m.id(i,j+1,k)])
                  ? 0.5*(U[c].y+U[m.id(i,j+1,k)].y)*m.Ay : 0.0;
        phiF[c] = (k<m.nz-1 && !m.solid[c] && !m.solid[m.id(i,j,k+1)])
                  ? 0.5*(U[c].z+U[m.id(i,j,k+1)].z)*m.Az : 0.0;
    }
    // Outlet tok
    for(int k=0;k<m.nz;++k) for(int j=0;j<m.ny;++j){
        int c=m.id(m.nx-1,j,k);
        phiE[c] = m.solid[c] ? 0.0 : U[c].x*m.Ax;
    }

    cgWork.ensure(m.N());
}

// ─────────────────────────────────────────────
// Hybnostní rovnice (div − nu·laplacian)
// ─────────────────────────────────────────────
static fvMatrix<Vec3> buildUEqn(const Mesh& m,
                                 const volScalar& phiE,
                                 const volScalar& phiN,
                                 const volScalar& phiF,
                                 double nu, double U_in,
                                 volVector& U)
{
    fvMatrix<Vec3> eq(m, U);
    const double gx  = nu*m.Ax/m.dx,  gy  = nu*m.Ay/m.dy,  gz  = nu*m.Az/m.dz;
    const double gbx = nu*m.Ax/(0.5*m.dx);
    const double gby = nu*m.Ay/(0.5*m.dy);
    const double gbz = nu*m.Az/(0.5*m.dz);
    const Vec3 Uin{U_in, 0.0, 0.0};

    // Buňkový přístup: každá buňka c počítá VŠECHNY své příspěvky sama
    // (difúze + konvekce), jen ČTE sousedy a toky. Žádné zápisy do sousedů
    // → žádný race condition → lze paralelizovat jedním pragma.
    //
    // Konvenice toků: phiE[c]=odtok přes východní stěnu c (kladný = ven na východ).
    // Pro západní stěnu c je tok = phiE[w] (odtok z w na východ = vtok do c),
    // takže z pohledu c je odtok na západ = -phiE[w].
    #pragma omp parallel for collapse(3)
    for(int k=0;k<m.nz;++k)
    for(int j=0;j<m.ny;++j)
    for(int i=0;i<m.nx;++i){
        int c = m.id(i,j,k);
        if (m.solid[c]){ eq.aP[c]=1.0; continue; }

        // ── Západ ──────────────────────────────────────
        if (i>0){ int w=m.id(i-1,j,k);
            if(!m.solid[w]){
                eq.aW[c]-=gx; eq.aP[c]+=gx;             // difúze
                double f=-phiE[w];                       // odtok z c na západ
                if(f>=0) eq.aP[c]+=f; else eq.aW[c]+=f; // konvekce upwind
            } else { eq.aP[c]+=gbx; }
        } else {
            eq.aP[c]  += gbx;
            eq.src[c]  = eq.src[c] + (gbx*Uin);          // difúze inlet
            eq.src[c]  = eq.src[c] + (U_in*m.Ax*Uin);   // konvekce inlet
        }

        // ── Východ ─────────────────────────────────────
        if (i<m.nx-1){ int e=m.id(i+1,j,k);
            if(!m.solid[e]){
                eq.aE[c]-=gx; eq.aP[c]+=gx;
                double f=phiE[c];                        // odtok z c na východ
                if(f>=0) eq.aP[c]+=f; else eq.aE[c]+=f;
            } else { eq.aP[c]+=gbx; }
        } else {
            // outlet (Neumann difúze) + konvekce
            double f=phiE[c];
            if(f>=0) eq.aP[c]+=f;
            else     eq.src[c]=eq.src[c]+(-f*Uin);
        }

        // ── Jih ────────────────────────────────────────
        if (j>0){ int s=m.id(i,j-1,k);
            if(!m.solid[s]){
                eq.aS[c]-=gy; eq.aP[c]+=gy;
                double f=-phiN[s];
                if(f>=0) eq.aP[c]+=f; else eq.aS[c]+=f;
            } else { eq.aP[c]+=gby; }
        } else { eq.aP[c]+=gby; }

        // ── Sever ──────────────────────────────────────
        if (j<m.ny-1){ int n=m.id(i,j+1,k);
            if(!m.solid[n]){
                eq.aN[c]-=gy; eq.aP[c]+=gy;
                double f=phiN[c];
                if(f>=0) eq.aP[c]+=f; else eq.aN[c]+=f;
            } else { eq.aP[c]+=gby; }
        } else { eq.aP[c]+=gby; }

        // ── Zadní (z-) ─────────────────────────────────
        if (k>0){ int b=m.id(i,j,k-1);
            if(!m.solid[b]){
                eq.aB[c]-=gz; eq.aP[c]+=gz;
                double f=-phiF[b];
                if(f>=0) eq.aP[c]+=f; else eq.aB[c]+=f;
            } else { eq.aP[c]+=gbz; }
        } else { eq.aP[c]+=gbz; }

        // ── Přední (z+) ────────────────────────────────
        if (k<m.nz-1){ int ff=m.id(i,j,k+1);
            if(!m.solid[ff]){
                eq.aF[c]-=gz; eq.aP[c]+=gz;
                double f=phiF[c];
                if(f>=0) eq.aP[c]+=f; else eq.aF[c]+=f;
            } else { eq.aP[c]+=gbz; }
        } else { eq.aP[c]+=gbz; }
    }

    return eq;
}

// ─────────────────────────────────────────────
// Tlaková rovnice laplacian(rAU, p) = div(HbyA flux)
// ─────────────────────────────────────────────
static fvMatrix<double> buildPEqn(const Mesh& m,
                                   const volScalar& rAU,
                                   volScalar& p)
{
    fvMatrix<double> L(m,p);
    auto rf=[&](int a,int b){return 0.5*(rAU[a]+rAU[b]);};

    #pragma omp parallel for collapse(3)
    for(int k=0;k<m.nz;++k)
    for(int j=0;j<m.ny;++j)
    for(int i=0;i<m.nx;++i){
        int c=m.id(i,j,k);
        if(m.solid[c]){L.aP[c]=1.0; continue;}

        if(i>0){  int w=m.id(i-1,j,k);
            if(!m.solid[w]){double g=rf(c,w)*m.Ax/m.dx; L.aW[c]+=g; L.aP[c]-=g;}}
        // i==0: Neumann

        if(i<m.nx-1){int e=m.id(i+1,j,k);
            if(!m.solid[e]){double g=rf(c,e)*m.Ax/m.dx; L.aE[c]+=g; L.aP[c]-=g;}
        } else {
            double g=rAU[c]*m.Ax/(0.5*m.dx);
            L.aP[c]-=g;  // outlet Dirichlet p=0
        }

        if(j>0){  int s=m.id(i,j-1,k);
            if(!m.solid[s]){double g=rf(c,s)*m.Ay/m.dy; L.aS[c]+=g; L.aP[c]-=g;}}
        if(j<m.ny-1){int n=m.id(i,j+1,k);
            if(!m.solid[n]){double g=rf(c,n)*m.Ay/m.dy; L.aN[c]+=g; L.aP[c]-=g;}}

        if(k>0){  int b=m.id(i,j,k-1);
            if(!m.solid[b]){double g=rf(c,b)*m.Az/m.dz; L.aB[c]+=g; L.aP[c]-=g;}}
        if(k<m.nz-1){int f=m.id(i,j,k+1);
            if(!m.solid[f]){double g=rf(c,f)*m.Az/m.dz; L.aF[c]+=g; L.aP[c]-=g;}}
    }
    return L;
}

// ─────────────────────────────────────────────
// Tok z HbyA (bez tlakové korekce)
// ─────────────────────────────────────────────
static void interpFlux(const Mesh& m, const volVector& H,
                        volScalar& pE, volScalar& pN, volScalar& pF)
{
    #pragma omp parallel for collapse(3)
    for(int k=0;k<m.nz;++k)
    for(int j=0;j<m.ny;++j)
    for(int i=0;i<m.nx;++i){
        int c=m.id(i,j,k);
        if(i<m.nx-1){
            int e=m.id(i+1,j,k);
            pE[c]=(!m.solid[c]&&!m.solid[e])?0.5*(H[c].x+H[e].x)*m.Ax:0.0;
        } else {
            pE[c]=m.solid[c]?0.0:H[c].x*m.Ax;
        }
        if(j<m.ny-1){
            int n=m.id(i,j+1,k);
            pN[c]=(!m.solid[c]&&!m.solid[n])?0.5*(H[c].y+H[n].y)*m.Ay:0.0;
        } else { pN[c]=0.0; }
        if(k<m.nz-1){
            int f=m.id(i,j,k+1);
            pF[c]=(!m.solid[c]&&!m.solid[f])?0.5*(H[c].z+H[f].z)*m.Az:0.0;
        } else { pF[c]=0.0; }
    }
}

// ─────────────────────────────────────────────
// Divergence toku (zahrnuje inlet/outlet)
// ─────────────────────────────────────────────
static volScalar divFlux(const Mesh& m,
                          const volScalar& phiE, const volScalar& phiN,
                          const volScalar& phiF, double U_in)
{
    volScalar d(m.N(),0.0);
    #pragma omp parallel for collapse(3)
    for(int k=0;k<m.nz;++k)
    for(int j=0;j<m.ny;++j)
    for(int i=0;i<m.nx;++i){
        int c=m.id(i,j,k);
        if(m.solid[c]) continue;
        double s=phiE[c];
        s -= (i>0) ? phiE[m.id(i-1,j,k)] : U_in*m.Ax;   // inlet
        s += phiN[c];
        s -= (j>0) ? phiN[m.id(i,j-1,k)] : 0.0;
        s += phiF[c];
        s -= (k>0) ? phiF[m.id(i,j,k-1)] : 0.0;
        d[c]=s;
    }
    return d;
}

// ─────────────────────────────────────────────
// Gradient tlaku (kanálové BC)
// ─────────────────────────────────────────────
static volVector gradP(const Mesh& m, const volScalar& p)
{
    volVector g(m.N());
    #pragma omp parallel for collapse(3)
    for(int k=0;k<m.nz;++k)
    for(int j=0;j<m.ny;++j)
    for(int i=0;i<m.nx;++i){
        int c=m.id(i,j,k);
        if(m.solid[c]) continue;

        double pe=(i<m.nx-1&&!m.solid[m.id(i+1,j,k)])?0.5*(p[c]+p[m.id(i+1,j,k)])
                  :(i==m.nx-1?0.0:p[c]);
        double pw=(i>0&&!m.solid[m.id(i-1,j,k)])?0.5*(p[c]+p[m.id(i-1,j,k)]):p[c];
        double pn=(j<m.ny-1&&!m.solid[m.id(i,j+1,k)])?0.5*(p[c]+p[m.id(i,j+1,k)]):p[c];
        double ps=(j>0&&!m.solid[m.id(i,j-1,k)])?0.5*(p[c]+p[m.id(i,j-1,k)]):p[c];
        double pf=(k<m.nz-1&&!m.solid[m.id(i,j,k+1)])?0.5*(p[c]+p[m.id(i,j,k+1)]):p[c];
        double pb=(k>0&&!m.solid[m.id(i,j,k-1)])?0.5*(p[c]+p[m.id(i,j,k-1)]):p[c];

        g[c].x=(pe-pw)*m.Ax/m.V;
        g[c].y=(pn-ps)*m.Ay/m.V;
        g[c].z=(pf-pb)*m.Az/m.V;
    }
    return g;
}

// ─────────────────────────────────────────────
// Hlavní SIMPLE smyčka
// ─────────────────────────────────────────────
void ChannelSolver::solve()
{
    double contInit=0.0;
    volVector Uprev = U;   // pole rychlosti z minulé iterace (pro ΔU)

    // Živý zápis reziduí do souboru (každou iteraci)
    std::FILE* convFile = nullptr;
    if(!convergencePath.empty()){
        convFile = std::fopen(convergencePath.c_str(), "w");
        if(convFile){
            std::fprintf(convFile,
                "# Historie konvergence SIMPLE (3D kanal + krychle)\n");
            std::fprintf(convFile,
                "# iter  cont_rel_L1   dU/U_rel_L2   dP/P_rel_L2\n");
        }
    }

    for(int it=0;it<maxIter;++it){

        Uprev = U;

        // 1. Hybnostní rovnice → U*
        auto UEqn=buildUEqn(m,phiE,phiN,phiF,nu,U_in,U);
        UEqn.relax(alphaU);
        auto gp=gradP(m,p);
        UEqn.solveWithExtra(gp,-m.V,40);
        #pragma omp parallel for
        for(int c=0;c<m.N();++c) if(m.solid[c]) U[c]=Vec3{0,0,0};

        // 2. A, HbyA
        volScalar A=UEqn.A();
        volVector Hf=UEqn.H();
        volScalar rAU(m.N(),1.0);
        volVector HbyA(m.N());
        #pragma omp parallel for
        for(int c=0;c<m.N();++c){
            if(m.solid[c]){rAU[c]=0.0; continue;}
            rAU[c]=1.0/A[c];
            HbyA[c]=rAU[c]*Hf[c];
        }

        // 3. HbyA toky
        volScalar pE(m.N(),0.0),pN(m.N(),0.0),pF2(m.N(),0.0);
        interpFlux(m,HbyA,pE,pN,pF2);

        // 4. Tlaková rovnice
        volScalar pOld=p;
        auto pEqn=buildPEqn(m,rAU,p);
        auto divH=divFlux(m,pE,pN,pF2,U_in);
        #pragma omp parallel for
        for(int c=0;c<m.N();++c) pEqn.src[c]+=divH[c];
        cgSolve(pEqn,cgWork,1000,1e-6);

        // 5. Korekce toků (Rhie–Chow)
        auto rf=[&](int a,int b){return 0.5*(rAU[a]+rAU[b]);};
        #pragma omp parallel for collapse(3)
        for(int k=0;k<m.nz;++k)
        for(int j=0;j<m.ny;++j)
        for(int i=0;i<m.nx;++i){
            int c=m.id(i,j,k);
            if(m.solid[c]){phiE[c]=phiN[c]=phiF[c]=0.0; continue;}

            if(i<m.nx-1){int e=m.id(i+1,j,k);
                phiE[c]=!m.solid[e]?pE[c]-rf(c,e)*m.Ax/m.dx*(p[e]-p[c]):0.0;
            } else {
                phiE[c]=pE[c]-rAU[c]*m.Ax/(0.5*m.dx)*(0.0-p[c]);
            }
            if(j<m.ny-1){int n=m.id(i,j+1,k);
                phiN[c]=!m.solid[n]?pN[c]-rf(c,n)*m.Ay/m.dy*(p[n]-p[c]):0.0;
            } else { phiN[c]=0.0; }
            if(k<m.nz-1){int f=m.id(i,j,k+1);
                phiF[c]=!m.solid[f]?pF2[c]-rf(c,f)*m.Az/m.dz*(p[f]-p[c]):0.0;
            } else { phiF[c]=0.0; }
        }

        // 6. Korekce rychlosti
        auto gpn=gradP(m,p);
        #pragma omp parallel for
        for(int c=0;c<m.N();++c){
            if(m.solid[c]){U[c]=Vec3{0,0,0}; continue;}
            U[c]=HbyA[c]-rAU[c]*gpn[c];
        }

        // 7. Podrelaxace tlaku
        #pragma omp parallel for
        for(int c=0;c<m.N();++c)
            p[c]=pOld[c]+alphaP*(p[c]-pOld[c]);

        // 8. Kontinuitní reziduum + změna polí mezi iteracemi
        resField=divFlux(m,phiE,phiN,phiF,U_in);
        double cont=0.0, dU2=0.0, dP2=0.0, U2=0.0, P2=0.0;
        #pragma omp parallel for reduction(+:cont,dU2,dP2,U2,P2)
        for(int c=0;c<m.N();++c){
            cont += std::fabs(resField[c]);
            dU2  += mag2(U[c]-Uprev[c]);          // ΔU = U_it − U_{it-1}
            dP2  += (p[c]-pOld[c])*(p[c]-pOld[c]);// Δp v této iteraci
            U2   += mag2(U[c]);
            P2   += p[c]*p[c];
        }
        if(it==0) contInit=cont+1e-30;

        // relativní L2 změny polí (bezrozměrné)
        double relU = std::sqrt(dU2)/(std::sqrt(U2)+1e-30);
        double relP = std::sqrt(dP2)/(std::sqrt(P2)+1e-30);

        // Živý zápis reziduí (KAŽDOU iteraci, ne jen každou 20.)
        if(convFile){
            std::fprintf(convFile, "%6d  %.6e  %.6e  %.6e\n",
                         it, cont/contInit, relU, relP);
            std::fflush(convFile);
        }

        if(it%20==0){
            std::printf("iter %5d   cont=%.3e   dU/U=%.3e   dP/P=%.3e\n",
                        it, cont/contInit, relU, relP);
            std::fflush(stdout);
        }

        // Snapshot polí v zadaných iteracích
        for(int si : snapIters) if(it==si){
            char path[64];
            std::snprintf(path,sizeof(path),"fields_%04d.dat",it);
            io::writeFields(path, m, U, p, resField);
            std::printf("  >> snapshot ulozen: %s\n", path);
            std::fflush(stdout);
        }

        // Callback pro GUI (zive updaty + moznost predcasneho ukonceni)
        if(onIteration && !onIteration(it, cont/contInit)){
            std::printf("Preruseno uzivatelem v iteraci %d\n", it);
            std::fflush(stdout);
            break;
        }

        if(cont/contInit<tol){
            std::printf("Konvergováno v iteraci %d  (rel.res=%.2e)\n",it,cont/contInit);
            if(!matrixExportPath.empty())
                io::writePMatrix(matrixExportPath, m, pEqn);
            if(!matrixExportPathU.empty()){
                // Přebuduj UEqn bez relaxace pro čistou fyzikální matici
                auto UEqnFinal = buildUEqn(m,phiE,phiN,phiF,nu,U_in,U);
                io::writeUMatrix(matrixExportPathU, m, UEqnFinal);
            }
            break;
        }
    }
    if(convFile) std::fclose(convFile);
}
