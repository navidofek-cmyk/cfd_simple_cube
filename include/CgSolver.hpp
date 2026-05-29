#pragma once
#include "fvMatrix.hpp"
#include <vector>
#include <cmath>
#include <omp.h>

struct CgWorkspace {
    std::vector<double> b, Ax, r, dir, Ap;
    void ensure(int N) {
        if(static_cast<int>(b.size())==N) return;
        b.resize(N); Ax.resize(N); r.resize(N); dir.resize(N); Ap.resize(N);
    }
};

inline void cgSolve(fvMatrix<double>& M, CgWorkspace& ws,
                    int maxIter=1000, double tol=1e-6) {
    const Mesh& m=M.m;
    const int   N=m.N();
    auto& x=M.psi.d;
    ws.ensure(N);

    // y = (-A)*v  — paralelní maticové násobení
    auto negAv=[&](const std::vector<double>& v, std::vector<double>& y){
        #pragma omp parallel for collapse(3)
        for(int k=0;k<m.nz;++k)
        for(int j=0;j<m.ny;++j)
        for(int i=0;i<m.nx;++i){
            int c=m.id(i,j,k);
            y[c]=-M.aP[c]*v[c];
            if(i<m.nx-1) y[c]-=M.aE[c]*v[m.id(i+1,j,k)];
            if(i>0)      y[c]-=M.aW[c]*v[m.id(i-1,j,k)];
            if(j<m.ny-1) y[c]-=M.aN[c]*v[m.id(i,j+1,k)];
            if(j>0)      y[c]-=M.aS[c]*v[m.id(i,j-1,k)];
            if(k<m.nz-1) y[c]-=M.aF[c]*v[m.id(i,j,k+1)];
            if(k>0)      y[c]-=M.aB[c]*v[m.id(i,j,k-1)];
        }
    };

    #pragma omp parallel for
    for(int c=0;c<N;++c) ws.b[c]=-M.src[c];

    std::fill(ws.Ax.begin(),ws.Ax.end(),0.0);
    negAv(x,ws.Ax);

    double rr=0.0;
    #pragma omp parallel for reduction(+:rr)
    for(int c=0;c<N;++c){ws.r[c]=ws.b[c]-ws.Ax[c]; rr+=ws.r[c]*ws.r[c];}
    ws.dir=ws.r;
    const double res0=std::sqrt(rr)+1e-30;

    for(int k=0;k<maxIter;++k){
        negAv(ws.dir,ws.Ap);

        double pAp=0.0;
        #pragma omp parallel for reduction(+:pAp)
        for(int c=0;c<N;++c) pAp+=ws.dir[c]*ws.Ap[c];
        double alpha=rr/(pAp+1e-30);

        #pragma omp parallel for
        for(int c=0;c<N;++c){x[c]+=alpha*ws.dir[c]; ws.r[c]-=alpha*ws.Ap[c];}

        double rr_new=0.0;
        #pragma omp parallel for reduction(+:rr_new)
        for(int c=0;c<N;++c) rr_new+=ws.r[c]*ws.r[c];
        if(std::sqrt(rr_new)/res0<tol) break;

        double beta=rr_new/(rr+1e-30);
        #pragma omp parallel for
        for(int c=0;c<N;++c) ws.dir[c]=ws.r[c]+beta*ws.dir[c];
        rr=rr_new;
    }
}
