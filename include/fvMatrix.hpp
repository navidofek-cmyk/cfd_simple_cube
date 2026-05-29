#pragma once
#include <vector>
#include <cmath>
#include <omp.h>
#include "Mesh.hpp"
#include "Field.hpp"

template<class T>
struct fvMatrix {
    const Mesh& m;
    Field<T>& psi;
    std::vector<double> aP, aE, aW, aN, aS, aF, aB;
    std::vector<T>      src;

    fvMatrix(const Mesh& m_, Field<T>& psi_)
        : m(m_), psi(psi_),
          aP(m_.N(),0.0), aE(m_.N(),0.0), aW(m_.N(),0.0),
          aN(m_.N(),0.0), aS(m_.N(),0.0),
          aF(m_.N(),0.0), aB(m_.N(),0.0),
          src(m_.N(), T{}) {}

    void relax(double alpha){
        #pragma omp parallel for
        for(int c=0;c<m.N();++c){
            if(m.solid[c]) continue;
            double aOld = aP[c];
            aP[c] /= alpha;
            src[c] += (aP[c]-aOld) * psi[c];
        }
    }

    volScalar A() const {
        volScalar a(m.N(), 1.0);
        #pragma omp parallel for
        for(int c=0;c<m.N();++c) a[c] = aP[c]/m.V;
        return a;
    }

    Field<T> H() const {
        Field<T> h(m.N());
        #pragma omp parallel for collapse(3)
        for(int k=0;k<m.nz;++k)
        for(int j=0;j<m.ny;++j)
        for(int i=0;i<m.nx;++i){
            int c=m.id(i,j,k); T sum=src[c];
            if(i<m.nx-1) sum=sum-(aE[c]*psi[m.id(i+1,j,k)]);
            if(i>0)      sum=sum-(aW[c]*psi[m.id(i-1,j,k)]);
            if(j<m.ny-1) sum=sum-(aN[c]*psi[m.id(i,j+1,k)]);
            if(j>0)      sum=sum-(aS[c]*psi[m.id(i,j-1,k)]);
            if(k<m.nz-1) sum=sum-(aF[c]*psi[m.id(i,j,k+1)]);
            if(k>0)      sum=sum-(aB[c]*psi[m.id(i,j,k-1)]);
            h[c] = (1.0/m.V)*sum;
        }
        return h;
    }

    // Jacobi iterace (paralelní) místo Gauss-Seidel.
    // Každé vlákno čte STARÉ hodnoty psi, zapisuje do tmp — žádný race condition.
    void solveWithExtra(const Field<T>& extra, double scale, int sweeps) {
        std::vector<T> tmp(m.N());
        double res0 = 1.0;
        for(int s=0;s<sweeps;++s){
            double res=0.0;
            #pragma omp parallel for collapse(3) reduction(+:res)
            for(int k=0;k<m.nz;++k)
            for(int j=0;j<m.ny;++j)
            for(int i=0;i<m.nx;++i){
                int c=m.id(i,j,k);
                T rhs=src[c];
                if(!m.solid[c]) rhs=rhs+(scale*extra[c]);
                if(i<m.nx-1) rhs=rhs-(aE[c]*psi[m.id(i+1,j,k)]);
                if(i>0)      rhs=rhs-(aW[c]*psi[m.id(i-1,j,k)]);
                if(j<m.ny-1) rhs=rhs-(aN[c]*psi[m.id(i,j+1,k)]);
                if(j>0)      rhs=rhs-(aS[c]*psi[m.id(i,j-1,k)]);
                if(k<m.nz-1) rhs=rhs-(aF[c]*psi[m.id(i,j,k+1)]);
                if(k>0)      rhs=rhs-(aB[c]*psi[m.id(i,j,k-1)]);
                T np=(1.0/aP[c])*rhs;
                res+=mag2(np-psi[c]);
                tmp[c]=np;
            }
            std::swap(psi.d, tmp);
            if(s==0) res0=std::sqrt(res)+1e-30;
            if(std::sqrt(res)/res0<1e-3) break;
        }
    }
};
