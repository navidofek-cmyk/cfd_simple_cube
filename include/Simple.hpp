#pragma once
#include <string>
#include <vector>
#include "Mesh.hpp"
#include "Field.hpp"
#include "CgSolver.hpp"

// SIMPLE solver pro 3D kanálové proudění kolem krychle.
//
// Geometrie (tunel čtvercového průřezu):
//   Inlet  (x=0):       Dirichlet U=(U_in,0,0), Neumann dp/dn=0
//   Outlet (x=Lx):      Neumann dU/dn=0,        Dirichlet p=0
//   Stěny  (±y, ±z):    no-slip U=0,             Neumann dp/dn=0
//   Krychle (uvnitř):   no-slip U=0,             Neumann dp/dn=0
//
// Re = U_in * D / nu  (D = strana krychle)
struct ChannelSolver {
    Mesh   m;
    double Re, U_in, nu;
    double alphaU, alphaP;
    int    maxIter;
    double tol;

    volVector U;
    volScalar p;
    volScalar phiE, phiN, phiF;   // toky přes E/N/F stěny
    volScalar resField;

    CgWorkspace cgWork;
    std::string matrixExportPath;    // tlakova matice
    std::string matrixExportPathU;   // hybnostni matice
    std::string convergencePath;     // soubor pro zivy zapis rezidui (kazdou iteraci)
    std::vector<int> snapIters;      // iterace, ve kterych ulozit snapshot pole

    ChannelSolver(int nx, int ny, int nz,
                  double Lx, double Ly, double Lz,
                  double cubeX0, double cubeX1,
                  double cubeY0, double cubeY1,
                  double cubeZ0, double cubeZ1,
                  double Re_,    double U_in_ = 1.0,
                  double aU = 0.7, double aP = 0.3,
                  int maxIter_ = 6000, double tol_ = 1e-6);

    void solve();
};
