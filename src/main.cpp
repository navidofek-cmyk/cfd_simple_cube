#include "Simple.hpp"
#include "io.hpp"
#include <cstdio>

int main()
{
    // Řádkové bufferování — progress je vidět živě i při přesměrování do souboru
    setvbuf(stdout, nullptr, _IOLBF, 0);

    // Hranatý tunel: 4×1×1 m, čtvercový průřez
    const int    NX=320, NY=80, NZ=80;
    const double LX=4.0, LY=1.0, LZ=1.0;

    // Krychle 0.2×0.2×0.2 m uprostřed tunelu
    const double CX0=0.9, CX1=1.1;
    const double CY0=0.4, CY1=0.6;
    const double CZ0=0.4, CZ1=0.6;
    const double D = CX1-CX0;

    // Re = U_in * D / nu = 20  →  laminární ustálené proudění
    const double RE=20.0, U_IN=1.0;

    std::printf("=== 3D Channel flow around cube ===\n");
    std::printf("Mesh : %d×%d×%d  (%.0fk buněk)  dx=dy=dz=%.4f\n",
                NX,NY,NZ,(double)NX*NY*NZ/1000,LX/NX);
    std::printf("Cube : x=[%.1f,%.1f] y=[%.1f,%.1f] z=[%.1f,%.1f]  D=%.2f\n",
                CX0,CX1,CY0,CY1,CZ0,CZ1,D);
    std::printf("Re=%.0f  U_in=%.1f  nu=%.4f\n\n",RE,U_IN,U_IN*D/RE);

    // Diagnostický běh: tol=1e-9 (neukončí předčasně), maxIter=410,
    // snapshoty polí v iter 250/300/400 pro porovnání.
    ChannelSolver sim(NX,NY,NZ, LX,LY,LZ,
                      CX0,CX1, CY0,CY1, CZ0,CZ1,
                      RE, U_IN, 0.7, 0.3, 410, 1e-9);
    sim.snapIters = {250, 300, 400};
    sim.convergencePath = "convergence.dat";   // živý zápis reziduí
    // Export matic jen pro malé sítě (velké by byly stovky MB)
    if(NX*NY*NZ <= 100000){
        sim.matrixExportPath  = "pmatrix.dat";
        sim.matrixExportPathU = "umatrix.dat";
    }
    sim.solve();

    io::writeVTK   ("channel_cube.vtk", sim.m, sim.U, sim.p, sim.resField);
    io::writeFields("fields.dat",       sim.m, sim.U, sim.p, sim.resField);
    io::writeMesh  ("mesh.dat",         sim.m);

    std::printf("\nVýstup: channel_cube.vtk  fields.dat  mesh.dat\n");
    return 0;
}
