#pragma once
#include <vector>

// Uniformní 3D kolokovaná mřížka s maskou pevné překážky.
struct Mesh {
    int nx, ny, nz;
    double dx, dy, dz, V;
    double Ax, Ay, Az;           // plochy stěn (E/W, N/S, F/B)
    std::vector<bool> solid;

    Mesh(int nx_, int ny_, int nz_, double Lx, double Ly, double Lz)
        : nx(nx_), ny(ny_), nz(nz_),
          dx(Lx/nx_), dy(Ly/ny_), dz(Lz/nz_),
          V(dx*dy*dz),
          Ax(dy*dz), Ay(dx*dz), Az(dx*dy),
          solid(nx_*ny_*nz_, false) {}

    int N() const { return nx*ny*nz; }
    int id(int i, int j, int k) const { return i + j*nx + k*nx*ny; }

    void markSolid(double x0, double x1,
                   double y0, double y1,
                   double z0, double z1) {
        auto rnd=[](double v){ return (int)(v+0.5); };
        int i0=rnd(x0/dx), i1=rnd(x1/dx);
        int j0=rnd(y0/dy), j1=rnd(y1/dy);
        int k0=rnd(z0/dz), k1=rnd(z1/dz);
        for (int k=k0; k<k1; ++k)
            for (int j=j0; j<j1; ++j)
                for (int i=i0; i<i1; ++i)
                    solid[id(i,j,k)] = true;
    }
};
