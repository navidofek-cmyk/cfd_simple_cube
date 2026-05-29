#include "io.hpp"
#include <fstream>
#include <cmath>

namespace io {

void writeVTK(const std::string& path, const Mesh& m,
              const volVector& U, const volScalar& p,
              const volScalar& res)
{
    std::ofstream f(path);
    f << "# vtk DataFile Version 3.0\n3D channel cube flow\nASCII\n"
         "DATASET STRUCTURED_POINTS\n";
    f << "DIMENSIONS " << m.nx << " " << m.ny << " " << m.nz << "\n";
    f << "ORIGIN " << 0.5*m.dx << " " << 0.5*m.dy << " " << 0.5*m.dz << "\n";
    f << "SPACING " << m.dx << " " << m.dy << " " << m.dz << "\n";
    f << "POINT_DATA " << m.N() << "\n";

    f << "VECTORS velocity double\n";
    for(int c=0;c<m.N();++c) f<<U[c].x<<" "<<U[c].y<<" "<<U[c].z<<"\n";

    f << "SCALARS pressure double 1\nLOOKUP_TABLE default\n";
    for(int c=0;c<m.N();++c) f<<p[c]<<"\n";

    f << "SCALARS solid double 1\nLOOKUP_TABLE default\n";
    for(int c=0;c<m.N();++c) f<<(m.solid[c]?1.0:0.0)<<"\n";

    f << "SCALARS continuity_res double 1\nLOOKUP_TABLE default\n";
    for(int c=0;c<m.N();++c) f<<res[c]<<"\n";
}

void writeFields(const std::string& path, const Mesh& m,
                 const volVector& U, const volScalar& p,
                 const volScalar& res)
{
    std::ofstream f(path);
    f << "# DIMS " << m.nx << " " << m.ny << " " << m.nz << "\n";
    f << "# i j k x y z ux uy uz p res solid\n";
    for(int k=0;k<m.nz;++k)
    for(int j=0;j<m.ny;++j)
    for(int i=0;i<m.nx;++i){
        int c=m.id(i,j,k);
        f<<i<<" "<<j<<" "<<k<<" "
         <<(i+0.5)*m.dx<<" "<<(j+0.5)*m.dy<<" "<<(k+0.5)*m.dz<<" "
         <<U[c].x<<" "<<U[c].y<<" "<<U[c].z<<" "
         <<p[c]<<" "<<res[c]<<" "<<(m.solid[c]?1:0)<<"\n";
    }
}

void writeMesh(const std::string& path, const Mesh& m)
{
    std::ofstream f(path);
    f << "# DIMS " << m.nx << " " << m.ny << " " << m.nz
      << " DX " << m.dx << " DY " << m.dy << " DZ " << m.dz << "\n";
    f << "# i j k x y z solid\n";
    for(int k=0;k<m.nz;++k)
    for(int j=0;j<m.ny;++j)
    for(int i=0;i<m.nx;++i){
        int c=m.id(i,j,k);
        f<<i<<" "<<j<<" "<<k<<" "
         <<(i+0.5)*m.dx<<" "<<(j+0.5)*m.dy<<" "<<(k+0.5)*m.dz<<" "
         <<(m.solid[c]?1:0)<<"\n";
    }
}

void writePMatrix(const std::string& path, const Mesh& m,
                  const fvMatrix<double>& mat)
{
    std::ofstream f(path);
    f << "# DIMS " << m.nx << " " << m.ny << " " << m.nz << "\n";
    f << "# i j k x y z aP aE aW aN aS aF aB src solid\n";
    for(int k=0;k<m.nz;++k)
    for(int j=0;j<m.ny;++j)
    for(int i=0;i<m.nx;++i){
        int c=m.id(i,j,k);
        f<<i<<" "<<j<<" "<<k<<" "
         <<(i+0.5)*m.dx<<" "<<(j+0.5)*m.dy<<" "<<(k+0.5)*m.dz<<" "
         <<mat.aP[c]<<" "<<mat.aE[c]<<" "<<mat.aW[c]<<" "
         <<mat.aN[c]<<" "<<mat.aS[c]<<" "<<mat.aF[c]<<" "<<mat.aB[c]<<" "
         <<mat.src[c]<<" "<<(m.solid[c]?1:0)<<"\n";
    }
}

void writeUMatrix(const std::string& path, const Mesh& m,
                  const fvMatrix<Vec3>& mat)
{
    std::ofstream f(path);
    f << "# DIMS " << m.nx << " " << m.ny << " " << m.nz << "\n";
    f << "# i j k x y z aP aE aW aN aS aF aB src_x src_y src_z solid\n";
    for(int k=0;k<m.nz;++k)
    for(int j=0;j<m.ny;++j)
    for(int i=0;i<m.nx;++i){
        int c=m.id(i,j,k);
        f<<i<<" "<<j<<" "<<k<<" "
         <<(i+0.5)*m.dx<<" "<<(j+0.5)*m.dy<<" "<<(k+0.5)*m.dz<<" "
         <<mat.aP[c]<<" "<<mat.aE[c]<<" "<<mat.aW[c]<<" "
         <<mat.aN[c]<<" "<<mat.aS[c]<<" "<<mat.aF[c]<<" "<<mat.aB[c]<<" "
         <<mat.src[c].x<<" "<<mat.src[c].y<<" "<<mat.src[c].z<<" "
         <<(m.solid[c]?1:0)<<"\n";
    }
}

} // namespace io
