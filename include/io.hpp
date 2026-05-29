#pragma once
#include <string>
#include "Mesh.hpp"
#include "Field.hpp"
#include "Vec3.hpp"
#include "fvMatrix.hpp"

namespace io {

void writeVTK(const std::string& path, const Mesh& m,
              const volVector& U, const volScalar& p,
              const volScalar& res);

// Textový výpis pro Python (x y z ux uy uz p solid)
void writeFields(const std::string& path, const Mesh& m,
                 const volVector& U, const volScalar& p,
                 const volScalar& res);

// Výpis geometrie sítě
void writeMesh(const std::string& path, const Mesh& m);

// Export koeficientů tlakové matice
void writePMatrix(const std::string& path, const Mesh& m,
                  const fvMatrix<double>& mat);

// Export koeficientů hybnostní matice (Vec3 — zapisuje scalar koeficienty + src.xyz)
void writeUMatrix(const std::string& path, const Mesh& m,
                  const fvMatrix<Vec3>& mat);

} // namespace io
