#pragma once
#include "auxsolver_base.hpp"

// Specify postprocessors that depend on one or more gridfunctions
namespace hephaestus
{
double calcMaxwellStressTensor(mfem::ParGridFunction * b_field, mfem::ParGridFunction * h_field, int face_attr, mfem::Coefficient & q);

double calcSurfaceForceDensity(mfem::ParGridFunction * b_field, mfem::ParGridFunction * h_field, int face_attr, mfem::Coefficient & q, mfem::Coefficient & mu);

// Class to calculate and store the flux of a vector GridFunction through a surface
// at each timestep, optionally scaled by a coefficient.
class MaxwellStressTensorAux : public AuxSolver
{

public:
  MaxwellStressTensorAux() = default;
  MaxwellStressTensorAux(std::string b_name, std::string h_name, int face_attr, std::string coef_name = "");

  ~MaxwellStressTensorAux() override = default;

  void Init(const hephaestus::GridFunctions & gridfunctions,
            hephaestus::Coefficients & coefficients) override;

  void Solve(double t = 0.0) override;

  void WriteForces(std::string fname, mfem::ParGridFunction & gf, int face_attr);

  std::string _b_name;  // name of the vector variable
  std::string _h_name;  // name of the vector variable
  std::string _coef_name; // name of the coefficient

  mfem::Array<double> _times;
  mfem::Array<double> _forces;

  mfem::ParGridFunction * _b_gf{nullptr};
  mfem::ParGridFunction * _h_gf{nullptr};
  mfem::ParGridFunction * _gf{nullptr};
  mfem::Coefficient * _mu_coef{nullptr};

  mfem::ParMesh * _mesh_parent{nullptr};

  int _face_attr;
};

} // namespace hephaestus
