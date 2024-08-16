#pragma once
#include "auxsolver_base.hpp"

// Specify postprocessors that depend on one or more gridfunctions
namespace hephaestus
{

double calcMaxwellStressTensor(mfem::GridFunction * b_field, mfem::GridFunction * h_field, int face_attr);

double calcMaxwellStressTensor(mfem::GridFunction * b_field, mfem::GridFunction * h_field, int face_attr, mfem::Coefficient & q);

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

  std::string _b_name;  // name of the vector variable
  std::string _h_name;  // name of the vector variable
  std::string _coef_name; // name of the coefficient

  mfem::Array<double> _times;
  mfem::Array<double> _forces;

  mfem::ParGridFunction * _b_gf{nullptr};
  mfem::ParGridFunction * _h_gf{nullptr};
  mfem::Coefficient * _coef{nullptr};
  int _face_attr;
};

} // namespace hephaestus
