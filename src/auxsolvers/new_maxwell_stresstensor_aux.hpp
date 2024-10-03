#pragma once
#include "auxsolver_base.hpp"

// Specify postprocessors that depend on one or more gridfunctions
namespace hephaestus
{
double calcMaxwellStressTensor(mfem::GridFunction * b_field, mfem::GridFunction * h_field, int face_attr, mfem::Coefficient & q);

double calcSurfaceForceDensity(mfem::GridFunction * b_field, mfem::GridFunction * h_field, int face_attr);

double calcSurfaceForceDensity(mfem::GridFunction * b_field, mfem::GridFunction * h_field, int face_attr, mfem::Coefficient & q);

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

  // Initialises the child submesh.
  void InitChildMesh();

  // Creates the relevant FE Collections and Spaces for the child submesh.
  void MakeFESpaces(int stage);

  // Creates the relevant GridFunctions for the child submesh.
  void MakeGridFunctions(int stage);

  std::string _b_name;  // name of the vector variable
  std::string _h_name;  // name of the vector variable
  std::string _coef_name; // name of the coefficient

  mfem::Array<double> _times;
  mfem::Array<double> _forces;

  mfem::ParGridFunction * _b_gf{nullptr};
  mfem::ParGridFunction * _h_gf{nullptr};
  // mfem::Coefficient * _coef{nullptr};
  mfem::ParGridFunction * _gf{nullptr};

  std::shared_ptr<mfem::ParGridFunction> _gf_child{nullptr};
  std::shared_ptr<mfem::ParGridFunction> _b_gf_child{nullptr};
  std::shared_ptr<mfem::ParGridFunction> _h_gf_child{nullptr};

  mfem::ParMesh * _mesh_parent{nullptr};
  std::unique_ptr<mfem::ParSubMesh> _mesh_child{nullptr};

  std::shared_ptr<mfem::ParFiniteElementSpace> _h1_fe_space_child{nullptr};
  std::unique_ptr<mfem::H1_FECollection> _h1_fe_space_fec_child{nullptr};

  std::shared_ptr<mfem::ParFiniteElementSpace> _h_div_fe_space_child{nullptr};
  std::unique_ptr<mfem::RT_FECollection> _h_div_fe_space_fec_child{nullptr};

  std::shared_ptr<mfem::ParFiniteElementSpace> _h_curl_fe_space_child{nullptr};
  std::unique_ptr<mfem::ND_FECollection> _h_curl_fe_space_fec_child{nullptr};

  std::shared_ptr<mfem::ParFiniteElementSpace> _l2_fe_space_child{nullptr};
  std::unique_ptr<mfem::L2_FECollection> _l2_fe_space_fec_child{nullptr};
  
  int _face_attr;
private:
  int _order_h1;
  int _order_hdiv;
  int _order_hcurl;
};

} // namespace hephaestus
