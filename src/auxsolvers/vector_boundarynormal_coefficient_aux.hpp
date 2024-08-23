#pragma once
#include "auxsolver_base.hpp"
#include "vector_coefficient_aux.hpp"

namespace hephaestus
{

  void inheritBdrAttributes(const mfem::ParMesh * parent_mesh, mfem::ParSubMesh * child_mesh);

// Interpolate a stored scalar Coefficient onto a (scalar) GridFunction
// by solving (p, q) = (λ, q)
class VectorBoundaryNormalCoefficientAux : public AuxSolver
{
public:
  VectorBoundaryNormalCoefficientAux(std::string gf_name,
                 std::string coef_name,
                 mfem::Array<int> boundary_attr,
                 hephaestus::InputParameters solver_options = hephaestus::InputParameters()
                 );

  ~VectorBoundaryNormalCoefficientAux() override = default;

  void Init(const hephaestus::GridFunctions & gridfunctions,
            hephaestus::Coefficients & coefficients) override;

  virtual void BuildBilinearForm();
  virtual void BuildLinearForm();
  // virtual void BuildLinearFormNormal();
  void Solve(double t = 0.0) override;

  // Initialises the child submesh.
  void InitChildMesh();

  // Creates the relevant FE Collections and Spaces for the child submesh.
  void MakeFESpaces();

  // Creates the relevant GridFunctions for the child submesh.
  void MakeGridFunctions();

protected:
  const std::string _gf_name;   // name of the variable
  const std::string _vec_coef_name; // name of the coefficient
  // std::vector<int> _attr_int; // int of attribute to limit boundary integration to


  mfem::ParMesh * _mesh_parent{nullptr};
  std::unique_ptr<mfem::ParSubMesh> _mesh_child{nullptr};
  std::shared_ptr<mfem::ParFiniteElementSpace> _h1_fe_space_child{nullptr};
  std::unique_ptr<mfem::H1_FECollection> _h1_fe_space_fec_child{nullptr};

  mfem::ParGridFunction * _gf{nullptr};
  std::shared_ptr<mfem::ParGridFunction> _gf_child{nullptr};
  mfem::VectorCoefficient * _vec_coef{nullptr};
  // mfem::Coefficient * _mass_coef{nullptr};
  std::shared_ptr<mfem::Coefficient> _mass_coef{nullptr};
  mfem::Array<int> _boundary_attr; // int of attribute to limit boundary integration to
  mfem::Array<int> _boundary_attr_marker; // int of attribute to limit boundary integration to


  // Pointer to store test FE space. Assumed to be same as trial FE space.
  mfem::ParFiniteElementSpace * _test_fes{nullptr};

  // Bilinear and linear forms
  std::unique_ptr<mfem::ParBilinearForm> _a{nullptr};
  std::unique_ptr<mfem::ParLinearForm> _b{nullptr};

private:
  int _order_h1;
  int _order_hcurl;
  int _order_hdiv;
  const hephaestus::InputParameters _solver_options;

  // Operator matrices
  std::unique_ptr<mfem::HypreParMatrix> _a_mat{nullptr};

  // Solver
  std::unique_ptr<hephaestus::DefaultJacobiPCGSolver> _solver{nullptr};
};

} // namespace hephaestus
