#pragma once
#include "auxsolver_base.hpp"
// #include "vector_boundarynormal_coefficient_aux.hpp"
// #include "boundary_coefficient_aux.hpp"

// Specify postprocessors that depend on one or more gridfunctions
namespace hephaestus
{

// The DevMaxwellStressTensorAuxCoefficient evaluates the surface force density as :
// F = [b_n^2 * (1/\mu_0 - 1/\mu) - h_t^2 * (\mu_0 - \mu)] * n * 1/2
// where b_n is the normal component of magnetic flux density
// h_t is the tangential component of the magnetic field at the interface
// \mu_0 and \mu are permeabilities of two materials
// n is the unit normal vector (of the surface)
// class DevMaxwellStressTensorAuxCoefficient : public mfem::Coefficient
class DevMaxwellStressTensorAuxCoefficient : public mfem::VectorCoefficient
{
private:
  // const mfem::ParGridFunction * _b_gf{nullptr};
  const std::shared_ptr<mfem::ParGridFunction> _b_gf{nullptr};
  const mfem::ParGridFunction * _h_gf{nullptr};

public:
  DevMaxwellStressTensorAuxCoefficient(
    // const mfem::ParGridFunction * b_gf,
    std::shared_ptr<mfem::ParGridFunction> b_gf,
    const mfem::ParGridFunction * h_gf
  )
    : mfem::VectorCoefficient(2), _b_gf{b_gf}, _h_gf{h_gf}
    // : _b_gf{b_gf}, _h_gf{h_gf}
  {
  }

  ~DevMaxwellStressTensorAuxCoefficient() override = default;

  void Eval(mfem::Vector & uxv,
            mfem::ElementTransformation & T,
            const mfem::IntegrationPoint & ip) override;
  // double Eval(mfem::ElementTransformation & T, const mfem::IntegrationPoint & ip) override;
};

// Auxsolver to project the dot product of two vector gridfunctions onto a third
// (scalar) GridFunction
class DevMaxwellStressTensorAux : public AuxSolver // : public VectorBoundaryNormalCoefficientAux
// class DevMaxwellStressTensorAux : public BoundaryCoefficientAux
{
private:
  mfem::ParGridFunction * _b_gf{nullptr};
  mfem::ParGridFunction * _h_gf{nullptr};

  const std::string _b_gf_name;
  const std::string _h_gf_name;

public:
  DevMaxwellStressTensorAux(const std::string & f_gf_name,
                            const std::string & f_coef_name,
                            std::string b_gf_name,
                            std::string h_gf_name,
                            mfem::Array<int> boundary_attr);

  // ~DevMaxwellStressTensorAux() override = default;

  void Init(const hephaestus::GridFunctions & gridfunctions,
            hephaestus::Coefficients & coefficients) override;

  virtual void BuildBilinearForm();
  virtual void BuildLinearForm();
  void Solve(double t = 0.0) override;

  // Initialises the child submesh.
  void InitChildMesh();

  // Creates the relevant FE Collections and Spaces for the child submesh.
  void MakeFESpaces(int stage);

  // Creates the relevant GridFunctions for the child submesh.
  void MakeGridFunctions(int stage);

protected:
  const std::string _gf_name;   // name of the variable
  const std::string _coef_name; // name of the coefficient

  mfem::ParMesh * _mesh_parent{nullptr};
  std::unique_ptr<mfem::ParSubMesh> _mesh_child{nullptr};
  std::shared_ptr<mfem::ParFiniteElementSpace> _h1_fe_space_child{nullptr};
  std::shared_ptr<mfem::ParFiniteElementSpace> _h_div_fe_space_child{nullptr};
  std::unique_ptr<mfem::H1_FECollection> _h1_fe_space_fec_child{nullptr};
  std::unique_ptr<mfem::RT_FECollection> _h_div_fe_space_fec_child{nullptr};
  
  mfem::ParGridFunction * _gf{nullptr};
  std::shared_ptr<mfem::ParGridFunction> _gf_child{nullptr};
  std::shared_ptr<mfem::ParGridFunction> _b_gf_child{nullptr};
  // mfem::ParGridFunction * _b_gf_child{nullptr};

  mfem::VectorCoefficient * _vec_coef{nullptr};
  mfem::Coefficient * _scalar_coef{nullptr};
  // mfem::Coefficient * _mass_coef{nullptr};
  std::shared_ptr<mfem::Coefficient> _mass_coef{nullptr};
  std::shared_ptr<mfem::Coefficient> _rt_boundary_coef{nullptr};
  mfem::Array<int> _boundary_attr; // int of attribute to limit boundary integration to
  mfem::Array<int> _boundary_attr_marker; // int of attribute to limit boundary integration to


  // Pointer to store test FE space. Assumed to be same as trial FE space.
  mfem::ParFiniteElementSpace * _test_fes{nullptr};
  mfem::ParFiniteElementSpace * _trial_fes{nullptr};

  // Bilinear and linear forms
  std::unique_ptr<mfem::ParBilinearForm> _a{nullptr};
  // std::unique_ptr<mfem::ParMixedBilinearForm> _a{nullptr};
  std::unique_ptr<mfem::ParLinearForm> _b{nullptr};

private:
  int _order_h1;
  int _order_hdiv;
  const hephaestus::InputParameters _solver_options;

  // Operator matrices
  std::unique_ptr<mfem::HypreParMatrix> _a_mat{nullptr};

  // Solver
  std::unique_ptr<hephaestus::DefaultJacobiPCGSolver> _solver{nullptr};
};

} // namespace hephaestus
