#pragma once
#include "auxsolver_base.hpp"

namespace hephaestus
{

// Scale a gridfunction in H(Curl) or H(Div) by a scalar Coefficient, and store
// the result. Suitable for solving for H(Div) or H(Curl) conforming fields for
// expressions like v = a*σ*u
class ScaledVectorGridFunctionAux : public AuxSolver
{
public:
  ScaledVectorGridFunctionAux(
      std::string input_gf_name,
      std::string scaled_gf_name,
      std::string coef_name,
      const double & aConst = 1.0,
      hephaestus::InputParameters solver_options = hephaestus::InputParameters());

  ~ScaledVectorGridFunctionAux() override = default;

  void Init(const hephaestus::GridFunctions & gridfunctions,
            hephaestus::Coefficients & coefficients) override;
  virtual void BuildBilinearForm();
  virtual void BuildMixedBilinearForm();
  void Solve(double t = 0.0) override;

protected:
  // Pointers to store trial and test FE spaces
  mfem::ParFiniteElementSpace * trial_fes{nullptr};
  mfem::ParFiniteElementSpace * test_fes{nullptr};

  // Bilinear forms
  std::unique_ptr<mfem::ParBilinearForm> a{nullptr};
  std::unique_ptr<mfem::ParMixedBilinearForm> a_mixed{nullptr};

  // Coefficient to scale input gridfunction by
  mfem::Coefficient * coef{nullptr};
  // Optional constant to scale input gridfunction by

private:
  const std::string _input_gf_name;
  const std::string _scaled_gf_name;
  const std::string _coef_name;
  const double _aConst;
  const hephaestus::InputParameters _solver_options;

  // Input gridfunction to be scaled by a scalar coefficient
  mfem::ParGridFunction * input_gf{nullptr};

  // Gridfunction in which to store result
  mfem::ParGridFunction * scaled_gf{nullptr};

  // Operator matrices
  std::unique_ptr<mfem::HypreParMatrix> a_mat{nullptr};
  std::unique_ptr<mfem::HypreParMatrix> mixed_mat{nullptr};

  // Solver
  std::unique_ptr<hephaestus::DefaultJacobiPCGSolver> solver{nullptr};
};
} // namespace hephaestus
