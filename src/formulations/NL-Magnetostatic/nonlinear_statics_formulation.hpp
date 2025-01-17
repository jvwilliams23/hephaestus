#pragma once
#include "../common/pfem_extras.hpp"
#include "formulation.hpp"
#include "inputs.hpp"
#include "sources.hpp"

namespace hephaestus
{

class NLStaticsFormulation : public SteadyStateEMFormulation
{
public:
  NLStaticsFormulation(std::string alpha_coef_name, std::string h_curl_var_name);

  ~NLStaticsFormulation() override = default;

  void ConstructJacobianPreconditioner() override;

  void ConstructJacobianSolver() override;

  void ConstructOperator() override;

  void RegisterGridFunctions() override;

  void RegisterCoefficients() override;

protected:
  const std::string _alpha_coef_name;
  const std::string _h_curl_var_name;
};

class NLStaticsOperator : public ProblemOperator
{
public:
  NLStaticsOperator(hephaestus::Problem & problem,
                    std::string h_curl_var_name,
                    std::string stiffness_coef_name);

  ~NLStaticsOperator() override = default;

  void SetGridFunctions() override;
  void Init(mfem::Vector & X) override;
  void Solve(mfem::Vector & X) override;

private:
  std::string _h_curl_var_name, _stiffness_coef_name;

  mfem::Coefficient * _stiff_coef{nullptr}; // Stiffness Material Coefficient
};

} // namespace hephaestus
