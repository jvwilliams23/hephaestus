#include "scaled_gridfunction_aux.hpp"

namespace hephaestus {

ScaledGridFunctionAuxSolver::ScaledGridFunctionAuxSolver(
    const hephaestus::InputParameters &params)
    : AuxSolver(), coef(nullptr), input_gf(nullptr), scaled_gf(nullptr),
      coef_name(params.GetParam<std::string>("CoefficientName")),
      input_gf_name(params.GetParam<std::string>("InputVariableName")),
      scaled_gf_name(params.GetParam<std::string>("ScaledVariableName")) {}

void ScaledGridFunctionAuxSolver::Init(
    const hephaestus::GridFunctions &gridfunctions,
    hephaestus::Coefficients &coefficients) {
  input_gf = gridfunctions.Get(input_gf_name);
  if (input_gf == NULL) {
    MFEM_ABORT("GridFunction "
               << input_gf_name
               << " not found when initializing CoefficientAuxSolver");
  }
  scaled_gf = gridfunctions.Get(scaled_gf_name);
  if (scaled_gf == NULL) {
    MFEM_ABORT("GridFunction "
               << scaled_gf_name
               << " not found when initializing CoefficientAuxSolver");
  }
  coef = dynamic_cast<mfem::ConstantCoefficient *>(
      coefficients.scalars.Get(coef_name));
  if (coef == NULL) {
    MFEM_ABORT("ConstantCoefficient "
               << coef_name
               << " not found when initializing CoefficientAuxSolver");
  }
}

void ScaledGridFunctionAuxSolver::Solve(double t) {
  coef->SetTime(t);
  *scaled_gf = *input_gf;
  *scaled_gf *= coef->constant;
}

} // namespace hephaestus
