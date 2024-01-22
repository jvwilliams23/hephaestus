#include "vector_coefficient_aux.hpp"

#include <utility>

namespace hephaestus
{

VectorCoefficientAux::VectorCoefficientAux(std::string gf_name, std::string vec_coef_name)
  : _gf_name(std::move(gf_name)), _vec_coef_name(std::move(vec_coef_name))
{
}

void
VectorCoefficientAux::Init(const hephaestus::GridFunctions & gridfunctions,
                           hephaestus::Coefficients & coefficients)
{
  _gf = gridfunctions.Get(_gf_name);
  if (_gf == nullptr)
  {
    MFEM_ABORT("GridFunction " << _gf_name << " not found when initializing VectorCoefficientAux");
  }
  _vec_coef = coefficients._vectors.Get(_vec_coef_name);
  if (_vec_coef == nullptr)
  {
    MFEM_ABORT("VectorCoefficient " << _vec_coef_name
                                    << " not found when initializing VectorCoefficientAux");
  }
}

void
VectorCoefficientAux::Solve(double t)
{
  _vec_coef->SetTime(t);
  _gf->ProjectCoefficient(*_vec_coef);
}

} // namespace hephaestus
