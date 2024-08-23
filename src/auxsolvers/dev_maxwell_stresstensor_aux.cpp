#include "dev_maxwell_stresstensor_aux.hpp"

#include <utility>

namespace hephaestus
{


void
DevMaxwellStressTensorAuxCoefficient::Eval(mfem::Vector & uxv,
                                                mfem::ElementTransformation & T,
                                                const mfem::IntegrationPoint & ip)
{
  double air_permeability = 1.25663706e-6;
  double sphere_permeability = 500*air_permeability;
  std::cout << "in Eval devMaxwell" << std::endl;

  mfem::Vector _b_gf_val;
  mfem::Vector _h_gf_val;
  _b_gf->GetVectorValue(T, ip, _b_gf_val);
  std::cout << "GetVectorValue " << _b_gf_val(0) << _b_gf_val(1) << std::endl;
}

DevMaxwellStressTensorAux::DevMaxwellStressTensorAux(
    const std::string & f_gf_name,
    const std::string & f_coef_name,
    std::string b_gf_name,
    std::string h_gf_name,
    mfem::Array<int> boundary_attr)
  : VectorBoundaryNormalCoefficientAux(f_gf_name, f_coef_name, boundary_attr),
    _b_gf_name(std::move(b_gf_name)),
    _h_gf_name(std::move(h_gf_name))
{
}

void
DevMaxwellStressTensorAux::Init(const hephaestus::GridFunctions & gridfunctions,
                                      hephaestus::Coefficients & coefficients)
{
  _b_gf = gridfunctions.Get(_b_gf_name);
  _h_gf = gridfunctions.Get(_h_gf_name);


  coefficients._vectors.Register(_vec_coef_name,
                                 std::make_shared<DevMaxwellStressTensorAuxCoefficient>(
                                     _b_gf, _h_gf));

  VectorBoundaryNormalCoefficientAux::Init(gridfunctions, coefficients);
}

} // namespace hephaestus
