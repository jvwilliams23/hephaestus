#include "dev_maxwell_stresstensor_aux.hpp"

#include <utility>

namespace hephaestus
{

// double
// DevMaxwellStressTensorAuxCoefficient::Eval(mfem::ElementTransformation & T,
//                                         const mfem::IntegrationPoint & ip)
// {
//   double air_permeability = 1.25663706e-6;
//   double sphere_permeability = 500*air_permeability;
//   std::random_device rd;  // Will be used to obtain a seed for the random number engine
//   std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
//   std::uniform_real_distribution<> dist(0, 10000);
  
//   mfem::Vector _b_gf_val;
//   mfem::Vector _h_gf_val;
//   _b_gf->GetVectorValue(T, ip, _b_gf_val);
//   _h_gf->GetVectorValue(T, ip, _h_gf_val);

//   return _b_gf_val * _h_gf_val;
//   // return dist(gen);///1000.0;
//   // return 1000.0;
// }

double
DevMaxwellStressTensorAuxCoefficient::Eval(mfem::ElementTransformation & T,
                                        const mfem::IntegrationPoint & ip)
{
  double air_permeability = 1.25663706e-6;
  double sphere_permeability = 500*air_permeability;
  std::random_device rd;  // Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> dist(0, 10000);
  
  mfem::Vector _b_gf_val;
  mfem::Vector _h_gf_val;
  // _b_gf->GetVectorValue(T, ip, _b_gf_val);
  // _h_gf->GetVectorValue(T, ip, _h_gf_val);

  // return _b_gf_val * _h_gf_val;
  // return dist(gen);///1000.0;
  return 1000.0;
}

// double
// DevMaxwellStressTensorAuxCoefficient::Eval(mfem::FaceElementTransformations & T,
//                                         const mfem::IntegrationPoint & ip)
// {
//   double air_permeability = 1.25663706e-6;
//   double sphere_permeability = 500*air_permeability;
//   return 1.0;
// }

DevMaxwellStressTensorAux::DevMaxwellStressTensorAux(
    const std::string & f_gf_name,
    const std::string & f_coef_name,
    std::string b_gf_name,
    std::string h_gf_name,
    mfem::Array<int> boundary_attr)
  : BoundaryCoefficientAux(f_gf_name, f_coef_name, boundary_attr),
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


  coefficients._scalars.Register(_coef_name,
                                 std::make_shared<DevMaxwellStressTensorAuxCoefficient>(
                                     _b_gf, _h_gf));

  BoundaryCoefficientAux::Init(gridfunctions, coefficients);
}

} // namespace hephaestus
