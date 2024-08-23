#pragma once
#include "vector_boundarynormal_coefficient_aux.hpp"

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
  const mfem::ParGridFunction * _b_gf{nullptr};
  const mfem::ParGridFunction * _h_gf{nullptr};

public:
  DevMaxwellStressTensorAuxCoefficient(const mfem::ParGridFunction * b_gf,
                                       const mfem::ParGridFunction * h_gf)
    : mfem::VectorCoefficient(3), _b_gf{b_gf}, _h_gf{h_gf}
    // : _b_gf{b_gf}, _h_gf{h_gf}
  {
  }

  ~DevMaxwellStressTensorAuxCoefficient() override = default;

  void Eval(mfem::Vector & uxv,
            mfem::ElementTransformation & T,
            const mfem::IntegrationPoint & ip) override;
  // double Eval(mfem::FaceElementTransformations & T, const mfem::IntegrationPoint & ip) override;
};

// Auxsolver to project the dot product of two vector gridfunctions onto a third
// (scalar) GridFunction
class DevMaxwellStressTensorAux : public VectorBoundaryNormalCoefficientAux
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

  ~DevMaxwellStressTensorAux() override = default;

  void Init(const hephaestus::GridFunctions & gridfunctions,
            hephaestus::Coefficients & coefficients) override;
};

} // namespace hephaestus
