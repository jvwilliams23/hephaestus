#pragma once
#include "kernel_base.hpp"

namespace hephaestus {

/*
(α∇×u, ∇×u')
*/
class CurlCurlKernel : public Kernel<mfem::ParBilinearForm> {
public:
  CurlCurlKernel(const hephaestus::InputParameters &params);
  virtual void
  Init(mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
       const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
       hephaestus::BCMap &bc_map,
       hephaestus::Coefficients &domain_properties) override;
  virtual void Apply(mfem::ParBilinearForm *blf) override;
  std::string coef_name;
  mfem::Coefficient *coef;
};

}; // namespace hephaestus
