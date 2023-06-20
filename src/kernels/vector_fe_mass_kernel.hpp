#pragma once
#include "kernel_base.hpp"

namespace hephaestus {

/*
(βu, u')
*/
class VectorFEMassKernel : public Kernel<mfem::ParBilinearForm> {
public:
  VectorFEMassKernel(const hephaestus::InputParameters &params);
  virtual void
  Init(mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
       const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
       hephaestus::BCMap &bc_map,
       hephaestus::DomainProperties &domain_properties) override;
  virtual void Apply(mfem::ParBilinearForm *blf) override;
  std::string coef_name;
  mfem::Coefficient *coef;
};

}; // namespace hephaestus
