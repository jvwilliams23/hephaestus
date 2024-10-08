#pragma once
#include "kernel_base.hpp"

namespace hephaestus
{

/*
(σ ∇ V, u')
*/
class MixedVectorGradientKernel : public Kernel<mfem::ParMixedBilinearForm>
{
public:
  MixedVectorGradientKernel(const hephaestus::InputParameters & params);

  ~MixedVectorGradientKernel() override = default;

  void Init(hephaestus::GridFunctions & gridfunctions,
            const hephaestus::FESpaces & fespaces,
            hephaestus::BCMap & bc_map,
            hephaestus::Coefficients & coefficients) override;
  void Apply(mfem::ParMixedBilinearForm * mblf) override;
  std::string _coef_name;
  mfem::Coefficient * _coef{nullptr};
};

} // namespace hephaestus
