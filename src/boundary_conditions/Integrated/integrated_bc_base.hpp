#pragma once
#include "boundary_condition_base.hpp"

namespace hephaestus
{

class IntegratedBC : public BoundaryCondition
{
public:
  IntegratedBC(const std::string & name_, mfem::Array<int> bdr_attributes_);
  IntegratedBC(const std::string & name_,
               mfem::Array<int> bdr_attributes_,
               mfem::LinearFormIntegrator * lfi_re_,
               mfem::LinearFormIntegrator * lfi_im_ = nullptr);

  // NB: assume ownership of pointers.
  std::unique_ptr<mfem::LinearFormIntegrator> lfi_re;
  std::unique_ptr<mfem::LinearFormIntegrator> lfi_im;

  void applyBC(mfem::LinearForm & b) override;
  void applyBC(mfem::ComplexLinearForm & b) override;
  void applyBC(mfem::ParComplexLinearForm & b) override;
};

} // namespace hephaestus
