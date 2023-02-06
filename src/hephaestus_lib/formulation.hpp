#pragma once
#include "../common/pfem_extras.hpp"
#include "auxkernels.hpp"
#include "equation_system.hpp"
#include "hephaestus_solvers.hpp"
#include "inputs.hpp"
#include "sources.hpp"

namespace hephaestus {
// Specifies output interfaces of a time-domain EM formulation.
class TransientFormulation : public mfem::TimeDependentOperator {

public:
  TransientFormulation(){};

  ~TransientFormulation(){};

  virtual void Init(mfem::Vector &X){};

  virtual void RegisterVariables() = 0;
  virtual void RegisterAuxKernels(hephaestus::AuxKernels &auxkernels){};

  mfem::Array<int> true_offsets;
};
} // namespace hephaestus
