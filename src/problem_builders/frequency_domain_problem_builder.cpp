#include "frequency_domain_problem_builder.hpp"

namespace hephaestus {

void FrequencyDomainProblemBuilder::InitializeKernels() {
  this->problem->preprocessors.Init(this->problem->gridfunctions,
                                    this->problem->coefficients);
  this->problem->sources.Init(this->problem->gridfunctions,
                              this->problem->fespaces, this->problem->bc_map,
                              this->problem->coefficients);
}
void FrequencyDomainProblemBuilder::ConstructOperator() {
  this->problem->fd_operator =
      std::make_unique<hephaestus::FrequencyDomainEquationSystemOperator>(
          *(this->problem->pmesh), this->problem->fespaces,
          this->problem->gridfunctions, this->problem->bc_map,
          this->problem->coefficients, this->problem->sources,
          this->problem->solver_options);
  this->problem->fd_operator->SetGridFunctions();
}

void FrequencyDomainProblemBuilder::ConstructState() {
  this->problem->F = new mfem::BlockVector(
      this->problem->fd_operator->true_offsets); // Vector of dofs
  this->problem->fd_operator->Init(
      *(this->problem->F)); // Set up initial conditions
}
} // namespace hephaestus
