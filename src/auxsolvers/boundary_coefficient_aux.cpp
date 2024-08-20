#include "boundary_coefficient_aux.hpp"

#include <utility>

namespace hephaestus
{

BoundaryCoefficientAux::BoundaryCoefficientAux(std::string gf_name,
                               std::string coef_name,
                               std::vector<int> attr_marker_int,
                               hephaestus::InputParameters solver_options
                               )
  : _gf_name(std::move(gf_name)),
    _coef_name(std::move(coef_name)),
    _attr_marker_int(attr_marker_int),
    _solver_options(std::move(solver_options))
{
}

void
BoundaryCoefficientAux::Init(const hephaestus::GridFunctions & gridfunctions,
                     hephaestus::Coefficients & coefficients)
{
  _gf = gridfunctions.Get(_gf_name);
  _coef = coefficients._scalars.Get(_coef_name);
  // mfem::ConstantCoefficient dummy_coef(0.0);
  // _coef = &dummy_coef;
  _test_fes = _gf->ParFESpace();

  // const mfem::FiniteElementSpace * fes = _gf->FESpace();
  // std::cout << "mesh attributes size " << fes->GetMesh()->attributes.Size() << " max " << fes->GetMesh()->attributes.Max() << std::endl;
  // std::cout << "check ? stuff " << (fes->GetMesh()->attributes.Size() ? fes->GetMesh()->attributes.Max() : 0) << std::endl;
  
  for (auto & element : _attr_marker_int) 
  {
    if (element != -1)
      _attr_marker.Append(element);
  }

  BuildBilinearForm();
  BuildLinearForm();
  _a_mat = std::unique_ptr<mfem::HypreParMatrix>(_a->ParallelAssemble());
  _solver = std::make_unique<hephaestus::DefaultJacobiPCGSolver>(_solver_options, *_a_mat);
}

void
BoundaryCoefficientAux::BuildBilinearForm()
{
  _a = std::make_unique<mfem::ParBilinearForm>(_test_fes);
  // _a->AddDomainIntegrator(new mfem::MassIntegrator()); // segfault in LinearForm
  _a->AddBoundaryIntegrator(new mfem::MassIntegrator(), _attr_marker); // segfault in LinearForm
  // _a->AddBoundaryIntegrator(new mfem::BoundaryMassIntegrator()); // does not compile
  // _a->AddDomainIntegrator(new mfem::MassIntegrator()); // segfault in LinearForm
  // _a->AddBoundaryIntegrator(new mfem::MassIntegrator(*_coef)); // segfault in BilinearForm
  // _a->AddBoundaryIntegrator(new mfem::BoundaryMassIntegrator()); // does not compile
  _a->Assemble();
  _a->Finalize();
}

void
BoundaryCoefficientAux::BuildLinearForm()
{
  _b = std::make_unique<mfem::ParLinearForm>(_test_fes);
  // _b = std::make_unique<mfem::ParLinearForm>(_test_fes;
  // _b = _test_fes;
  // _b->AddDomainIntegrator(new mfem::DomainLFIntegrator(*_coef));
  // _b->AddDomainIntegrator(new mfem::BoundaryLFIntegrator(*_coef));
  if (_attr_marker.Size() > 0)
  {
    _b->AddBoundaryIntegrator(new mfem::BoundaryLFIntegrator(*_coef), _attr_marker);
  }
  else
  {
    _b->AddBoundaryIntegrator(new mfem::BoundaryLFIntegrator(*_coef));
  }

  // _b.AddBoundaryIntegrator(new mfem::BoundaryLFIntegrator(*_coef));
  // _b->AddBoundaryIntegrator(new mfem::BoundaryNormalLFIntegrator(*_coef));
  /*if (_attr_marker.Size() > 0)
  {
    _b->AddBoundaryIntegrator(new mfem::BoundaryLFIntegrator(*_coef), _attr_marker);
  }
  else
  {
    _b->AddBoundaryIntegrator(new mfem::BoundaryLFIntegrator(*_coef));
  }
  */
  /*std::cout << "_attr_marker size " << _attr_marker.Size() << std::endl;
  if (_attr_marker.Size() > 0)
  {
    std::cout << "_attr_marker[0] " << _attr_marker[0] << std::endl;
    _b->AddDomainIntegrator(new mfem::DomainLFIntegrator(*_coef), _attr_marker);
  }
  else
  {
    _b->AddDomainIntegrator(new mfem::DomainLFIntegrator(*_coef));
  }*/
  _b->Assemble();
  // _b.Assemble();
}

void
BoundaryCoefficientAux::Solve(double t)
{
  /**/
  mfem::Vector x(_test_fes->GetTrueVSize()); // Gridfunction true DOFs
  _gf->ProjectCoefficient(*_coef);           // Initial condition
  _gf->GetTrueDofs(x);

  // Reassemble in case coef has changed
  _b->Update();
  _b->Assemble();

  _solver->Mult(*_b, x);
  _gf->SetFromTrueDofs(x);
  /**/
}

} // namespace hephaestus
