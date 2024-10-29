#include "boundary_coefficient_aux.hpp"
#include "utils.hpp"

#include <utility>

namespace hephaestus
{

BoundaryCoefficientAux::BoundaryCoefficientAux(std::string gf_name,
                               std::string coef_name,
                               mfem::Array<int> boundary_attr,
                               hephaestus::InputParameters solver_options
                               )
  : _gf_name(std::move(gf_name)),
    _coef_name(std::move(coef_name)),
    _boundary_attr(boundary_attr),
    _solver_options(std::move(solver_options))
{
}

void
BoundaryCoefficientAux::Init(const hephaestus::GridFunctions & gridfunctions,
                     hephaestus::Coefficients & coefficients)
{
  _gf = gridfunctions.Get(_gf_name);
  _coef = coefficients._scalars.Get(_coef_name);

  _mesh_parent = _gf->ParFESpace()->GetParMesh();
  hephaestus::AttrToMarker(_boundary_attr, _boundary_attr_marker, _mesh_parent->attributes.Max());

  _test_fes = _gf->ParFESpace();

  BuildBilinearForm();
  BuildLinearForm();
  _a_mat = std::unique_ptr<mfem::HypreParMatrix>(_a->ParallelAssemble());
  _solver = std::make_unique<hephaestus::DefaultJacobiPCGSolver>(_solver_options, *_a_mat);
}

void
BoundaryCoefficientAux::BuildBilinearForm()
{
  _a = std::make_unique<mfem::ParBilinearForm>(_test_fes);
  if (_test_fes->FEColl()->GetRangeType(3) == mfem::FiniteElement::SCALAR)
  {
    std::cout << "IsSCALAR" << std::endl;
  }
  else
  {
    std::cout << "IsVECTOR" << std::endl;
  }
  
  _a->AddDomainIntegrator(new mfem::MassIntegrator());
  // _a->AddBoundaryIntegrator(new mfem::MassIntegrator(*_coef), _boundary_attr_marker); 
  _a->AddBoundaryIntegrator(new mfem::MassIntegrator()); 
  _a->Assemble();
  _a->Finalize();
}

void
BoundaryCoefficientAux::BuildLinearForm()
{
  _b = std::make_unique<mfem::ParLinearForm>(_test_fes);
  // _b->AddBoundaryIntegrator(new mfem::VectorFEBoundaryFluxLFIntegrator(*_coef), _boundary_attr_marker);
  _b->AddBoundaryIntegrator(new mfem::VectorFEBoundaryFluxLFIntegrator(*_coef));
  // if (_test_fes->FEColl()->GetRangeType(3) == mfem::FiniteElement::SCALAR)
  // {
  //   _b->AddDomainIntegrator(new mfem::VectorDomainLFIntegrator(*_coef));
  // }
  // else
  // {
  //   _b->AddDomainIntegrator(new mfem::VectorFEDomainLFIntegrator(*_coef));
  // }
  _b->Assemble();
}

void
BoundaryCoefficientAux::Solve(double t)
{
  mfem::Vector x(_test_fes->GetTrueVSize()); // Gridfunction true DOFs
  x = 0.0;

  _gf->ProjectBdrCoefficient(*_coef, _boundary_attr_marker);
  _gf->GetTrueDofs(x);

  std::cout << "GetTrueVSize " << _test_fes->GetTrueVSize() 
    << " _gf.Size() " << _gf->Size() << std::endl;

  // Reassemble in case coef has changed
  _b->Update();
  _b->Assemble();

  _solver->Mult(*_b, x);

  _gf->SetFromTrueDofs(x);
}

} // namespace hephaestus
