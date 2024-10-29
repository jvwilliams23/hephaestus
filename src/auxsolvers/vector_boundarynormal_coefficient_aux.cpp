#include "vector_boundarynormal_coefficient_aux.hpp"
#include "utils.hpp"

#include <utility>

namespace hephaestus
{

VectorBoundaryNormalCoefficientAux::VectorBoundaryNormalCoefficientAux(std::string gf_name,
                               std::string coef_name,
                               mfem::Array<int> boundary_attr,
                               hephaestus::InputParameters solver_options
                               )
  : _gf_name(std::move(gf_name)),
    _vec_coef_name(std::move(coef_name)),
    _boundary_attr(boundary_attr),
    _solver_options(std::move(solver_options))
{
}

void
VectorBoundaryNormalCoefficientAux::Init(const hephaestus::GridFunctions & gridfunctions,
                     hephaestus::Coefficients & coefficients)
{
  _gf = gridfunctions.Get(_gf_name);
  _vec_coef = coefficients._vectors.Get(_vec_coef_name);
  // _hcurl_coef = coefficients..Get(_vec_coef_name);
  _rt_boundary_coef = std::make_shared<mfem::ConstantCoefficient>(1.0);
  // _mass_coef = coefficients._scalars.Get("one");
  _mass_coef = std::make_shared<mfem::ConstantCoefficient>(1.0);
  // _mass_coef = mfem::ConstantCoefficient(1.0);

  _mesh_parent = _gf->ParFESpace()->GetParMesh();
  
  std::cout << "Attributes max " << _mesh_parent->bdr_attributes.Max() << std::endl;
  hephaestus::AttrToMarker(_boundary_attr, _boundary_attr_marker, _mesh_parent->bdr_attributes.Max());

  
  // InitChildMesh();
  // MakeFESpaces();
  // MakeGridFunctions();

  // _test_fes = _gf_child->ParFESpace();
  _test_fes = _gf->ParFESpace();

  BuildBilinearForm();
  BuildLinearForm();
  _a_mat = std::unique_ptr<mfem::HypreParMatrix>(_a->ParallelAssemble());
  _solver = std::make_unique<hephaestus::DefaultJacobiPCGSolver>(_solver_options, *_a_mat);
}

void
VectorBoundaryNormalCoefficientAux::BuildBilinearForm()
{
  _a = std::make_unique<mfem::ParBilinearForm>(_test_fes);
  // _a->AddDomainIntegrator(new mfem::MixedVectorProductIntegrator(*_vec_coef));
  _a->AddDomainIntegrator(new mfem::DivDivIntegrator());//*_mass_coef));
  // _a->AddBoundaryIntegrator(new mfem::MassIntegrator()); 
  _a->Assemble();
  _a->Finalize();
}

void
VectorBoundaryNormalCoefficientAux::BuildLinearForm()
{
  _b = std::make_unique<mfem::ParLinearForm>(_test_fes);
  // _b->AddBoundaryIntegrator(new mfem::VectorFEBoundaryFluxLFIntegrator(), _boundary_attr_marker);
  _b->AddDomainIntegrator(new mfem::VectorFEDomainLFIntegrator(*_vec_coef));
  _b->Assemble();
}

void
VectorBoundaryNormalCoefficientAux::Solve(double t)
{
  // when all this stuff is commented, there is no NaN
  /* */
  mfem::Vector x(_test_fes->GetTrueVSize()); // Gridfunction true DOFs
  // x = 0.0;

  //  Array<int> ess_bdr(pmesh->bdr_attributes.Max());
  //  ess_bdr = 1;
  //  mfem::Array<int> ess_dof;
  //  _test_fes->GetEssentialTrueDofs(_attr_marker, ess_dof);
  // mfem::Vector x(ess_dof.Size()); // Gridfunction true DOFs

  std::cout << "GetTrueVSize " << _test_fes->GetTrueVSize() 
    << " _gf.Size() " << _gf->Size() << std::endl;

  // mfem::Array<int> new_boundary_markers;
  // new_boundary_markers.Append(1);
  // TODO: Why are these the same size?
  std::cout << "Set initial cond" << std::endl;
  // _gf_child->ProjectCoefficient(*_vec_coef);           // Initial condition
  _gf->ProjectBdrCoefficient(*_vec_coef, _boundary_attr_marker);
  // _gf_child->ProjectBdrCoefficient(*_coef, _boundary_attr_marker);
  std::cout << "Set x" << std::endl;
  _gf->GetTrueDofs(x);

  // Reassemble in case coef has changed
  _b->Update();
  _b->Assemble();

  std::cout << "Solve" << std::endl;
  _solver->Mult(*_b, x);

  _gf->SetFromTrueDofs(x);


  std::ostringstream mesh_name, fes_name, sub_mesh_name, sub_fes_name;
  int myid = mfem::Mpi::WorldRank();
  mesh_name << "mesh." << std::setfill('0') << std::setw(6) << myid;
  fes_name << "field." << std::setfill('0') << std::setw(6) << myid;
  // sub_mesh_name << "sub_mesh." << std::setfill('0') << std::setw(6) << myid;
  // sub_fes_name << "sub_field." << std::setfill('0') << std::setw(6) << myid;
  std::ofstream mesh_ofs(mesh_name.str().c_str());
  // std::ofstream sub_mesh_ofs(sub_mesh_name.str().c_str());
  std::ofstream fes_ofs(fes_name.str().c_str());
  // std::ofstream sub_fes_ofs(sub_fes_name.str().c_str());
  _mesh_parent->Print(mesh_ofs);
  // _mesh_child->Print(sub_mesh_ofs);
  _gf->Save(fes_ofs);
  // _gf_child->Save(sub_fes_ofs);
  // for (int i=0; i<_gf_child->Size();++i)

}

} // namespace hephaestus
