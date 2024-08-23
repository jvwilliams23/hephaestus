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
  // _mass_coef = coefficients._scalars.Get("one");
  _mass_coef = std::make_shared<mfem::ConstantCoefficient>(1.0);
  // _mass_coef = mfem::ConstantCoefficient(1.0);

  _mesh_parent = _gf->ParFESpace()->GetParMesh();
  hephaestus::AttrToMarker(_boundary_attr, _boundary_attr_marker, _mesh_parent->attributes.Max());
  
  InitChildMesh();
  MakeFESpaces();
  MakeGridFunctions();

  _test_fes = _gf_child->ParFESpace();

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
  _a->AddDomainIntegrator(new mfem::MassIntegrator(*_mass_coef));
  // _a->AddBoundaryIntegrator(new mfem::MassIntegrator(*_vec_coef)); 
  _a->Assemble();
  _a->Finalize();
}

void
VectorBoundaryNormalCoefficientAux::BuildLinearForm()
{
  _b = std::make_unique<mfem::ParLinearForm>(_test_fes);
  // _b->AddDomainIntegrator(new mfem::DomainLFIntegrator(*_vec_coef));
  // _b->AddBoundaryIntegrator(new mfem::BoundaryNormalLFIntegrator(*_vec_coef));
  // _b->AddBoundaryIntegrator(new mfem::BoundaryNormalLFIntegrator(*_vec_coef));
  _b->AddDomainIntegrator(new mfem::BoundaryNormalLFIntegrator(*_vec_coef));
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
    << " _gf.Size() " << _gf_child->Size() << std::endl;

  mfem::Array<int> new_boundary_markers;
  new_boundary_markers.Append(1);
  // TODO: Why are these the same size?
  // _gf_child->ProjectCoefficient(*_vec_coef);           // Initial condition
  _gf_child->ProjectBdrCoefficient(*_vec_coef, new_boundary_markers);
  // _gf_child->ProjectBdrCoefficient(*_coef, _boundary_attr_marker);
  _gf_child->GetTrueDofs(x);

  // Reassemble in case coef has changed
  _b->Update();
  _b->Assemble();

  _solver->Mult(*_b, x);

  _gf_child->SetFromTrueDofs(x);

  if (_gf)
    _mesh_child->Transfer(*_gf_child, *_gf);
  /* */

  std::ostringstream mesh_name, fes_name, sub_mesh_name, sub_fes_name;
  int myid = mfem::Mpi::WorldRank();
  mesh_name << "mesh." << std::setfill('0') << std::setw(6) << myid;
  fes_name << "field." << std::setfill('0') << std::setw(6) << myid;
  sub_mesh_name << "sub_mesh." << std::setfill('0') << std::setw(6) << myid;
  sub_fes_name << "sub_field." << std::setfill('0') << std::setw(6) << myid;
  std::ofstream mesh_ofs(mesh_name.str().c_str());
  std::ofstream sub_mesh_ofs(sub_mesh_name.str().c_str());
  std::ofstream fes_ofs(fes_name.str().c_str());
  std::ofstream sub_fes_ofs(sub_fes_name.str().c_str());
  _mesh_parent->Print(mesh_ofs);
  _mesh_child->Print(sub_mesh_ofs);
  _gf->Save(fes_ofs);
  _gf_child->Save(sub_fes_ofs);
  // for (int i=0; i<_gf_child->Size();++i)

}


void
VectorBoundaryNormalCoefficientAux::InitChildMesh()
{
  if (_mesh_child == nullptr)
  {
    _mesh_child = std::make_unique<mfem::ParSubMesh>(
        mfem::ParSubMesh::CreateFromBoundary(*_mesh_parent, _boundary_attr));
  }
}


void
VectorBoundaryNormalCoefficientAux::MakeFESpaces()
{ 
  if (_h1_fe_space_child == nullptr)
  {
    int dim = _mesh_parent->Dimension();
    int dim_child = _mesh_child->Dimension();
    std::cout << "parent dim = " << dim << " child dim " << dim_child << std::endl;
    _order_h1 = _gf->ParFESpace()->FEColl()->GetOrder();
    _h1_fe_space_fec_child =
        std::make_unique<mfem::H1_FECollection>(_order_h1, dim_child);
    _h1_fe_space_child = std::make_shared<mfem::ParFiniteElementSpace>(
        _mesh_child.get(), _h1_fe_space_fec_child.get());
  }

  // if (_h_curl_fe_space_child == nullptr)
  // {
  //   _h_curl_fe_space_fec_child =
  //       std::make_unique<mfem::ND_FECollection>(_order_hcurl, _mesh_child->Dimension());
  //   _h_curl_fe_space_child = std::make_shared<mfem::ParFiniteElementSpace>(
  //       _mesh_child.get(), _h_curl_fe_space_fec_child.get());
  // }

  // if (_source_current_density && _h_div_fe_space_child == nullptr)
  // {
  //   _h_div_fe_space_fec_child =
  //       std::make_unique<mfem::RT_FECollection>(_order_hdiv - 1, _mesh_child->Dimension());
  //   _h_div_fe_space_child = std::make_shared<mfem::ParFiniteElementSpace>(
  //       _mesh_child.get(), _h_div_fe_space_fec_child.get());
  // }
}

void 
VectorBoundaryNormalCoefficientAux::MakeGridFunctions()
{
  if (_gf_child == nullptr){
    std::cout << "setting _gf_child" << std::endl;
    _gf_child = std::make_shared<mfem::ParGridFunction>(_h1_fe_space_child.get());

  }

}

} // namespace hephaestus
