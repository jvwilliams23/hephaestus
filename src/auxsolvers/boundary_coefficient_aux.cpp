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
BoundaryCoefficientAux::BuildBilinearForm()
{
  _a = std::make_unique<mfem::ParBilinearForm>(_test_fes);
  // _a->AddBoundaryIntegrator(new mfem::MassIntegrator(*_coef)); 
  _a->AddDomainIntegrator(new mfem::MassIntegrator());
  _a->AddBoundaryIntegrator(new mfem::MassIntegrator(*_coef)); 
  /*if (_boundary_attr.Size() > 0)
  {
    std::cout << "constraining bilinearform boundary integration" << std::endl;
    // _a->AddBoundaryIntegrator(new mfem::MassIntegrator(*_coef), _boundary_attr_marker); 
    // _a->AddBoundaryIntegrator(new mfem::MassIntegrator(*_coef)); 
    _a->AddDomainIntegrator(new mfem::MassIntegrator(*_coef)); 
  }
  else
  {
    _a->AddBoundaryIntegrator(new mfem::MassIntegrator()); 
  }*/
  _a->Assemble();
  _a->Finalize();
}

void
BoundaryCoefficientAux::BuildLinearForm()
{
  _b = std::make_unique<mfem::ParLinearForm>(_test_fes);
  // _b->AddBoundaryIntegrator(new mfem::BoundaryLFIntegrator(*_coef));
  _b->AddDomainIntegrator(new mfem::DomainLFIntegrator(*_coef));
  /*if (_boundary_attr.Size() > 0)
  {
    std::cout << "constraining linearform boundary integration" << std::endl;
    // _b->AddBoundaryIntegrator(new mfem::BoundaryLFIntegrator(*_coef), _boundary_attr_marker);
    // _b->AddBoundaryIntegrator(new mfem::BoundaryLFIntegrator(*_coef));
    _b->AddDomainIntegrator(new mfem::DomainLFIntegrator(*_coef));
  }
  else
  {
    std::cout << "boundary attr size = 0, no constraint on integration" << std::endl;
    _b->AddBoundaryIntegrator(new mfem::BoundaryLFIntegrator(*_coef));
  }*/
  _b->Assemble();
}

void
BoundaryCoefficientAux::Solve(double t)
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


  mfem::Array<int> new_boundary_markers;
  new_boundary_markers.Append(1);
  // TODO: Why are these the same size?
  _gf_child->ProjectCoefficient(*_coef);           // Initial condition
  _gf_child->ProjectBdrCoefficient(*_coef, new_boundary_markers);
  // _gf_child->ProjectBdrCoefficient(*_coef, _boundary_attr_marker);
  _gf_child->GetTrueDofs(x);

  std::cout << "GetTrueVSize " << _test_fes->GetTrueVSize() 
    // << " EssentialTrueDofs size " << ess_dof.Size() 
    << " _gf.Size() " << _gf_child->Size() << std::endl;
  // Reassemble in case coef has changed
  _b->Update();
  _b->Assemble();

  // for (int i=0; i<_gf_child->Size();++i)
  // {
  //   std::cout << "init..   x[i="<<i<<"] = "<< x[i] << " "
  //     << " _b = " << _b->GetData()[i]
  //     << " _gf_child[i] " << _gf_child->GetData()[i]
  //     << std::endl;
  // }

  _solver->Mult(*_b, x);


  // for (int i = 0; i<_gf_child->Size(); ++i)
  // {
  //   std::cout << "BE = " << i 
  //     << " gf " << _gf_child->GetData()[i]
  // }
  _gf_child->SetFromTrueDofs(x);
  std::ostringstream mesh_name, port_name;
  int myid = mfem::Mpi::WorldRank();
  mesh_name << "port_mesh." << std::setfill('0') << std::setw(6) << myid;
  port_name << "port_mode." << std::setfill('0') << std::setw(6) << myid;
  std::ofstream mesh_ofs(mesh_name.str().c_str());
  std::ofstream port_ofs(port_name.str().c_str());
  _mesh_child->Print(mesh_ofs);
  _gf_child->Save(port_ofs);
  // for (int i=0; i<_gf_child->Size();++i)
  // {
  //   std::cout << "post-Mult.   x[i="<<i<<"] = "<< x[i] << " "
  //     << " _b = " << _b->GetData()[i]
  //     << " _gf_child[i] " << _gf_child->GetData()[i]
  //     << std::endl;
  // }

  if (_gf)
    _mesh_child->Transfer(*_gf_child, *_gf);
  /* */
}


void
BoundaryCoefficientAux::InitChildMesh()
{
  if (_mesh_child == nullptr)
  {
    // for (int i = 0; i < _attr_marker.Size(); ++i)
    //   std::cout << "_attr_marker[i="<<i<<"] "  << _attr_marker[i]  << std::endl; 
    _mesh_child = std::make_unique<mfem::ParSubMesh>(
        mfem::ParSubMesh::CreateFromBoundary(*_mesh_parent, _boundary_attr));
        // mfem::ParSubMesh::CreateFromDomain(*_mesh_parent, domain_marker));
  }
}


void
BoundaryCoefficientAux::MakeFESpaces()
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
BoundaryCoefficientAux::MakeGridFunctions()
{
  if (_gf_child == nullptr){
    std::cout << "setting _gf_child" << std::endl;
    _gf_child = std::make_shared<mfem::ParGridFunction>(_h1_fe_space_child.get());

  }

}

} // namespace hephaestus
