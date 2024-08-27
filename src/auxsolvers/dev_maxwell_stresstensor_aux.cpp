#include "dev_maxwell_stresstensor_aux.hpp"
#include "utils.hpp"

#include <utility>

namespace hephaestus
{


void
DevMaxwellStressTensorAuxCoefficient::Eval(mfem::Vector & uxv,
                                                mfem::ElementTransformation & T,
                                                const mfem::IntegrationPoint & ip)
{
  double air_mu = 1.25663706e-6;
  double sphere_mu = 500*air_mu;
  double air_nu = 1.0/air_mu;
  double sphere_nu = 1.0/sphere_mu;

  std::cout << "in DevMaxwellStressTensorAuxCoefficient::Eval" << std::endl;

  _b_gf->GetVectorValue(T, ip, uxv);

  std::cout << "b GetVectorValue " << uxv(0) << " " << uxv(1) << std::endl;

}

DevMaxwellStressTensorAux::DevMaxwellStressTensorAux(
    const std::string & f_gf_name,
    const std::string & f_coef_name,
    std::string b_gf_name,
    std::string h_gf_name,
    mfem::Array<int> boundary_attr)
  : //VectorBoundaryNormalCoefficientAux(f_gf_name, f_coef_name, boundary_attr),
  // : BoundaryCoefficientAux(f_gf_name, f_coef_name, boundary_attr),
    _gf_name(std::move(f_gf_name)),
    _coef_name(std::move(f_coef_name)),
    _b_gf_name(std::move(b_gf_name)),
    _h_gf_name(std::move(h_gf_name)),
    _boundary_attr(boundary_attr)
{
}

void
DevMaxwellStressTensorAux::Init(const hephaestus::GridFunctions & gridfunctions,
                                      hephaestus::Coefficients & coefficients)
{
  _b_gf = gridfunctions.Get(_b_gf_name);
  _h_gf = gridfunctions.Get(_h_gf_name);

  _mesh_parent = _b_gf->ParFESpace()->GetParMesh();
  InitChildMesh();
  MakeFESpaces(0);
  MakeGridFunctions(0);

  // FIXME: This could be tricky part. Can we give a different name? i.e. is this coef the output field?
  coefficients._vectors.Register(_coef_name,
                                std::make_shared<DevMaxwellStressTensorAuxCoefficient>(
                                _b_gf_child));

  _gf = gridfunctions.Get(_gf_name);
  _vec_coef = coefficients._vectors.Get(_coef_name);
  _mass_coef = std::make_shared<mfem::ConstantCoefficient>(1.0);

  MakeFESpaces(1);
  MakeGridFunctions(1);
  
  std::cout << "Attributes max " << _mesh_child->bdr_attributes.Max() << std::endl;
  hephaestus::AttrToMarker(_boundary_attr, _boundary_attr_marker, _mesh_child->bdr_attributes.Max());
  std::cout << "_boundary_attr size = " << _boundary_attr.Size() << std::endl; 

  _test_fes = _gf_child->ParFESpace();

  BuildBilinearForm();
  BuildLinearForm();
  _a_mat = std::unique_ptr<mfem::HypreParMatrix>(_a->ParallelAssemble());
  _solver = std::make_unique<hephaestus::DefaultJacobiPCGSolver>(_solver_options, *_a_mat);
}


void
DevMaxwellStressTensorAux::BuildBilinearForm()
{
  mfem::ConstantCoefficient k(1.0);
  _a = std::make_unique<mfem::ParBilinearForm>(_test_fes);
  // _a = std::make_unique<mfem::ParMixedBilinearForm>(_test_fes, _trial_fes);
  // _a->AddDomainIntegrator(new mfem::MixedVectorProductIntegrator(*_vec_coef));
  // _a->AddDomainIntegrator(new mfem::DivDivIntegrator());//*_mass_coef));
  // _a->AddDomainIntegrator(new mfem::VectorFEMassIntegrator(k)); 
  // _a->AddBoundaryIntegrator(new mfem::VectorFEMassIntegrator(k)); 
  _a->AddDomainIntegrator(new mfem::MassIntegrator());
  _a->AddBoundaryIntegrator(new mfem::MassIntegrator()); // think this is correct
  // _a->AddBoundaryIntegrator(new mfem::BoundaryMassIntegrator(k));
  // _a->AddBdrFaceIntegrator(new mfem::BoundaryMassIntegrator(k));
  // _a->AddBoundaryIntegrator(new mfem::MixedVectorProductIntegrator(*_vec_coef)); // FIXME: need different _vec_coef
  _a->Assemble();
  _a->Finalize();
}

void
DevMaxwellStressTensorAux::BuildLinearForm()
{
  _b = std::make_unique<mfem::ParLinearForm>(_test_fes);
  // _b->AddBoundaryIntegrator(new mfem::VectorFEBoundaryFluxLFIntegrator(_scalar_coef), _boundary_attr_marker); // does not compile
  // _b->AddDomainIntegrator(new mfem::BoundaryNormalLFIntegrator(*_vec_coef));
  // _b->AddBoundaryIntegrator(new mfem::BoundaryNormalLFIntegrator(*_vec_coef), _boundary_attr_marker);
  _b->AddBoundaryIntegrator(new mfem::BoundaryNormalLFIntegrator(*_vec_coef), _boundary_attr_marker);
  // _b->AddBdrFaceIntegrator(new mfem::BoundaryNormalLFIntegrator(*_vec_coef), _boundary_attr_marker);
  // _b->AddDomainIntegrator(new mfem::VectorFEDomainLFIntegrator(*_vec_coef));
  _b->Assemble();
}

void
DevMaxwellStressTensorAux::Solve(double t)
{
  // when all this stuff is commented, there is no NaN
  /* */
  mfem::Vector x(_test_fes->GetTrueVSize()); // Gridfunction true DOFs
  // x = 0.0;

  std::cout << "GetTrueVSize " << _test_fes->GetTrueVSize() 
    << " _gf.Size() " << _gf->Size() << std::endl;

  if (_gf_child)
  {
    _gf_child->GetTrueDofs(x);
    std::cout 
      << "before transfer, _b_gf[0] = " << _b_gf->GetData()[0] << " "
      << "before transfer, _b_gf_child[0] = " << _b_gf_child->GetData()[0] << " "
      << std::endl;
    _mesh_child->Transfer(*_b_gf, *_b_gf_child);
    std::cout 
      << "after transfer, _b_gf[0] = " << _b_gf->GetData()[0] << " "
      << "after transfer, _b_gf_child[0] = " << _b_gf_child->GetData()[0] << " "
      << std::endl;
  }
  else
  {
    _gf->GetTrueDofs(x);
  }

  // Reassemble in case coef has changed
  _b->Update();
  _b->Assemble();

  std::cout << "Solve" << std::endl;
  std::cout 
    << "before solve, b_[0] = " << _b->GetData()[0] << " "
    << "before solve, x_[0] = " << x.GetData()[0] << " "
    << std::endl;
  _solver->Mult(*_b, x);
  
  if (_gf_child)
  {
    _gf_child->SetFromTrueDofs(x);
  }
  else
  {
    _gf->SetFromTrueDofs(x);
  }

  std::cout 
    << "after solve, b_[0] = " << _b->GetData()[0] << " "
    << "after solve, x_[0] = " << x.GetData()[0] << " "
    << std::endl;


  if (_gf)
    _mesh_child->Transfer(*_gf_child, *_gf);

  // Do IO
  std::ostringstream mesh_name, fes_name;
  int myid = mfem::Mpi::WorldRank();
  mesh_name << "mesh." << std::setfill('0') << std::setw(6) << myid;
  fes_name << "field." << std::setfill('0') << std::setw(6) << myid;
  std::ofstream mesh_ofs(mesh_name.str().c_str());
  std::ofstream fes_ofs(fes_name.str().c_str());
  _mesh_parent->Print(mesh_ofs);
  _gf->Save(fes_ofs);
}

void
DevMaxwellStressTensorAux::InitChildMesh()
{
  mfem::Array<int> domain_marker;
  domain_marker.Append(100);
  if (_mesh_child == nullptr)
  {
    _mesh_child = std::make_unique<mfem::ParSubMesh>(
        mfem::ParSubMesh::CreateFromDomain(*_mesh_parent, domain_marker));
        // mfem::ParSubMesh::CreateFromBoundary(*_mesh_parent, _boundary_attr));
  }
}


void
DevMaxwellStressTensorAux::MakeFESpaces(int stage)
{ 
  if (_h1_fe_space_child == nullptr && stage == 1)
  {
    std::cout << "Define _h1_child" << std::endl;
    int dim = _mesh_parent->Dimension();
    int dim_child = _mesh_child->Dimension();
    std::cout << "parent dim = " << dim << " child dim " << dim_child << std::endl;
    _order_h1 = _gf->ParFESpace()->FEColl()->GetOrder();
    _h1_fe_space_fec_child =
        std::make_unique<mfem::H1_FECollection>(_order_h1, dim_child);
    _h1_fe_space_child = std::make_shared<mfem::ParFiniteElementSpace>(
        _mesh_child.get(), _h1_fe_space_fec_child.get());
  }
  if (_h_div_fe_space_child == nullptr && stage == 0)
  {
    std::cout << "Define _h_div_child" << std::endl;
    _order_hdiv = _b_gf->ParFESpace()->FEColl()->GetOrder();
    _h_div_fe_space_fec_child =
        std::make_unique<mfem::RT_FECollection>(_order_hdiv - 1, _mesh_child->Dimension());
    _h_div_fe_space_child = std::make_shared<mfem::ParFiniteElementSpace>(
        _mesh_child.get(), _h_div_fe_space_fec_child.get());
  }
}

void 
DevMaxwellStressTensorAux::MakeGridFunctions(int stage)
{
  if (_gf_child == nullptr && stage == 1){
    std::cout << "setting _gf_child" << std::endl;
    _gf_child = std::make_shared<mfem::ParGridFunction>(_h1_fe_space_child.get());
  }
  if (_b_gf_child == nullptr && stage == 0){
    std::cout << "setting _b_gf_child" << std::endl;
    _b_gf_child = std::make_shared<mfem::ParGridFunction>(_h_div_fe_space_child.get());
    // *_b_gf_child = mfem::ParGridFunction(_h_div_fe_space_child.get());
  }

}

} // namespace hephaestus
