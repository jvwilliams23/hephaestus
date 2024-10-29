#include "dev_maxwell_stresstensor_aux.hpp"
#include "utils.hpp"

#include <utility>

namespace hephaestus
{


void
DevMaxwellStressVectorAuxCoefficient::Eval(mfem::Vector & uxv,
                                                mfem::ElementTransformation & T,
                                                const mfem::IntegrationPoint & ip)
{
  double air_mu = 1.25663706e-6;
  double sphere_mu = 500*air_mu;
  double air_nu = 1.0/air_mu;
  double sphere_nu = 1.0/sphere_mu;

  std::cout << "in DevMaxwellStressVectorAuxCoefficient::Eval (vector)" << std::endl;

  _b_gf->GetVectorValue(T, ip, uxv);

  std::cout << "b GetVectorValue " 
    << uxv(0) << " " 
    << uxv(1) << " " 
    << uxv(2) << " " 
    << std::endl;
}

double
DevMaxwellStressScalarAuxCoefficient::Eval(mfem::ElementTransformation & T,
                                           const mfem::IntegrationPoint & ip)
{
  double air_mu = 1.25663706e-6;
  double sphere_mu = 500*air_mu;
  double air_nu = 1.0/air_mu;
  double sphere_nu = 1.0/sphere_mu;

  mfem::Vector local_dofs, normal_vec;
  mfem::Array<int> dof_ids;

  _b_gf->ParFESpace()->GetBdrElementDofs(T.ElementNo, dof_ids);
  _b_gf->GetSubVector(dof_ids, local_dofs);
  assert (local_dofs.Size() == 1);
  /*
  std::cout << "local_dofs " ;
  for (int j=0; j<local_dofs.Size();++j)
  {
    std::cout << local_dofs(j) << " ";
  }
  std::cout << std::endl;
  */ 
  return 0.5*(air_nu - sphere_nu)*local_dofs(0)*local_dofs(0);
  // return 1000.0;
}

double 
DevMaxwellStressScalarHFieldAuxCoefficient::Eval(mfem::ElementTransformation & T,
                                           const mfem::IntegrationPoint & ip)
{
  double air_mu = 1.25663706e-6;
  double sphere_mu = 500*air_mu;
  double air_nu = 1.0/air_mu;
  double sphere_nu = 1.0/sphere_mu;

  mfem::Vector local_dofs, normal_vec;
  mfem::Array<int> dof_ids;

  int space_dim = 3;
  // f_tr->Face->GetSpaceDim()

  _h_gf->ParFESpace()->GetBdrElementDofs(T.ElementNo, dof_ids);
  _h_gf->GetSubVector(dof_ids, local_dofs);
  assert (local_dofs.Size() == 1);
  // mfem::CalcOrtho(T.Face->Jacobian(), normal_vec);

  mfem::Vector h_vec(space_dim);
  mfem::Vector h_tangent_vec(space_dim);
  _h_gf->GetVectorValue(T,ip,h_vec);
  for (int k = 0; k<space_dim;++k)
  {
    h_tangent_vec(k) = h_vec(k) - normal_vec(k)*local_dofs(0);
  }
  double h_norm;
  h_norm = h_tangent_vec.Norml2();

  return 0.5*(air_mu - sphere_mu)*h_norm;
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

  _use_scalar_coef = true;

  // FIXME: This could be tricky part. Can we give a different name? i.e. is this coef the output field?
  _gf = gridfunctions.Get(_gf_name);
  if (_use_scalar_coef)
  {
    std::string _bcoef_name, _hcoef_name;
    _bcoef_name = "bcoef";
    _hcoef_name = "hcoef";
    coefficients._scalars.Register(_bcoef_name,
                                  std::make_shared<DevMaxwellStressScalarAuxCoefficient>(
                                  _b_gf_child));
    _scalar_b_coef = coefficients._scalars.Get(_bcoef_name);

    // coefficients._scalars.Register(_hcoef_name,
    //                               std::make_shared<DevMaxwellStressScalarHFieldAuxCoefficient>(
    //                               _h_gf_child));
    // _scalar_h_coef = coefficients._scalars.Get(_hcoef_name);
  }
  else
  {
    coefficients._vectors.Register(_coef_name,
                                  std::make_shared<DevMaxwellStressVectorAuxCoefficient>(
                                  _b_gf_child));
    _vec_coef = coefficients._vectors.Get(_coef_name);

  }

  _mass_coef = std::make_shared<mfem::ConstantCoefficient>(1.0);

  MakeFESpaces(1);
  MakeGridFunctions(1);
  
  std::cout << "Attributes max " << _mesh_child->bdr_attributes.Max() << std::endl;
  hephaestus::AttrToMarker(_boundary_attr, _boundary_attr_marker, _mesh_child->bdr_attributes.Max());
  std::cout << "_boundary_attr size = " << _boundary_attr.Size() << std::endl; 

  _test_fes_bf = _b_gf_child->ParFESpace();
  _trial_fes = _gf_child->ParFESpace();

  BuildBilinearForm();
  BuildLinearForm();
  _a_mat = std::unique_ptr<mfem::HypreParMatrix>(_a->ParallelAssemble());
  _solver = std::make_unique<hephaestus::DefaultJacobiPCGSolver>(_solver_options, *_a_mat);
}

inline void fFun(const mfem::Vector & x, mfem::Vector & f) {  f = 0.0; }

void
DevMaxwellStressTensorAux::BuildBilinearForm()
{
  mfem::Coefficient *alpha = new mfem::ConstantCoefficient(1.0);
  mfem::Coefficient *beta = new mfem::ConstantCoefficient(1.0);

  // MixedScalarWeakGradientIntegrator
  /* */
  // _a = std::make_unique<mfem::ParBilinearForm>(_test_fes_bf);
  // _a->AddDomainIntegrator(new mfem::VectorFEMassIntegrator(*beta)); 
  // _a = std::make_unique<mfem::ParMixedBilinearForm>(_trial_fes, _test_fes_bf);
  _a = std::make_unique<mfem::ParMixedBilinearForm>(_test_fes_bf, _trial_fes);
  mfem::VectorFunctionCoefficient fcoeff = mfem::VectorFunctionCoefficient(3, *fFun);
  _a->AddDomainIntegrator(new mfem::MixedVectorProductIntegrator(fcoeff)); 
  _a->Assemble();
  _a->Finalize();
}


void
DevMaxwellStressTensorAux::BuildLinearForm()
{
  // std::make_shared<mfem::ParGridFunction>
  // std::make_unique<mfem::ParGridFunction> j_full(_test_fes);
  // _jfunc = std::make_shared<mfem::ParGridFunction>(_test_fes);
  // mfem::ParGridFunction j_full(_test_fes);
  // mfem::VectorGridFunctionCoefficient fcoeff(_b_gf_child);
  // j_full = 0.0;

  // mfem::VectorConstantCoefficient fcoeff(mfem::Vector(3));
  // mfem::VectorFunctionCoefficient fcoeff(3, fFun);
  auto fcoeff = std::make_unique<mfem::VectorFunctionCoefficient>(3, *fFun);
  // mfem::VectorArrayCoefficient fcoeff(3);
  // _b = std::make_unique<mfem::ParLinearForm>(_test_fes_bf);
  _b = std::make_unique<mfem::ParLinearForm>(_test_fes_bf);
  // _h_div_fe_space_child gives b not Equal cols
  _b->AddBoundaryIntegrator(new mfem::VectorFEBoundaryFluxLFIntegrator(*_scalar_coef), _boundary_attr_marker); // does not compile
  // _b->AddDomainIntegrator(new mfem::VectorFEDomainLFIntegrator(*fcoeff.release())); // taken from ex5
  // _b->AddBoundaryIntegrator(new mfem::BoundaryNormalLFIntegrator(*_vec_coef), _boundary_attr_marker);
  _b->Assemble();
}

void
DevMaxwellStressTensorAux::Solve(double t)
{
  // when all this stuff is commented, there is no NaN
  /* */
  std::cout << "Starting solve" << std::endl;
  // mfem::Vector x(_test_fes->GetNDofs());
  mfem::ParGridFunction x(_trial_fes);
  x = 0.0;
  // mfem::Vector x(_h1_fe_space_child.get()->GetTrueVSize()); // Gridfunction true DOFs
  // x = 0.0;

  std::cout << "GetTrueVSize " << _test_fes_bf->GetTrueVSize() 
    << " _gf.Size() " << _gf_child->Size() 
    << " _h_div_fe_space_child.ndofs() " << _h_div_fe_space_child.get()->GetNDofs()
    << " _h1_fe_space_child.ndofs() " << _h1_fe_space_child.get()->GetNDofs()
    << " _h1_fe_space_child.truevsize() " << _h1_fe_space_child.get()->GetTrueVSize()
    << " x.size() " << x.Size()
    << " _b.size() " << _b->Size()
    << " _a.NumCols() " << _a->NumCols()
    << " _a.NumRows() " << _a->NumRows()
    << std::endl;

  if (_gf_child)
  {
    // _gf_child->GetTrueDofs(x);
    std::cout 
      << "before transfer, _b_gf[0] = " << _b_gf->GetData()[0] << " "
      << "before transfer, _b_gf_child[0] = " << _b_gf_child->GetData()[0] << " "
      << std::endl;
    _mesh_child->Transfer(*_b_gf, *_b_gf_child);
    _mesh_child->Transfer(*_h_gf, *_h_gf_child);
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
  // _a_mat->Mult(*_b, x);
  // _a_mat->Mult(x, *_b);
  // _solver->Mult(x, *_b); // shape is correct, but capacity not correct
  
  // mfem::GridFunction gftrial(_trial_fes);
  // _a_mat->Mult(gftrial, x);
  // _a_mat->Mult(x, gftrial);
  
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
  // _mesh_parent->Print(mesh_ofs);
  // _gf->Save(fes_ofs);
  _mesh_child->Print(mesh_ofs);
  _gf_child->Save(fes_ofs);
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
  if (_h_curl_fe_space_child == nullptr && stage == 0)
  {
    std::cout << "Define _h_curl_child" << std::endl;
    _order_hcurl = _h_gf->ParFESpace()->FEColl()->GetOrder();
    _h_curl_fe_space_fec_child =
      std::make_unique<mfem::ND_FECollection>(_order_hcurl, _mesh_child->Dimension());
    _h_curl_fe_space_child = std::make_shared<mfem::ParFiniteElementSpace>(
      _mesh_child.get(), _h_curl_fe_space_fec_child.get());
  }
}

void 
DevMaxwellStressTensorAux::MakeGridFunctions(int stage)
{
  if (_gf_child == nullptr && stage == 1){
    std::cout << "setting _gf_child" << std::endl;
    // if (_use_scalar_coef)
    // {
    //   _gf_child = std::make_shared<mfem::ParGridFunction>(_h_div_fe_space_child.get());
    // }
    // else
    // {
    // }
    // _gf_child = std::make_shared<mfem::ParGridFunction>(_h_div_fe_space_child.get());
    // _gf_child = std::make_shared<mfem::ParGridFunction>(_h_div_fe_space_child.get());
    _gf_child = std::make_shared<mfem::ParGridFunction>(_h1_fe_space_child.get());

  }
  if (_b_gf_child == nullptr && stage == 0){
    std::cout << "setting _b_gf_child" << std::endl;
    _b_gf_child = std::make_shared<mfem::ParGridFunction>(_h_div_fe_space_child.get());
    // *_b_gf_child = mfem::ParGridFunction(_h_div_fe_space_child.get());
  }
  if (_h_gf_child == nullptr && stage == 0)
  {
    std::cout << "setting _h_gf_child" << std::endl;
    _h_gf_child = std::make_shared<mfem::ParGridFunction>(_h_curl_fe_space_child.get());
  }

}

} // namespace hephaestus
