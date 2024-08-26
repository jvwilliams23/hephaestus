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
  // mfem::ConstantCoefficient air_mu_coef(air_mu);
  // mfem::ConstantCoefficient sphere_mu_coef(sphere_mu);
  // mfem::ConstantCoefficient air_nu_coef(air_nu);
  // mfem::ConstantCoefficient sphere_nu_coef(sphere_nu);
  // std::cout << "in Eval devMaxwell" << std::endl;

  // mfem::ConstantCoefficient nu_diff (air_nu - sphere_nu);
  // mfem::ConstantCoefficient mu_diff (air_mu - sphere_mu);
  // mfem::Coefficient b_gf_val;
  _b_gf->GetVectorValue(T, ip, uxv);
  uxv *= (air_mu - sphere_mu);

  // mfem::ProductCoefficient b_gf_sqr(b_gf_val, b_gf_val);
  // mfem::Vector b_gf_sqr = b_gf_val * b_gf_val;
  // mfem::ProductCoefficient b_gf_val(b_gf_val, nu_diff);
  // std::cout << "b GetVectorValue " << b_gf_val(0) << " " << b_gf_val(1) << std::endl;
  // mfem::Vector h_gf_val;
  // _h_gf->GetVectorValue(T, ip, h_gf_val);
  // // mfem::ProductCoefficient h_gf_val(h_gf_val, h_gf_val);
  // // mfem::ProductCoefficient h_gf_val(h_gf_val, mu_diff);
  // std::cout << "h GetVectorValue " << h_gf_val(0) << " " << h_gf_val(1) << std::endl;
  // uxv = b_gf_val - h_gf_val;
}

// double
// DevMaxwellStressTensorAuxCoefficient::Eval(mfem::ElementTransformation & T,
//                                            const mfem::IntegrationPoint & ip)
// {
//   double air_mu = 1.25663706e-6;
//   double sphere_mu = 500*air_mu;
//   std::cout << "in Eval devMaxwell" << std::endl;

//   mfem::Vector _b_gf_val;
//   mfem::Vector _h_gf_val;
//   _b_gf->GetVectorValue(T, ip, _b_gf_val);
//   std::cout << "GetVectorValue " << _b_gf_val(0) << _b_gf_val(1) << std::endl;
//   // return _b_gf_val;
//   return 1000.0;
// }

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
    _h_gf_name(std::move(h_gf_name))
{
}

void
DevMaxwellStressTensorAux::Init(const hephaestus::GridFunctions & gridfunctions,
                                      hephaestus::Coefficients & coefficients)
{
  _b_gf = gridfunctions.Get(_b_gf_name);
  _h_gf = gridfunctions.Get(_h_gf_name);


  coefficients._vectors.Register(_coef_name,
                                 std::make_shared<DevMaxwellStressTensorAuxCoefficient>(
                                     _b_gf, _h_gf));

  _gf = gridfunctions.Get(_gf_name);
  _vec_coef = coefficients._vectors.Get(_coef_name);
  // _hcurl_coef = coefficients..Get(_coef_name);
  _rt_boundary_coef = std::make_shared<mfem::ConstantCoefficient>(1.0);
  // _mass_coef = coefficients._scalars.Get("one");
  _mass_coef = std::make_shared<mfem::ConstantCoefficient>(1.0);
  // _mass_coef = mfem::ConstantCoefficient(1.0);

  _mesh_parent = _gf->ParFESpace()->GetParMesh();
  
  std::cout << "Attributes max " << _mesh_parent->bdr_attributes.Max() << std::endl;
  hephaestus::AttrToMarker(_boundary_attr, _boundary_attr_marker, _mesh_parent->bdr_attributes.Max());

  _test_fes = _gf->ParFESpace();

  BuildBilinearForm();
  BuildLinearForm();
  _a_mat = std::unique_ptr<mfem::HypreParMatrix>(_a->ParallelAssemble());
  _solver = std::make_unique<hephaestus::DefaultJacobiPCGSolver>(_solver_options, *_a_mat);

}


void
DevMaxwellStressTensorAux::BuildBilinearForm()
{
  _a = std::make_unique<mfem::ParBilinearForm>(_test_fes);
  // _a->AddDomainIntegrator(new mfem::MixedVectorProductIntegrator(*_vec_coef));
  _a->AddDomainIntegrator(new mfem::DivDivIntegrator());//*_mass_coef));
  // _a->AddBoundaryIntegrator(new mfem::MassIntegrator()); 
  _a->Assemble();
  _a->Finalize();
}

void
DevMaxwellStressTensorAux::BuildLinearForm()
{
  _b = std::make_unique<mfem::ParLinearForm>(_test_fes);
  // _b->AddBoundaryIntegrator(new mfem::VectorFEBoundaryFluxLFIntegrator(), _boundary_attr_marker);
  _b->AddDomainIntegrator(new mfem::VectorFEDomainLFIntegrator(*_vec_coef));
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

  _gf->ProjectBdrCoefficient(*_vec_coef, _boundary_attr_marker);
  _gf->GetTrueDofs(x);

  // Reassemble in case coef has changed
  _b->Update();
  _b->Assemble();

  std::cout << "Solve" << std::endl;
  _solver->Mult(*_b, x);

  _gf->SetFromTrueDofs(x);

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

} // namespace hephaestus
