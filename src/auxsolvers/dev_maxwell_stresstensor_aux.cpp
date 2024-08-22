#include "dev_maxwell_stresstensor_aux.hpp"

#include <utility>

namespace hephaestus
{

double
DevMaxwellStressTensorAuxCoefficient::Eval(mfem::ElementTransformation & T,
                                        const mfem::IntegrationPoint & ip)
{
  // double coef_value;
  // coef_value = _coef.Eval(T, ip);
  // mfem::Vector b;
  // mfem::Vector h;

  // _b_gf->GetVectorValue(T, ip, b);
  // _h_gf->GetVectorValue(T, ip, h);

  double air_permeability = 1.25663706e-6;
  double sphere_permeability = 500*air_permeability;

  // F = [b_n^2 * (1/\mu_0 - 1/\mu) - h_t^2 * (\mu_0 - \mu)] * n * 1/2
  // return 0.5 * (
  //   ((_b*_b) * (1.0/air_permeability - 1.0/sphere_permeability)) 
  //   - ((_h*_h) * (air_permeability - sphere_permeability))
  // )*;

  // const mfem::FiniteElementSpace * fes = _b_gf->FESpace();
  // mfem::Mesh * mesh = fes->GetMesh();
  
  // mfem::Vector local_dofs, normal_vec;
  // mfem::DenseMatrix dshape;
  // mfem::Array<int> dof_ids;

  // int face_attr(101); // TODO: Pass as arg to Eval

  // for (int i = 0; i < mesh->GetNBE(); i++)
  // {

  //   if (mesh->GetBdrAttribute(i) != face_attr)
  //     continue;

  //   mfem::FaceElementTransformations * f_tr =
  //       mesh->GetFaceElementTransformations(mesh->GetBdrElementFaceIndex(i));
  //   if (f_tr == nullptr)
  //     continue;

  //   const mfem::FiniteElement & elem = *fes->GetFE(f_tr->Elem1No);
  //   f_tr->Attribute = mesh->GetAttribute(f_tr->Elem1No);
  //   const int int_order = 2 * elem.GetOrder() + 3;
  //   const mfem::IntegrationRule & ir = mfem::IntRules.Get(f_tr->FaceGeom, int_order);

  //   fes->GetElementDofs(f_tr->Elem1No, dof_ids);
  //   _b_gf->GetSubVector(dof_ids, local_dofs);
  //   const int space_dim = f_tr->Face->GetSpaceDim();
  //   normal_vec.SetSize(space_dim);
  //   dshape.SetSize(elem.GetDof(), space_dim);
  // }

  // mfem::Vector
  // return _b_gf->GetValue(T, ip); 
  return 1.0;
}

DevMaxwellStressTensorAux::DevMaxwellStressTensorAux(
    const std::string & f_gf_name,
    const std::string & f_coef_name,
    std::string b_gf_name,
    std::string h_gf_name,
    mfem::Array<int> boundary_attr)
  : BoundaryCoefficientAux(f_gf_name, f_coef_name, boundary_attr),
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


  coefficients._scalars.Register(_coef_name,
                                 std::make_shared<DevMaxwellStressTensorAuxCoefficient>(
                                     _b_gf, _h_gf));

  BoundaryCoefficientAux::Init(gridfunctions, coefficients);
}

} // namespace hephaestus
