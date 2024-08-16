#include "new_maxwell_stresstensor_aux.hpp"

namespace hephaestus
{
double
calcMaxwellStressTensor(mfem::GridFunction * b_field, mfem::GridFunction * h_field, int face_attr, mfem::Coefficient & q)
{

  double flux = 0.0;
  double force = 0.0;
  double total_force = 0.0;
  double force_density = 0.0;
  double area = 0.0;

  double air_permeability = M_PI * 4.0e-7; //1.25663706e-6;
  // double air_permeability = 1.0;
  // double air_permeability = 1.25663706e-6;
  double sphere_permeability = 500*air_permeability;

  mfem::FiniteElementSpace * b_fes = b_field->FESpace();
  mfem::Mesh * mesh = b_fes->GetMesh(); // Should be same for B and H?
  mfem::FiniteElementSpace * h_fes = h_field->FESpace();
  // mfem::Mesh * mesh = b_fes->GetMesh(); // Should be same for B and H?

  mfem::Vector b_local_dofs, h_local_dofs, normal_vec, tangent_vec;
  mfem::DenseMatrix b_dshape, h_dshape;
  mfem::Array<int> b_dof_ids, h_dof_ids;

  mfem::Vector ones;
  
  std::cout << "Looping GetNBE " << std::endl;
  for (int i = 0; i < mesh->GetNBE(); i++)
  {

    if (mesh->GetBdrAttribute(i) != face_attr)
      continue;

    mfem::FaceElementTransformations * f_tr =
        mesh->GetFaceElementTransformations(mesh->GetBdrElementFaceIndex(i));

    if (f_tr == nullptr)
      continue;

    const mfem::FiniteElement & b_elem = *b_fes->GetFE(f_tr->Elem1No);
    const mfem::FiniteElement & h_elem = *h_fes->GetFE(f_tr->Elem1No);
    f_tr->Attribute = mesh->GetAttribute(f_tr->Elem1No);
    const int int_order = 2 * b_elem.GetOrder() + 3;
    // const int h_int_order = 2 * h_elem.GetOrder() + 3;
    std::cout << "Get integration rule " << std::endl;
    const mfem::IntegrationRule & ir = mfem::IntRules.Get(f_tr->FaceGeom, int_order);

    b_fes->GetElementDofs(f_tr->Elem1No, b_dof_ids);
    h_fes->GetElementDofs(f_tr->Elem1No, h_dof_ids);
    b_field->GetSubVector(b_dof_ids, b_local_dofs);
    h_field->GetSubVector(h_dof_ids, h_local_dofs);
    const int space_dim = f_tr->Face->GetSpaceDim();
    normal_vec.SetSize(space_dim);
    tangent_vec.SetSize(space_dim);
    
    b_dshape.SetSize(b_elem.GetDof(), space_dim);
    h_dshape.SetSize(h_elem.GetDof(), space_dim);

    double force_i = 0.0;
    double area_i = 0.0;
    double force_density_i = 0.0;

    for (int j = 0; j < ir.GetNPoints(); j++)
    {
      std::cout << "\n\nLoop over NPoints " << j << std::endl;
      double force_j(0.0);
      const mfem::IntegrationPoint & ip = ir.IntPoint(j);
      mfem::IntegrationPoint eip;
      f_tr->Loc1.Transform(ip, eip);
      f_tr->Face->SetIntPoint(&ip);
      double face_weight = f_tr->Face->Weight();
      // double b_normal_val = 0.0;
      // double h_tangent_val = 0.0;
      f_tr->Elem1->SetIntPoint(&eip);
      b_elem.CalcVShape(*f_tr->Elem1, b_dshape);
      h_elem.CalcVShape(*f_tr->Elem1, h_dshape);
      mfem::CalcOrtho(f_tr->Face->Jacobian(), normal_vec);
      normal_vec /= face_weight;

      /*
      mfem::Vector b_normal_vec_unit(space_dim);
      mfem::Vector b_normal_vec(space_dim);
      mfem::Vector h_normal_vec_unit(space_dim);
      mfem::Vector h_normal_vec(space_dim);
      mfem::Vector h_tangent_vec(space_dim);

      // get b unit normal vector, and normal component (value) 
      b_normal_vec_unit = b_local_dofs;
      b_normal_vec_unit /= b_local_dofs.Norml2();
      b_normal_vec = b_normal_vec_unit;
      b_normal_vec *= b_local_dofs;
      double b_normal_val = b_normal_vec.Norml2();

      // get h unit normal and tangential
      h_normal_vec_unit = h_local_dofs;
      h_normal_vec_unit /= h_local_dofs.Norml2();
      h_normal_vec = h_normal_vec_unit;
      h_normal_vec *= h_local_dofs;
      h_tangent_vec = h_local_dofs - h_normal_vec;
      double h_tangent_val = h_tangent_vec.Norml2();
      */
      /**/
      double b_normal_val_j = b_dshape.InnerProduct(normal_vec, b_local_dofs);
      double b_normal_val = b_normal_val_j;
      
      // find h_{t} (tangential part) as norm(h - h_\perp * n)
      // double h_normal_val_j = 0.0;
      double h_normal_val_j = h_dshape.InnerProduct(normal_vec, h_local_dofs);
      mfem::Vector h_normal_vec(space_dim);
      mfem::Vector h_tangent_vec(space_dim);
      // h_normal_vec = normal_vec;
      // h_normal_vec *= h_normal_val_j;
      for (int k = 0; k<space_dim;++k)
      {
        h_normal_vec(k) = normal_vec(k) * h_normal_val_j;
        std::cout << "k = " << k << " h_normal_vec(k): " << h_normal_vec(k) << " normal_vec(k): " << normal_vec(k) << std::endl;
        h_tangent_vec(k) = h_local_dofs(k) - h_normal_vec(k);
      }
      double h_tangent_val = h_tangent_vec.Norml2();
      // // get tangent vec
      // mfem::Vector tangent_vec(space_dim);
      // tangent_vec = h_tangent_vec;
      // tangent_vec /= h_tangent_vec.Norml2();
      // h_tangent_val += h_dshape.ProjectBdrCoefficientTangent(h_local_dofs)
      /**/
      std::cout << "normal_vec ([" 
        << normal_vec.GetData()[0] << ", " 
        << normal_vec.GetData()[1] << ", " 
        << normal_vec.GetData()[2] << "]) " 
        << std::endl;
      std::cout << "b_local_dofs ([" 
        << b_local_dofs.GetData()[0] << ", " 
        << b_local_dofs.GetData()[1] << ", " 
        << b_local_dofs.GetData()[2] << "]) " 
        << std::endl;
      std::cout << "h_local_dofs ([" 
        << h_local_dofs.GetData()[0] << ", " 
        << h_local_dofs.GetData()[1] << ", " 
        << h_local_dofs.GetData()[2] << "]) " 
        << std::endl;
      std::cout << "h_normal_vec ([" 
        << h_normal_vec.GetData()[0] << ", " 
        << h_normal_vec.GetData()[1] << ", " 
        << h_normal_vec.GetData()[2] << "]) " 
        << std::endl;
      std::cout << "h_tangent_vec ([" 
        << h_tangent_vec.GetData()[0] << ", " 
        << h_tangent_vec.GetData()[1] << ", " 
        << h_tangent_vec.GetData()[2] << "]) " 
        << std::endl;
      /*
      std::cout << "h_normal_vec " 
        << h_normal_vec.GetData()[0] << " " 
        << h_normal_vec.GetData()[1] << " " 
        << h_normal_vec.GetData()[2] << " " 
        << std::endl;
      std::cout << "h_tangent_vec " 
        << h_tangent_vec.GetData()[0] << " " 
        << h_tangent_vec.GetData()[1] << " " 
        << h_tangent_vec.GetData()[2] << " " 
        << std::endl;
      */
      std::cout << "b_normal_val " << b_normal_val << std::endl;

      double term_1(0.0);
      double term_2(0.0);

      term_1 = (b_normal_val * b_normal_val) * (1.0/air_permeability - 1.0/sphere_permeability);
      term_2 = (h_tangent_val * h_tangent_val) * (air_permeability - sphere_permeability);

      // Measure the area of the boundary
      area += ip.weight * face_weight;
      area_i += ip.weight * face_weight;

      // Integrate alpha * n.Grad(x) + beta * x
      double force_density_j = ((term_1 + term_2) / 2.0);
      std::cout 
        << "q.Eval() " << q.Eval(*f_tr, ip) 
        << " ip.weight " << ip.weight 
        << " face_weight " << face_weight 
        << std::endl;
      force_j = q.Eval(*f_tr, ip) * force_density_j * ip.weight * face_weight;
      force_density_i += force_density_j;
      force_i += force_j;
      force += force_j;
      force_density += force_density_j;
      // flux += q.Eval(*f_tr, ip) * b_val * ip.weight * face_weight * (1.0/air_permeability - 1.0/sphere_permeability);
      std::cout << "force_j " << force_j << " force_density_j " << force_density_j << std::endl;
    }
    std::cout << "force_i: " << force_i << ", force_density_i: " << force_density_i << ", area_i: " << area_i <<  std::endl;
  }

  std::cout << "\n\nforce: " << force << ", force_density: " << force_density << ", area: " << area <<  std::endl;

  MPI_Allreduce(&force, &total_force, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  return total_force;
  // return area;
}

double
calcMaxwellStressTensor(mfem::GridFunction * b_field, mfem::GridFunction * h_field, int face_attr)
{
  mfem::ConstantCoefficient one_coef(1.0);
  return calcMaxwellStressTensor(b_field, h_field, face_attr, one_coef);
}

MaxwellStressTensorAux::MaxwellStressTensorAux(std::string b_name, std::string h_name, int face_attr, std::string coef_name)
  : _b_name(std::move(b_name)), _h_name(std::move(h_name)), _coef_name(std::move(coef_name)), _face_attr(face_attr)
{
}

void
MaxwellStressTensorAux::Init(const hephaestus::GridFunctions & gridfunctions,
                     hephaestus::Coefficients & coefficients)
{
  _b_gf = gridfunctions.Get(_b_name);
  _h_gf = gridfunctions.Get(_h_name);
  if (coefficients._scalars.Has(_coef_name))
  {
    _coef = coefficients._scalars.Get(_coef_name);
  }
}

void
MaxwellStressTensorAux::Solve(double t)
{
  double force;
  // TODO: _gf should be _b_gf and _h_gf
  if (_coef != nullptr)
  {
    force = calcMaxwellStressTensor(_b_gf, _h_gf, _face_attr, *_coef);
  }
  else
  {
    force = calcMaxwellStressTensor(_b_gf, _h_gf, _face_attr);
  }

  _times.Append(t);
  _forces.Append(force);
}

} // namespace hephaestus
