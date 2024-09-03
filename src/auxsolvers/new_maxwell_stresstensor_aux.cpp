#include "new_maxwell_stresstensor_aux.hpp"

namespace hephaestus
{
double
calcMaxwellStressTensor(mfem::GridFunction * b_field, mfem::GridFunction * h_field, int face_attr, mfem::GridFunction & gf)
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

  mfem::FiniteElementSpace * gf_fes = gf.FESpace();
  mfem::FiniteElementSpace * b_fes = b_field->FESpace();
  mfem::Mesh * mesh = b_fes->GetMesh(); // Should be same for B and H?
  mfem::FiniteElementSpace * h_fes = h_field->FESpace();
  // mfem::Mesh * mesh = b_fes->GetMesh(); // Should be same for B and H?

  mfem::Vector b_local_dofs, h_local_dofs, normal_vec;
  mfem::DenseMatrix b_dshape, h_dshape, dshape;
  mfem::Vector el_shape;
  mfem::Array<int> b_dof_ids, h_dof_ids, g_dof_ids, b_bdr_dof_ids, gf_bdr_dof_ids;

  int elndofs;
  mfem::Array<int> v_dofs, dofs;

  int bdim = b_fes->GetVDim();
  std::cout << "bdim " << bdim 
    << " hdim " << h_fes->GetVDim() 
    << " gfdim " << gf_fes->GetVDim() 
    << std::endl;

  mfem::Vector vals;
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

    // const mfem::FiniteElement & b_elem = *b_fes->GetFE(f_tr->Elem1No);
    // const mfem::FiniteElement & h_elem = *h_fes->GetFE(f_tr->Elem1No);
    // const int int_order = 2 * b_elem.GetOrder() + 3;
    // const mfem::IntegrationRule & ir = mfem::IntRules.Get(f_tr->FaceGeom, int_order);
    gf_fes->GetElementVDofs(i, dofs);
    const mfem::FiniteElement *el = gf_fes->GetFE(i);
    elndofs = el->GetDof();
    gf_fes->GetBdrElementDofs(f_tr->Elem1No, gf_bdr_dof_ids);
    int dim = el->GetDim();
    vals.SetSize(dofs.Size());
    dshape.SetSize(elndofs, dim);
    // const mfem::FiniteElement & h_elem = *fes->GetFE(f_tr->Elem1No);
    // f_tr->Attribute = mesh->GetAttribute(f_tr->Elem1No);
    
    // // const int h_int_order = 2 * h_elem.GetOrder() + 3;
    // std::cout << "Get integration rule " << std::endl;

    b_fes->GetElementDofs(f_tr->Elem1No, b_dof_ids);
    b_fes->GetBdrElementDofs(f_tr->Elem1No, b_bdr_dof_ids);
    b_field->GetSubVector(b_dof_ids, b_local_dofs);

    const int space_dim = f_tr->Face->GetSpaceDim();
    normal_vec.SetSize(space_dim);

    // std::cout << "b_local_dofs size " << b_local_dofs.Size() 
    //   << " b_bdr_dofs size" << b_bdr_dof_ids.Size() 
    //   << " elndofs " << elndofs  
    //   << " gf_bdr_dof_ids size " <<  gf_bdr_dof_ids.Size()
    //   << std::endl;

    for (int dof = 0; dof < elndofs; ++dof)
    {
      // get normal vec
      double face_weight = f_tr->Face->Weight();
      // std::cout << "CalcOrtho" << std::endl;
      // mfem::CalcOrtho(f_tr->Face->Jacobian(), normal_vec);
      // normal_vec /= face_weight;

      // Project
      // std::cout << "Project" << std::endl;
      const mfem::IntegrationPoint &ip = el->GetNodes().IntPoint(dof);
      f_tr->SetIntPoint(&ip);

      // std::cout << "Calcshape" << std::endl;
      el->CalcShape(f_tr->GetIntPoint(), el_shape);
      // grad_hat.SetSize(vdim, dim);
      // mfem::DenseMatrix loc_data_mat(b_local_dofs.GetData(), elndofs, vdim);
    
      // std::cout << "Inner product" << std::endl;
      // double b_normal_val = mfem::InnerProduct(normal_vec, b_local_dofs);
      // std::cout << "normal vec " << normal_vec.Size()<< " b_dof " << b_local_dofs.Size() << std::endl;
      b_field->GetSubVector(b_bdr_dof_ids, b_local_dofs);
      assert (b_local_dofs.Size() == 1);
      double b_normal_val = b_local_dofs(0);

      // std::cout << "b_normal_val " << b_normal_val << std::endl;

      double term_1(0.0);
      double term_2(0.0);

      term_1 = (b_normal_val * b_normal_val) * (1.0/air_permeability - 1.0/sphere_permeability);

      // Measure the area of the boundary
      area += ip.weight * face_weight;

      vals(dof) = term_1;
    }
    // Accumulate values in all dofs
    // for (int j = 0; j < gf_bdr_dof_ids.Size(); j++)
    for (int j = 0; j < dofs.Size(); j++)
    {
        // int ldof = gf_bdr_dof_ids[j];
        int ldof = dofs[j];
        gf(ldof) += vals[j];
        std::cout << "vals[j] = " << vals[j] << " at dof " << ldof << std::endl;
    }
  }
  std::cout << "end of loop " << std::endl;

    /*
    const mfem::IntegrationRule & ir = mfem::IntRules.Get(f_tr->FaceGeom, int_order);
    g_fes->GetElementDofs(f_tr->Elem1No, g_dof_ids);
    b_fes->GetElementDofs(f_tr->Elem1No, b_dof_ids);
    b_field->GetSubVector(b_dof_ids, b_local_dofs);
    h_fes->GetElementDofs(f_tr->Elem1No, h_dof_ids);
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
      f_tr->Elem1->SetIntPoint(&eip);
      b_elem.CalcVShape(*f_tr->Elem1, b_dshape);
      h_elem.CalcVShape(*f_tr->Elem1, h_dshape);
      mfem::CalcOrtho(f_tr->Face->Jacobian(), normal_vec);
      normal_vec /= face_weight;

      double b_normal_val_j = b_dshape.InnerProduct(normal_vec, b_local_dofs);
      double b_normal_val = b_normal_val_j;

      std::cout << "b_normal_val " << b_normal_val << std::endl;

      double term_1(0.0);
      double term_2(0.0);

      term_1 = (b_normal_val * b_normal_val) * (1.0/air_permeability - 1.0/sphere_permeability);
      term_2 = 0.0;///(h_tangent_val * h_tangent_val) * (air_permeability - sphere_permeability);

      // Measure the area of the boundary
      area += ip.weight * face_weight;
      area_i += ip.weight * face_weight;

      // Integrate alpha * n.Grad(x) + beta * x
      double force_density_j = ((term_1 + term_2) / 2.0);
      // std::cout 
      //   // << "q.Eval() " << q.Eval(*f_tr, ip) 
      //   << " ip.weight " << ip.weight 
      //   << " face_weight " << face_weight 
      //   << std::endl;
      force_j = 
        // q.Eval(*f_tr, ip) * 
        force_density_j * ip.weight * face_weight;
      force_density_i += force_density_j;
      force_i += force_j;
      force += force_j;
      force_density += force_density_j;
      // flux += q.Eval(*f_tr, ip) * b_val * ip.weight * face_weight * (1.0/air_permeability - 1.0/sphere_permeability);
      std::cout << "force_j " << force_j << " force_density_j " << force_density_j << " dof " << g_dof_ids[j] << std::endl;
      gf(g_dof_ids[j]) += force_j;
    }
    std::cout << "force_i: " << force_i << ", force_density_i: " << force_density_i << ", area_i: " << area_i <<  std::endl;
  }
  */

  std::cout << "\n\nforce: " << force << ", force_density: " << force_density << ", area: " << area <<  std::endl;

  MPI_Allreduce(&force, &total_force, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  std::cout << "end of loop " << std::endl;
  return total_force;
  // return area;
}

double
calcMaxwellStressTensor(mfem::GridFunction * b_field, mfem::GridFunction * h_field, int face_attr)
{
  mfem::ConstantCoefficient one_coef(1.0);
  // return calcMaxwellStressTensor(b_field, h_field, face_attr, one_coef);
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
  std::cout << "Finding " << _coef_name << " " 
    << coefficients._scalars.Has(_coef_name) << " " 
    << gridfunctions.Has(_coef_name) << " " 
    << std::endl;
  if (gridfunctions.Has(_coef_name))
  {
    _gf = gridfunctions.Get(_coef_name);
  }
}

void
MaxwellStressTensorAux::Solve(double t)
{
  double force;
  // TODO: _gf should be _b_gf and _h_gf
  if (_gf != nullptr)
  {
    std::cout << "Passing a coef to calc" << std::endl;
    force = calcMaxwellStressTensor(_b_gf, _h_gf, _face_attr, *_gf);
  }
  else
  {
    std::cout << "Passing a dummy coef to calc" << std::endl;
    force = calcMaxwellStressTensor(_b_gf, _h_gf, _face_attr);
  }

  _times.Append(t);
  _forces.Append(force);
}

} // namespace hephaestus
