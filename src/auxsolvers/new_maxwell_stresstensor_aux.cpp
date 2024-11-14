#include "new_maxwell_stresstensor_aux.hpp"

namespace hephaestus
{
double
calcMaxwellStressTensor(mfem::ParGridFunction * b_field, mfem::ParGridFunction * h_field, int face_attr, mfem::ParGridFunction & gf)
{
  return 0.0;
}

double
calcSurfaceForceDensity(mfem::ParGridFunction * b_field, mfem::ParGridFunction * h_field, int face_attr, mfem::ParGridFunction & gf, mfem::Coefficient & mu)
{
  double flux = 0.0;
  double force = 0.0;
  double total_force = 0.0;
  double force_density = 0.0;
  double area = 0.0;

  double air_permeability = M_PI * 4.0e-7;
  // double sphere_permeability = 500*air_permeability;

  mfem::ParFiniteElementSpace * gf_fes = gf.ParFESpace();
  mfem::ParFiniteElementSpace * b_fes = b_field->ParFESpace();
  mfem::ParFiniteElementSpace * h_fes = h_field->ParFESpace();

  mfem::ParMesh * mesh = gf_fes->GetParMesh();

  mfem::Vector normal_vec, unit_normal_vec;
  mfem::Array<int> g_dof_ids;

  mfem::FaceElementTransformations * f_tr = NULL;
  
  for (int i = 0; i < mesh->GetNBE(); i++)
  {
    if (mesh->GetBdrAttribute(i) != face_attr)
      continue;

    f_tr = mesh->GetFaceElementTransformations(
      mesh->GetBdrElementFaceIndex(i)
    );
    const mfem::FiniteElement &elem = *b_fes->GetBE(i);
    const mfem::IntegrationRule *ir = NULL;
    if (ir == NULL)
    {
      const int order = 2 * elem.GetOrder() + 3;
      ir = &mfem::IntRules.Get(f_tr->FaceGeom, order);
    }
    const int space_dim = 3;

    normal_vec.SetSize(space_dim);
    unit_normal_vec.SetSize(space_dim);

    const mfem::FiniteElement *el = gf_fes->GetFE(i);

    double force_i = 0.0;
    double area_i = 0.0;
    double force_density_i = 0.0;
    double force_density_j = 0.0;

    mfem::Vector force_vec_elem(space_dim);
    force_vec_elem *= 0.0;

    // get dofs for writing to gridfunction
    // gf_fes->GetBdrElementVDofs(i, g_dof_ids);
    gf_fes->GetElementVDofs(f_tr->Elem1No, g_dof_ids);

    for (int j = 0; j < ir->GetNPoints(); j++)
    {
      mfem::Vector force_vec_ip(space_dim);
      const mfem::IntegrationPoint & ip = ir->IntPoint(j);
      double face_weight(0.0);
      mfem::IntegrationPoint eip;
      f_tr->Loc1.Transform(ip, eip);
      f_tr->Face->SetIntPoint(&ip);
      mfem::CalcOrtho(f_tr->Face->Jacobian(), normal_vec);
      face_weight = f_tr->Face->Weight();
      f_tr->Elem1->SetIntPoint(&eip);

      // setup empty vectors
      mfem::Vector b_vec(space_dim);
      mfem::Vector h_vec(space_dim);
      mfem::Vector h_tang(space_dim);
      
      // get vector values at integration point
      b_field->GetVectorValue(*f_tr->Elem1, ip, b_vec);
      h_field->GetVectorValue(*f_tr->Elem1, ip, h_vec);
      double sphere_permeability = mu.Eval(*f_tr->Elem1, ip);

      // compute b normal component
      unit_normal_vec.Set(1.0/face_weight, normal_vec);
      double b_normal_val = b_vec * unit_normal_vec;
      double h_normal_val = h_vec * unit_normal_vec;
      for (int k = 0; k < space_dim; ++k){
        h_tang(k) = h_vec(k) - (unit_normal_vec(k)*h_normal_val);
      }

      double term_1(0.0);
      double term_2(0.0);
      term_1 = (b_normal_val * b_normal_val) * (1.0/air_permeability - 1.0/sphere_permeability);
      term_2 = (h_tang * h_tang) * (air_permeability - sphere_permeability);

      // Measure the area of the boundary
      area += ip.weight * face_weight;

      force_density_j = ((term_1 - term_2) / 2.0);
      // take y component of force for hollow sphere/levitation example
      force_vec_ip.Set(((term_1 - term_2) / 2.0) * ip.weight * face_weight, unit_normal_vec);
      double force_j = force_vec_ip(1);
      force_vec_elem += force_vec_ip;

      force_i += force_j;
      force_density_i += force_density_j;
    }
    force_density += force_density_i;
    // write force, works for Vector_L2. Does not work with Vector_H1.
    for (int j = 0; j < g_dof_ids.Size(); j++)
    {
      int ldof = g_dof_ids[j];
      gf(ldof) = force_vec_elem(j);
    }
    force += force_i;
  }


  MPI_Allreduce(&force, &total_force, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  std::cout 
    << "force: " << total_force
    << ", force_density: " << force_density 
    << ", area: " << area 
    <<  std::endl;

  return total_force;
}

// ************************************************************************** //

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
  _mu_coef = coefficients._scalars.Get("magnetic_permeability");

  _mesh_parent = _b_gf->ParFESpace()->GetParMesh();

  std::cout << "Finding " << _coef_name << " " 
    << coefficients._scalars.Has(_coef_name) << " " 
    << gridfunctions.Has(_coef_name) << " " 
    << std::endl;
  if (gridfunctions.Has(_coef_name))
  {
    _gf = gridfunctions.Get(_coef_name);
  }
  // init field
  *_gf = 0.0;
}

// ************************************************************************** //

void
MaxwellStressTensorAux::Solve(double t)
{
  double force;

  if (_gf != nullptr)
  {
    std::cout << "Passing a gf to calc" << std::endl;
    force = calcSurfaceForceDensity(_b_gf, _h_gf, 101, *_gf, *_mu_coef);
    // calcMaxwellStressTensor(_b_gf, _h_gf, 101, *_gf);

    // WriteForces("gf_coords.txt", *_gf, _face_attr);
  }

  _times.Append(t);
  _forces.Append(force);
}

// ************************************************************************** //

void 
MaxwellStressTensorAux::WriteForces(std::string fname, mfem::ParGridFunction & gf, int face_attr)
{
  // for post-proc
  mfem::ParFiniteElementSpace * gf_fes = gf.ParFESpace();
  mfem::ParMesh * mesh = gf_fes->GetParMesh();

  mfem::Array<int> g_dof_ids;
  std::ofstream Rfs(fname, std::ofstream::out);

  mfem::ElementTransformation *eltrans = NULL;
  mfem::FaceElementTransformations * f_tr = NULL;  
  
  for (int i = 0; i < mesh->GetNBE(); i++)
  {
    if (mesh->GetBdrAttribute(i) != face_attr)
      continue;

    f_tr =
        mesh->GetFaceElementTransformations(mesh->GetBdrElementFaceIndex(i));

    // get dofs for writing to gridfunction
    gf_fes->GetBdrElementDofs(i, g_dof_ids);
    // eltrans = gf_fes->GetBdrElementTransformation(i);
    // // get coordinates for outputting angle vs force (post-proc)
    mfem::Element * be = mesh->GetBdrElement(i);
    mfem::Array<int> vertices;
    be->GetVertices(vertices);
    mfem::real_t * coords1 = mesh->GetVertex(vertices[0]);
    double x_coord = coords1[0];
    double y_coord = coords1[1];
    double z_coord = coords1[2];
    double rad = std::sqrt(x_coord*x_coord + y_coord*y_coord + z_coord*z_coord);
    double theta_coord = std::atan(y_coord/x_coord);
    double phi_coord = std::acos(y_coord / rad);
    // write cartesian and radial coords to file
    if (x_coord != 0.0){
      Rfs 
        << x_coord << " "
        << y_coord << " "
        << z_coord << " "
        << rad << " "
        << theta_coord << " "
        << phi_coord //<< " "
      ;
    }
    Rfs << "\n";
  }

  Rfs.close();
}

// ************************************************************************** //

} // namespace hephaestus
