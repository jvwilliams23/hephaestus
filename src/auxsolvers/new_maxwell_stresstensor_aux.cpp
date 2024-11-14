#include "new_maxwell_stresstensor_aux.hpp"

/*
TODO:
  how to get `h_tangent_val`?

Questions:
  - Is surf force density = $terms * ip.weight * face_weight$? and force is $terms * ip.weight$ only?
*/

namespace hephaestus
{
double
calcMaxwellStressTensor(mfem::ParGridFunction * b_field, mfem::ParGridFunction * h_field, int face_attr, mfem::ParGridFunction & gf)
{
  return 0.0;
}


double
calcSurfaceForceDensity(mfem::ParGridFunction * b_field, mfem::ParGridFunction * h_field, int face_attr, mfem::ParGridFunction & gf)
{
  double flux = 0.0;
  double force = 0.0;
  double total_force = 0.0;
  double force_density = 0.0;
  double area = 0.0;

  double air_permeability = M_PI * 4.0e-7;
  double sphere_permeability = 500*air_permeability;

  mfem::ParFiniteElementSpace * gf_fes = gf.ParFESpace();
  mfem::ParFiniteElementSpace * b_fes = b_field->ParFESpace();
  mfem::ParFiniteElementSpace * h_fes = h_field->ParFESpace();

  mfem::ParMesh * mesh = gf_fes->GetParMesh();

  mfem::Vector normal_vec, unit_normal_vec;
  mfem::Array<int> g_dof_ids;

  mfem::ElementTransformation *eltrans = NULL;
  mfem::FaceElementTransformations * f_tr = NULL;
  bool use_eltrans(false);
  
  for (int i = 0; i < mesh->GetNBE(); i++)
  {
    if (mesh->GetBdrAttribute(i) != face_attr)
      continue;

    // get dofs for writing to gridfunction
    gf_fes->GetBdrElementDofs(i, g_dof_ids);

    if (use_eltrans)
    {
      eltrans = b_fes->GetBdrElementTransformation(i);
    }
    else
    {
      f_tr =
          mesh->GetFaceElementTransformations(mesh->GetBdrElementFaceIndex(i));
    }
    const mfem::FiniteElement &elem = *b_fes->GetBE(i);
    const mfem::IntegrationRule *ir = NULL;
    if (ir == NULL)
    {
      if (use_eltrans)
      {
        const int order = 2*elem.GetOrder() + eltrans->OrderW(); // <-----
        ir = &mfem::IntRules.Get(eltrans->GetGeometryType(), order);
      }
      else
      {
        const int order = 2 * elem.GetOrder() + 3;
        ir = &mfem::IntRules.Get(f_tr->FaceGeom, order);
      }
    }
    const int space_dim = 3;

    // // get coordinates for outputting angle vs force (post-proc)
    mfem::Element * be = mesh->GetBdrElement(i);
    mfem::Array<int> vertices;
    be->GetVertices(vertices);
    normal_vec.SetSize(space_dim);
    unit_normal_vec.SetSize(space_dim);

    const mfem::FiniteElement *el = gf_fes->GetFE(i);

    double force_i = 0.0;
    double area_i = 0.0;
    double force_density_i = 0.0;

    for (int j = 0; j < ir->GetNPoints(); j++)
    {
      const mfem::IntegrationPoint & ip = ir->IntPoint(j);
      double face_weight(0.0);
      mfem::IntegrationPoint eip;
      if (use_eltrans)
      {
        eltrans->SetIntPoint(&ip);
        mfem::CalcOrtho(eltrans->Jacobian(), normal_vec);
        face_weight = normal_vec.Norml2();
      }
      else
      {
        f_tr->Loc1.Transform(ip, eip);
        f_tr->Face->SetIntPoint(&ip);
        mfem::CalcOrtho(f_tr->Face->Jacobian(), normal_vec);
        face_weight = f_tr->Face->Weight();
        f_tr->Elem1->SetIntPoint(&eip);
      }
      // setup empty vectors
      mfem::Vector b_vec(space_dim);
      mfem::Vector h_vec(space_dim);
      mfem::Vector h_tang(space_dim);
      
      // get vector values at integration point
      if (use_eltrans)
      {
        b_field->GetVectorValue(*eltrans, ip, b_vec);
        h_field->GetVectorValue(*eltrans, ip, h_vec);
      }
      else
      {
        mfem::Vector loc_data;
        b_field->GetVectorValue(*f_tr->Elem1, ip, b_vec);
        h_field->GetVectorValue(*f_tr->Elem1, ip, h_vec);
        // int vdim = mfem::VectorDim();
        // mfem::DenseMatrix vshape(dof, vdim);
        // elem.CalcVShape(*f_tr->Elem1, vshape);
        // val.SetSize(vdim);
        // vshape.MultTranspose(loc_data, val);

      }

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

      double force_density_j = ((term_1 - term_2) / 2.0);
      // take y component of force for hollow sphere/levitation example
      double force_j = ((term_1 - term_2) / 2.0) * ip.weight * face_weight * unit_normal_vec(1);

      force_i += force_j;
      force_density_i += force_density_j;
    }
    force_density += force_density_i;
    // write force_density_i to all g_dofs (not sure this is correct)
    for (int j = 0; j < g_dof_ids.Size(); j++)
    {
      int ldof = g_dof_ids[j];
      gf(ldof) = force_density_i;
    }
    force += force_i;
  }


  
  // Rfs.close();
  std::cout 
    << "force: " << force 
    << ", force_density: " << force_density 
    << ", area: " << area 
    <<  std::endl;

  MPI_Allreduce(&force, &total_force, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  // std::cout << "end of loop " << std::endl;
  return total_force;
  // return area;
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

  // InitChildMesh();
  // MakeFESpaces(0);
  // MakeGridFunctions(0);
  // MakeFESpaces(1);
  // MakeGridFunctions(1);

  // _mesh_child->Transfer(*_gf, *_gf_child);
}

// ************************************************************************** //

void
MaxwellStressTensorAux::Solve(double t)
{
  double force;

  // _mesh_child->Transfer(*_b_gf, *_b_gf_child);
  // _mesh_child->Transfer(*_h_gf, *_h_gf_child);
  if (_gf != nullptr)
  {
    std::cout << "Passing a gf to calc" << std::endl;
    // force = calcSurfaceForceDensity(_b_gf_child.get(), _h_gf_child.get(), 101, *_gf_child.get());
    force = calcSurfaceForceDensity(_b_gf, _h_gf, 101, *_gf);
    
    // std::cout << "\n\nCalculating stresses on inner surface" << std::endl;
    // calcMaxwellStressTensor(_b_gf_child.get(), _h_gf_child.get(), 102, *_gf_child.get());
    // std::cout << "\n\nCalculating stresses on outer surface" << std::endl;
    // calcMaxwellStressTensor(_b_gf_child.get(), _h_gf_child.get(), 101, *_gf_child.get());
    // calcMaxwellStressTensor(_b_gf, _h_gf, 101, *_gf);
    // if (_gf)
    //   _mesh_child->Transfer(*_gf_child, *_gf);

    // WriteForces("gf_coords.txt", *_gf, _face_attr);

    std::ostringstream mesh_name, fes_name, b_fes_name, h_fes_name;
    int myid = mfem::Mpi::WorldRank();
    mesh_name << "mesh." << std::setfill('0') << std::setw(6) << myid;
    fes_name << "gf_field." << std::setfill('0') << std::setw(6) << myid;
    std::ofstream mesh_ofs(mesh_name.str().c_str());
    std::ofstream fes_ofs(fes_name.str().c_str());
    _mesh_parent->Print(mesh_ofs);
    _gf->Save(fes_ofs);
  }
  else
  {
    std::cout << "Passing a dummy coef to calc" << std::endl;
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
    double force_density_i = gf.GetVectorValue(*eltrans);
    Rfs << " " << force_density_i;
    Rfs << "\n";
  }


  
  Rfs.close();
}

// ************************************************************************** //

void
MaxwellStressTensorAux::InitChildMesh()
{
  mfem::Array<int> domain_marker;
  domain_marker.Append(100);
  if (_mesh_child == nullptr)
  {
    _mesh_child = std::make_unique<mfem::ParSubMesh>(
        mfem::ParSubMesh::CreateFromDomain(*_mesh_parent, domain_marker));
  }
}


void
MaxwellStressTensorAux::MakeFESpaces(int stage)
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
MaxwellStressTensorAux::MakeGridFunctions(int stage)
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
