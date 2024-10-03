#include "new_maxwell_stresstensor_aux.hpp"

/*
TODO:
  how to get `h_tangent_val`?

Questions:
  - what about inner sphere surface? I not really considered this
  - Is surf force density = $terms * ip.weight * face_weight$? and force is $terms * ip.weight$ only?
*/

namespace hephaestus
{
double
calcMaxwellStressTensor(mfem::ParGridFunction * b_field, mfem::ParGridFunction * h_field, int face_attr, mfem::ParGridFunction & gf)
{
  double flux = 0.0;
  double force_mag_mst = 0.0;
  double force_mag_jm = 0.0;
  double force_mag_qm = 0.0;
  double total_force = 0.0;
  double area = 0.0;

  double air_permeability = M_PI * 4.0e-7;
  double mu0_ = air_permeability;
  double muR_ = 500.0;
  double muRSqr_ = muR_ * muR_;
  double sphere_permeability = 500*air_permeability;

  mfem::ParFiniteElementSpace * gf_fes = gf.ParFESpace();
  mfem::ParFiniteElementSpace * b_fes = b_field->ParFESpace();
  mfem::ParFiniteElementSpace * h_fes = h_field->ParFESpace();

  mfem::ParFiniteElementSpace * fes_to_use = b_field->ParFESpace();
  // mfem::ParFiniteElementSpace * fes_to_use = gf.ParFESpace();

  mfem::ParMesh * mesh = b_fes->GetParMesh();

  mfem::Vector normal_vec;
  mfem::Vector unit_normal_vec;
  mfem::Vector tangent_vec;

  const int space_dim = 3;
  normal_vec.SetSize(space_dim);
  tangent_vec.SetSize(space_dim);
  unit_normal_vec.SetSize(space_dim);

  mfem::Array<int> g_dof_ids;

  mfem::ElementTransformation *eltrans = NULL;

  int logvar(0);
  
  for (int i = 0; i < mesh->GetNBE(); i++)
  {
    if (mesh->GetBdrAttribute(i) != face_attr)
      continue;

    eltrans = fes_to_use->GetBdrElementTransformation(i);
    const mfem::FiniteElement &gf_elem = *fes_to_use->GetBE(i);
    const mfem::IntegrationRule *ir = NULL;
    if (ir == NULL)
    {
        const int order = 2*gf_elem.GetOrder() + eltrans->OrderW(); // <-----
        ir = &mfem::IntRules.Get(eltrans->GetGeometryType(), order);
    }

    for (int j = 0; j < ir->GetNPoints(); j++)
    {
      const mfem::IntegrationPoint & ip = ir->IntPoint(j);
      eltrans->SetIntPoint(&ip);
      
      // compute b normal component
      mfem::CalcOrtho(eltrans->Jacobian(), normal_vec);
      double face_weight = normal_vec.Norml2();
      unit_normal_vec.Set(1.0/face_weight, normal_vec);
      
      mfem::Vector b_vec(space_dim);
      mfem::Vector h_vec(space_dim);
      b_field->GetVectorValue(*eltrans, ip, b_vec);
      double b_normal_val = b_vec * unit_normal_vec;

      mfem::Vector h_tang(space_dim);
      h_field->GetVectorValue(*eltrans, ip, h_vec);
      double h_normal_val = h_vec * unit_normal_vec;
      // divide h_normal_val by face_weight as it is multiplied by normal vec 
      for (int k = 0; k < space_dim; ++k){
        h_tang(k) = h_vec(k) - (unit_normal_vec(k)*h_normal_val);
      }
      double h_tangent_val = h_tang.Norml2();

      mfem::Vector term_1(space_dim);
      mfem::Vector term_2(space_dim);
      mfem::Vector f_mst(space_dim);
      mfem::Vector f_jm(space_dim);
      mfem::Vector f_qm(space_dim);
      term_1 = 0.0;
      term_2 = 0.0;
      f_mst = 0.0;
      f_jm = 0.0;
      f_qm = 0.0;
      
      // maxwell stress tensor method (eq 5)
      double mag_term_1 = ip.weight*b_normal_val;
      double mag_term_2 = 0.5*ip.weight*b_normal_val * b_normal_val / mu0_;
      double mag_term_3 = -0.5*ip.weight * mu0_ * (h_tang * h_tang);\
      f_mst.Set(mag_term_1, h_tang);
      f_mst.Add(mag_term_2, unit_normal_vec);
      f_mst.Add(mag_term_3, unit_normal_vec);
      /**/
      /*
      force_vec.Set(ip.weight * b_normal_val * b_normal_val / air_permeability, normal_vec);
      force_vec.Add(ip.weight * face_weight * b_normal_val, h_tang);
      force_vec.Add(-0.5 * ip.weight * (air_permeability * (h_tang * h_tang) + b_normal_val * b_normal_val / air_permeability), normal_vec);
      */
      // magnetising current method (eq 6)
      f_jm.Set(ip.weight * b_normal_val * (1.0-muR_), h_tang);
      f_jm.Add(ip.weight * (0.5*mu0_) * (muRSqr_-1.0) * b_normal_val * b_normal_val, unit_normal_vec);
      // magnetic charge method (eq 7)
      f_qm.Set(ip.weight * face_weight * b_normal_val * (1.0-1.0/muR_), h_tang);
      f_qm.Add(ip.weight * (0.5/mu0_) * b_normal_val * b_normal_val * (1.0-1.0/muRSqr_), normal_vec);

      force_mag_mst += f_mst.Norml2();
      force_mag_jm += f_jm.Norml2();
      force_mag_qm += f_qm.Norml2();
    }
  }


  
  // Rfs.close();
  std::cout 
    << "\n\ntotal forces - f_mst: " << force_mag_mst
    << ", f_jm: " << force_mag_jm
    << ", f_qm: " << force_mag_qm 
    <<  std::endl;

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

  double air_permeability = M_PI * 4.0e-7; //1.25663706e-6;````
  double sphere_permeability = 500*air_permeability;

  // std::cout << "Get FES" << std::endl;
  mfem::ParFiniteElementSpace * gf_fes = gf.ParFESpace();
  // std::cout << "Get FES" << std::endl;
  mfem::ParFiniteElementSpace * b_fes = b_field->ParFESpace();
  mfem::ParFiniteElementSpace * h_fes = h_field->ParFESpace();
  // h_field->ProjectGridFunction(*b_field);

  mfem::ParMesh * mesh = gf_fes->GetParMesh();

  mfem::Vector normal_vec;
  mfem::Array<int> g_dof_ids;

  double max_flux (-1.0);
  double min_flux (1.0);
  // for post-proc
  // std::ofstream Rfs("gf_coords.txt", std::ofstream::out);

  mfem::ElementTransformation *eltrans = NULL;
  
  // std::cout << "Looping GetNBE " << std::endl;
  for (int i = 0; i < mesh->GetNBE(); i++)
  {
    if (mesh->GetBdrAttribute(i) != face_attr)
      continue;

    // get dofs for writing to gridfunction
    gf_fes->GetBdrElementDofs(i, g_dof_ids);

    eltrans = b_fes->GetBdrElementTransformation(i);
    const mfem::FiniteElement &gf_elem = *b_fes->GetBE(i);
    const mfem::IntegrationRule *ir = NULL;
    if (ir == NULL)
    {
        const int order = 2*gf_elem.GetOrder() + eltrans->OrderW(); // <-----
        ir = &mfem::IntRules.Get(eltrans->GetGeometryType(), order);
    }
    const int space_dim = 3;

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
    // if (x_coord != 0.0){
    //   Rfs 
    //     << x_coord << " "
    //     << y_coord << " "
    //     << z_coord << " "
    //     << rad << " "
    //     << theta_coord << " "
    //     << phi_coord //<< " "
    //   ;
    // }

    normal_vec.SetSize(space_dim);

    const mfem::FiniteElement *el = gf_fes->GetFE(i);

    double force_i = 0.0;
    double area_i = 0.0;
    double force_density_i = 0.0;

    /**/
    for (int j = 0; j < ir->GetNPoints(); j++)
    {
      const mfem::IntegrationPoint & ip = ir->IntPoint(j);
      eltrans->SetIntPoint(&ip);
      
      // compute b normal component
      mfem::CalcOrtho(eltrans->Jacobian(), normal_vec);
      double face_weight = normal_vec.Norml2();
      mfem::Vector b_vec(space_dim);
      mfem::Vector h_vec(space_dim);
      b_field->GetVectorValue(*eltrans, ip, b_vec);
      double b_normal_val = b_vec * normal_vec / face_weight;


      /* 
      std::cout << "B-vector (" 
        << b_vec(0) <<  " " 
        << b_vec(1) <<  " " 
        << b_vec(2) <<  ") " 
        << " normal val "
        << b_normal_val << " "
        << " normal vec ("
        << normal_vec(0)/face_weight <<  " " 
        << normal_vec(1)/face_weight <<  " " 
        << normal_vec(2)/face_weight <<  ") " 
        << std::endl;
      */
      mfem::Vector h_tang(space_dim);
      h_field->GetVectorValue(*eltrans, ip, h_vec);
      double h_normal_val = h_vec * normal_vec / face_weight;
      for (int k = 0; k < space_dim; ++k){
        h_tang(k) = h_vec(k) - (normal_vec(k)*h_normal_val);
      }
      double h_tangent_val = h_tang.Norml2();
      /*
      std::cout << "H-vector (" 
        << h_vec(0) <<  " " 
        << h_vec(1) <<  " " 
        << h_vec(2) <<  ") " 
        << " normal val "
        << h_normal_val << " "
        << " tangent val "
        << h_tangent_val << " "
        << std::endl;
      */
      if (b_normal_val > max_flux) { max_flux = b_normal_val; }
      if (b_normal_val < min_flux) { min_flux = b_normal_val; }

      double term_1(0.0);
      double term_2(0.0);

      term_1 = (b_normal_val * b_normal_val) * (1.0/air_permeability - 1.0/sphere_permeability);// / (face_weight*face_weight*ip.weight*ip.weight);
      term_2 = (h_tangent_val * h_tangent_val) * (air_permeability - sphere_permeability);

      // Measure the area of the boundary
      area += ip.weight * face_weight;

      // Integrate alpha * n.Grad(x) + beta * x
      double force_density_j = ((term_1 - term_2) / 2.0) * ip.weight ;
      double force_j = ((term_1 - term_2) / 2.0) * ip.weight * face_weight ;
      // std::cout 
      //   << "force density " << force_density_j
      //   << " ip.weight " << ip.weight 
      //   << " face_weight " << face_weight 
      //   << std::endl;
      force_i += force_j;
      force_density_i += force_density_j;
    }
    force_density += force_density_i;
    // write force_density_i to all g_dofs (not sure this is correct)
    for (int j = 0; j < g_dof_ids.Size(); j++)
    {
      int ldof = g_dof_ids[j];
      // std::cout << "adding to dof=" << ldof 
      //   << " (existing val is) " << gf(ldof)
      //   << " adding " 
      //   << force_density_i
      //   << std::endl;
      // gf(ldof) += force_density_i;
      gf(ldof) = force_density_i;
      // gf(ldof) = y_coord;
    }
    force += force_i;
    // Rfs << " " << force_i << " " << force_density_i;
    // Rfs << "\n";
  }


  
  // Rfs.close();
  std::cout 
    << "force: " << force 
    << ", force_density: " << force_density 
    << ", area: " << area 
    << " (expected area is 0.0314) "
    << " min flux " << min_flux 
    << " max flux " << max_flux 
    <<  std::endl;

  MPI_Allreduce(&force, &total_force, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  // std::cout << "end of loop " << std::endl;
  return total_force;
  // return area;
}

double
calcSurfaceForceDensity(mfem::GridFunction * b_field, mfem::GridFunction * h_field, int face_attr)
{
  mfem::ConstantCoefficient one_coef(1.0);
  // return calcSurfaceForceDensity(b_field, h_field, face_attr, one_coef);
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

  InitChildMesh();
  MakeFESpaces(0);
  MakeGridFunctions(0);
  MakeFESpaces(1);
  MakeGridFunctions(1);
}

void
MaxwellStressTensorAux::Solve(double t)
{
  double force;
  // TODO: _gf should be _b_gf and _h_gf

  _mesh_child->Transfer(*_b_gf, *_b_gf_child);
  _mesh_child->Transfer(*_h_gf, *_h_gf_child);
  if (_gf != nullptr)
  {
    std::cout << "Passing a gf to calc" << std::endl;
    // force = calcSurfaceForceDensity(_b_gf, _h_gf, _face_attr, *_gf);
    force = calcSurfaceForceDensity(_b_gf_child.get(), _h_gf_child.get(), 101, *_gf_child.get());
    
    std::cout << "\n\nCalculating stresses on inner surface" << std::endl;
    calcMaxwellStressTensor(_b_gf_child.get(), _h_gf_child.get(), 102, *_gf_child.get());
    std::cout << "\n\nCalculating stresses on outer surface" << std::endl;
    calcMaxwellStressTensor(_b_gf_child.get(), _h_gf_child.get(), 101, *_gf_child.get());
    // force = calcSurfaceForceDensity(_b_gf_child.get(), _h_gf_child.get(), 102, *_gf_child.get());
    if (_gf)
      _mesh_child->Transfer(*_gf_child, *_gf);

    WriteForces("gf_coords.txt", *_gf_child.get(), _face_attr);

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
  else
  {
    std::cout << "Passing a dummy coef to calc" << std::endl;
    // force = calcSurfaceForceDensity(_b_gf, _h_gf, _face_attr);
  }

  _times.Append(t);
  _forces.Append(force);
}

void 
MaxwellStressTensorAux::WriteForces(std::string fname, mfem::ParGridFunction & gf, int face_attr)
{
  // for post-proc
  mfem::ParFiniteElementSpace * gf_fes = gf.ParFESpace();
  mfem::ParMesh * mesh = gf_fes->GetParMesh();

  mfem::Array<int> g_dof_ids;
  std::ofstream Rfs(fname, std::ofstream::out);


  mfem::ElementTransformation *eltrans = NULL;
  
  
  for (int i = 0; i < mesh->GetNBE(); i++)
  {
    if (mesh->GetBdrAttribute(i) != face_attr)
      continue;

    // get dofs for writing to gridfunction
    gf_fes->GetBdrElementDofs(i, g_dof_ids);
    eltrans = gf_fes->GetBdrElementTransformation(i);
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
    // double force_density_i = 0.0;
    double force_density_i = gf.GetValue(*eltrans);
    // write force_density_i to all g_dofs (not sure this is correct)
    // for (int j = 0; j < g_dof_ids.Size(); j++)
    // {
    //   int ldof = g_dof_ids[j];
    //   force_density_i = gf(ldof);
    // }
    // force_density_i /= g_dof_ids.Size();
    Rfs << " " << force_density_i;
    Rfs << "\n";
  }


  
  Rfs.close();
}


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
