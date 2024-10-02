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
  double force = 0.0;
  double total_force = 0.0;
  double force_density = 0.0;
  double area = 0.0;

  double air_permeability = M_PI * 4.0e-7; //1.25663706e-6;````
  // double air_permeability = 1.0;
  // double air_permeability = 1.25663706e-6;
  double sphere_permeability = 500*air_permeability;

  std::cout << "Get FES" << std::endl;
  mfem::FiniteElementSpace * gf_fes = gf.ParFESpace();
  std::cout << "Get FES" << std::endl;
  mfem::FiniteElementSpace * b_fes = b_field->ParFESpace();
  mfem::FiniteElementSpace * h_fes = h_field->ParFESpace();
  // h_field->ProjectGridFunction(*b_field);
  // mfem::Mesh * mesh = b_fes->GetMesh(); // Should be same for B and H?

  std::cout << "Get mesh" << std::endl;
  mfem::Mesh * mesh = gf_fes->GetMesh();
  // mfem::Mesh * mesh = b_fes->GetMesh(); // Should be same for B and H?

  mfem::Vector b_local_dofs, h_local_dofs, normal_vec;
  mfem::DenseMatrix b_dshape, h_dshape, dshape;
  mfem::Vector el_shape;
  mfem::Array<int> b_dof_ids, h_dof_ids, g_dof_ids, b_bdr_dof_ids, h_bdr_dof_ids, gf_bdr_dof_ids;

  int elndofs;
  mfem::Array<int> v_dofs, dofs;


  std::cout << "Get vdim" << std::endl;
  int bdim = b_fes->GetVDim();
  std::cout << "bdim " << bdim 
    << " hdim " << h_fes->GetVDim() 
    << " gfdim " << gf_fes->GetVDim() 
    << std::endl;

  mfem::Vector vals;
  mfem::Vector ones;

  double max_flux (-1.0);
  double min_flux (1.0);
  // for post-proc
  mfem::GridFunction gf_coords(gf_fes);
  mfem::DenseMatrix coords, coords_t;
  mfem::Array<int> gf_vdofs;
  std::ofstream Rfs("gf_coords.txt", std::ofstream::out);
  // for (int i = 0; i < mesh->GetNE(); i++)
  // {
  //   // get coordinates 
  //   const mfem::FiniteElement & gf_elem = *gf_fes->GetFE(i);
  //   const mfem::IntegrationRule &gf_nodes = gf_elem.GetNodes();
  //   mfem::ElementTransformation *T = mesh->GetElementTransformation(i);
  //   // f_tr->Transform(gf_nodes, coords);
  //   T->Transform(gf_nodes, coords);
  //   coords_t.Transpose(coords);
  //   gf_fes->GetElementVDofs(i, gf_vdofs);
  //   gf_coords.SetSubVector(gf_vdofs, coords_t.GetData());


  //   // std::cout << "dim " << coords_t.Size() << "coord is " << coords_t.GetData()[0] << " " << coords_t.GetData()[1] << " " << coords_t.GetData()[2] << std::endl;
  //   Rfs 
  //     << coords_t.GetData()[0] << " "
  //     << coords_t.GetData()[1] << " "
  //     << coords_t.GetData()[2] << " "
  //   ;
  //   Rfs << "\n";
  // }
  mfem::ElementTransformation *eltrans = NULL;
  
  std::cout << "Looping GetNBE " << std::endl;
  for (int i = 0; i < mesh->GetNBE(); i++)
  {
    if (mesh->GetBdrAttribute(i) != face_attr)
      continue;

    // mfem::FaceElementTransformations * f_tr =
    //     mesh->GetFaceElementTransformations(mesh->GetBdrElementFaceIndex(i));

    eltrans = b_fes->GetBdrElementTransformation(i);
    const mfem::FiniteElement &gf_elem = *b_fes->GetBE(i);

    // if (f_tr == nullptr)
    //   continue;

    // std::cout << "f_tr->OrderW() = " << f_tr->OrderW() << std::endl;
    // const mfem::FiniteElement & gf_elem = *gf_fes->GetFE(f_tr->Elem1No);
    // const mfem::FiniteElement & h_elem = *h_fes->GetFE(f_tr->Elem1No);
    // const mfem::FiniteElement & b_elem = *b_fes->GetFE(f_tr->Elem1No);
    const mfem::IntegrationRule *ir = NULL;
    if (ir == NULL)
    {
        const int order = 2*gf_elem.GetOrder() + eltrans->OrderW(); // <-----
        ir = &mfem::IntRules.Get(eltrans->GetGeometryType(), order);
    }
    // const int int_order = 2 * gf_elem.GetOrder() + f_tr->OrderW();
    // const mfem::IntegrationRule & ir = mfem::IntRules.Get(f_tr->FaceGeom, int_order);
    // gf_fes->GetElementDofs(f_tr->Elem1No, g_dof_ids);
    // b_fes->GetElementDofs(f_tr->Elem1No, b_dof_ids);
    // b_field->GetSubVector(b_dof_ids, b_local_dofs);
    // h_fes->GetElementDofs(f_tr->Elem1No, h_dof_ids);
    // h_field->GetSubVector(h_dof_ids, h_local_dofs);
    const int space_dim = 3;//f_tr->Face->GetSpaceDim();


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

    // std::cout << "f_tr " << f_tr->ElementNo << " " << f_tr->Elem1No << " " << i << std::endl;
    // b_fes->GetBdrElementDofs(i, b_bdr_dof_ids);
    // b_field->GetSubVector(b_bdr_dof_ids, b_local_dofs);
    // assert (b_local_dofs.Size() == 1);
    // h_fes->GetBdrElementDofs(i, h_bdr_dof_ids);
    // h_field->GetSubVector(h_bdr_dof_ids, h_local_dofs);
    
    // b_dshape.SetSize(b_elem.GetDof(), space_dim);
    // h_dshape.SetSize(h_elem.GetDof(), space_dim);
    normal_vec.SetSize(space_dim);

    const mfem::FiniteElement *el = gf_fes->GetFE(i);

    double force_i = 0.0;
    double area_i = 0.0;
    double force_density_i = 0.0;

    /**/
    for (int j = 0; j < ir->GetNPoints(); j++)
    {
      const mfem::IntegrationPoint & ip = ir->IntPoint(j);
      // mfem::IntegrationPoint eip;
      // f_tr->Loc1.Transform(ip, eip);
      // f_tr->Face->SetIntPoint(&ip);
      eltrans->SetIntPoint(&ip);
      // double face_weight = f_tr->Face->Weight();
      // double force_j(0.0);
      // std::cout << "\n\nLoop over NPoints " << j << " out of " << ir.GetNPoints() << std::endl;
      // const mfem::IntegrationPoint &ip = el->GetNodes().IntPoint(dof);
      // // const mfem::IntegrationPoint & ip = ir.IntPoint(j);
      // mfem::IntegrationPoint eip;
      // f_tr->Loc1.Transform(ip, eip);
      // f_tr->Face->SetIntPoint(&ip);
      // double face_weight = f_tr->Face->Weight();
      // f_tr->Elem1->SetIntPoint(&eip);
      
      // compute b normal component
      mfem::CalcOrtho(eltrans->Jacobian(), normal_vec);
      double face_weight = normal_vec.Norml2();
      mfem::Vector b_vec(space_dim);
      mfem::Vector h_vec(space_dim);
      // double b_normal_val = b_dshape.InnerProduct(normal_vec, b_local_dofs) / face_weight; /// face_weight;
      b_field->GetVectorValue(*eltrans, ip, b_vec);
      double b_normal_val = b_vec * normal_vec / face_weight;


      std::cout << "B-vector (" 
        << b_vec(0) <<  " " 
        << b_vec(1) <<  " " 
        << b_vec(2) <<  ") " 
        << " normal val "
        << b_normal_val << " "
        << " normal vec ("
        << normal_vec(0) <<  " " 
        << normal_vec(1) <<  " " 
        << normal_vec(2) <<  ") " 
        << std::endl;

      // double h_normal_val = h_dshape.InnerProduct(normal_vec, h_local_dofs) / face_weight;
      // mfem::Vector h_vec(space_dim);
      /*mfem::Vector h_tang(space_dim);
      h_field->GetVectorValue(*eltrans, ip, h_vec);
      double h_normal_val = h_vec * normal_vec / face_weight;
      for (int k = 0; k < space_dim; ++k){
        h_tang(k) = h_vec(k) - (normal_vec(k)*h_normal_val);
      }
      double h_tangent_val = h_tang.Norml2();
      std::cout << "H-vector " 
        << h_vec(0) <<  " " 
        << h_vec(1) <<  " " 
        << h_vec(2) <<  " " 
        << " normal val "
        << h_normal_val << " "
        << " tangent val "
        << h_tangent_val << " "
        << std::endl;
      */
      if (b_normal_val > max_flux) { max_flux = b_normal_val; }
      if (b_normal_val < min_flux) { min_flux = b_normal_val; }


      // h_normal_val += h_dshape.InnerProduct(normal_vec, h_local_dofs) / face_weight;
      // h_tangent_val += h_dshape.InnerProduct(normal_vec, h_local_dofs) / face_weight;


      // std::cout << "Jacobian size " << f_tr->Face->Jacobian().NumCols() <<  " x " << f_tr->Face->Jacobian().NumRows() << std::endl;
      // std::cout << "Jacobian vals " 
      //   << f_tr->Face->Jacobian()(0,0) <<  " " 
      //   << f_tr->Face->Jacobian()(1,0) <<  " " 
      //   << f_tr->Face->Jacobian()(2,0) <<  " " 
      //   << f_tr->Face->Jacobian()(0,1) <<  " " 
      //   << f_tr->Face->Jacobian()(1,1) <<  " " 
      //   << f_tr->Face->Jacobian()(2,1) <<  " " 
      //   << std::endl;
      // for (int counter = 0; counter < h_local_dofs.Size(); ++counter)
      //   std::cout << counter << " h_dof " << h_local_dofs(counter) << std::endl;

      double term_1(0.0);
      double term_2(0.0);
      std::cout << "b_normal_val " << b_normal_val << std::endl;

      term_1 = (b_normal_val * b_normal_val) * (1.0/air_permeability - 1.0/sphere_permeability);// / (face_weight*face_weight*ip.weight*ip.weight);
      // term_2 = (h_tangent_val * h_tangent_val) * (air_permeability - sphere_permeability);

      // Measure the area of the boundary
      // area += face_weight;
      area += ip.weight * face_weight;

      // Integrate alpha * n.Grad(x) + beta * x
      double force_density_j = ((term_1 - term_2) / 2.0) * ip.weight ;
      double force_j = ((term_1 - term_2) / 2.0) * ip.weight  ;
      std::cout 
        << "force density " << force_density_j
        << " ip.weight " << ip.weight 
        << " face_weight " << face_weight 
        << std::endl;
      force_i += force_j;
      force_density_i += force_density_j;
    }
    force_density += force_density_i;
    // write force_density_i to all g_dofs (not sure this is correct)
    for (int j = 0; j < g_dof_ids.Size(); j++)
    {
      int ldof = g_dof_ids[j];
      // gf(ldof) += force_density_i;
      gf(ldof) = y_coord;
    }
    force += force_i;
    Rfs << " " << force_i << " " << force_density_i;
    /*
    vals.SetSize(g_dof_ids.Size());
    for (int dof = 0; dof < g_dof_ids.Size(); ++dof)
    {
      double force_j(0.0);
      // std::cout << "\n\nLoop over NPoints " << j << " out of " << ir.GetNPoints() << std::endl;
      const mfem::IntegrationPoint &ip = el->GetNodes().IntPoint(dof);
      // const mfem::IntegrationPoint & ip = ir.IntPoint(j);
      mfem::IntegrationPoint eip;
      f_tr->Loc1.Transform(ip, eip);
      f_tr->Face->SetIntPoint(&ip);
      double face_weight = f_tr->Face->Weight();
      f_tr->Elem1->SetIntPoint(&eip);
      
      double b_normal_val = b_local_dofs(0);

      if (b_normal_val > max_flux) { max_flux = b_normal_val; }
      if (b_normal_val < min_flux) { min_flux = b_normal_val; }


      // std::cout << "h_local_dofs size " << h_local_dofs.Size() << std::endl;
      // for (int counter = 0; counter < h_local_dofs.Size(); ++counter)
      //   std::cout << counter << " h_dof " << h_local_dofs(counter) << std::endl;

      double term_1(0.0);
      double term_2(0.0);
      std::cout << "b_normal_val " << b_normal_val << std::endl;

      term_1 = (b_normal_val * b_normal_val) * (1.0/air_permeability - 1.0/sphere_permeability) / (face_weight*face_weight);
      // term_2 = (h_tangent_val * h_tangent_val) * (air_permeability - sphere_permeability);

      // Measure the area of the boundary
      area += face_weight;
      // area_i += ip.weight * face_weight;

      // Integrate alpha * n.Grad(x) + beta * x
      double force_density_j = ((term_1 - term_2) / 2.0);
      std::cout 
        << "force density " << force_density_j
        << " ip.weight " << ip.weight 
        << " eip.weight " << eip.weight 
        << " face_weight " << face_weight 
        << std::endl;
      // FIXME: Why ip.weight == 0? ignore for now
      force_j = force_density_j * face_weight ;
      //   // q.Eval(*f_tr, ip) * 
      // force_density_i += force_density_j;
      // force_i += force_j;
      // force += force_j;
      // force_density += force_density_j;
      // // flux += q.Eval(*f_tr, ip) * b_val * ip.weight * face_weight * (1.0/air_permeability - 1.0/sphere_permeability);
      // std::cout << "force_j " << force_j << " force_density_j " << force_density_j << " dof " << g_dof_ids[j] << std::endl;
      // gf(g_dof_ids[j]) += force_j;


      vals(dof) = force_density_j;
      // el.CalcShape(ip, el_shape);
      // gf.add(force_j, el_shape);
      force_density += force_density_j;
      force += force_j;
    }
    for (int j = 0; j < g_dof_ids.Size(); j++)
    {
      int ldof = g_dof_ids[j];
      gf(ldof) += vals[j];
      if (x_coord != 0.0)
      {
        Rfs << " " << gf(ldof);
      }
    }
    */
    Rfs << "\n";
    // std::cout << "force_i: " << force_i << ", force_density_i: " << force_density_i << ", area_i: " << area_i <<  std::endl;
  }


  
  Rfs.close();
  std::cout 
    << "\n\nforce: " << force 
    << ", force_density: " << force_density 
    << ", area: " << area 
    << " (expected area is 0.0314) "
    << " min flux " << min_flux 
    << " max flux " << max_flux 
    <<  std::endl;

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

  _mesh_parent = _b_gf->ParFESpace()->GetParMesh();

  std::cout << "Finding " << _coef_name << " " 
    << coefficients._scalars.Has(_coef_name) << " " 
    << gridfunctions.Has(_coef_name) << " " 
    << std::endl;
  if (gridfunctions.Has(_coef_name))
  {
    _gf = gridfunctions.Get(_coef_name);
  }

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
  if (_gf != nullptr)
  {
    std::cout << "Passing a gf to calc" << std::endl;
    // force = calcMaxwellStressTensor(_b_gf, _h_gf, _face_attr, *_gf);
    force = calcMaxwellStressTensor(_b_gf_child.get(), _h_gf_child.get(), _face_attr, *_gf_child.get());
  }
  else
  {
    std::cout << "Passing a dummy coef to calc" << std::endl;
    // force = calcMaxwellStressTensor(_b_gf, _h_gf, _face_attr);
  }

  _times.Append(t);
  _forces.Append(force);
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
