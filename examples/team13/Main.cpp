#include "hephaestus.hpp"

const char * DATA_DIR = "../../data/";

hephaestus::Coefficients
defineCoefficients()
{
  hephaestus::Subdomain air("air", 11);
  air._scalar_coefficients.Register("electrical_conductivity",
                                    std::make_shared<mfem::ConstantCoefficient>(1.0));
  air._scalar_coefficients.Register("magnetic_permeability",
                                    std::make_shared<mfem::ConstantCoefficient>(M_PI * 4.0e-7));

  hephaestus::Subdomain plate("plate", 10);
  plate._scalar_coefficients.Register("electrical_conductivity",
                                      std::make_shared<mfem::ConstantCoefficient>(3.526e7));
  plate._scalar_coefficients.Register("magnetic_permeability",
                                      std::make_shared<mfem::ConstantCoefficient>(M_PI * 4.0e-7));

  hephaestus::Subdomain coil1("coil1", 1);
  coil1._scalar_coefficients.Register("electrical_conductivity",
                                      std::make_shared<mfem::ConstantCoefficient>(1.0));
  coil1._scalar_coefficients.Register("magnetic_permeability",
                                      std::make_shared<mfem::ConstantCoefficient>(M_PI * 4.0e-7));

  hephaestus::Subdomain coil2("coil2", 2);
  coil2._scalar_coefficients.Register("electrical_conductivity",
                                      std::make_shared<mfem::ConstantCoefficient>(1.0));
  coil2._scalar_coefficients.Register("magnetic_permeability",
                                      std::make_shared<mfem::ConstantCoefficient>(M_PI * 4.0e-7));

  hephaestus::Subdomain coil3("coil3", 3);
  coil3._scalar_coefficients.Register("electrical_conductivity",
                                      std::make_shared<mfem::ConstantCoefficient>(1.0));
  coil3._scalar_coefficients.Register("magnetic_permeability",
                                      std::make_shared<mfem::ConstantCoefficient>(M_PI * 4.0e-7));

  hephaestus::Subdomain coil4("coil4", 4);
  coil4._scalar_coefficients.Register("electrical_conductivity",
                                      std::make_shared<mfem::ConstantCoefficient>(1.0));
  coil4._scalar_coefficients.Register("magnetic_permeability",
                                      std::make_shared<mfem::ConstantCoefficient>(M_PI * 4.0e-7));

  hephaestus::Coefficients coefficients(
      std::vector<hephaestus::Subdomain>({air, plate, coil1, coil2, coil3, coil4}));

  coefficients._scalars.Register("I", std::make_shared<mfem::ConstantCoefficient>(1000));

  return coefficients;
}

hephaestus::Sources
defineSources()
{
  hephaestus::InputParameters div_free_source_params;

  // This vector of subdomains will form the coil that we pass to
  // ClosedCoilSolver
  int order = 1;
  int electrode_attr = 11;
  std::string coil_attr = "1 2 3 4";
  mfem::Array<int> coil_domains;
  std::stringstream ss(coil_attr);
  int att;

  while (ss >> att)
    coil_domains.Append(att);

  hephaestus::Sources sources;
  sources.Register("source",
                   std::make_shared<hephaestus::ClosedCoilSolver>("source_grad_phi",
                                                                  "HCurl",
                                                                  "H1",
                                                                  "I",
                                                                  "electrical_conductivity",
                                                                  coil_domains,
                                                                  electrode_attr,
                                                                  true));
  return sources;
}

hephaestus::Outputs
defineOutputs()
{
  hephaestus::Outputs outputs;
  outputs.Register("ParaViewDataCollection",
                   std::make_shared<mfem::ParaViewDataCollection>("Team13ParaView"));
  return outputs;
}

int
main(int argc, char * argv[])
{
  mfem::OptionsParser args(argc, argv);
  args.AddOption(
      &DATA_DIR, "-dataDir", "--data_directory", "Directory storing input data for tests.");
  args.Parse();
  MPI_Init(&argc, &argv);

  hephaestus::logger.set_level(spdlog::level::info);

  // Create Formulation
  auto problem_builder = std::make_unique<hephaestus::NLMagnetostaticFormulation>(
      "magnetic_reluctivity", "magnetic_permeability", "magnetic_vector_potential");

  // Set Mesh
  mfem::Mesh mesh((std::string(DATA_DIR) + std::string("./vac_team13.e")).c_str(), 1, 1);
  auto pmesh = std::make_shared<mfem::ParMesh>(MPI_COMM_WORLD, mesh);
  mfem::H1_FECollection fecm(1, 3);
  mfem::ParFiniteElementSpace pfespace(pmesh.get(), &fecm, 3);
  // Necessary, in case the nodal FE space is not set on the pmesh because it is lowest order.
  pmesh->SetNodalFESpace(&pfespace);

  int par_ref_lvl = 0;
  for (int l = 0; l < par_ref_lvl; ++l)
    pmesh->UniformRefinement();

  problem_builder->SetMesh(pmesh);
  problem_builder->AddFESpace(std::string("H1"), std::string("H1_3D_P1"));
  problem_builder->AddFESpace(std::string("HCurl"), std::string("ND_3D_P1"));
  problem_builder->AddFESpace(std::string("HDiv"), std::string("RT_3D_P0"));
  problem_builder->AddGridFunction(std::string("magnetic_vector_potential"), std::string("HCurl"));
  problem_builder->AddGridFunction(std::string("source_grad_phi"), std::string("HCurl"));
  problem_builder->AddGridFunction(std::string("magnetic_flux_density"), std::string("HDiv"));
  problem_builder->RegisterMagneticFluxDensityAux("magnetic_flux_density");

  hephaestus::Coefficients coefficients = defineCoefficients();
  problem_builder->SetCoefficients(coefficients);

  hephaestus::Sources sources = defineSources();
  problem_builder->SetSources(sources);

  hephaestus::Outputs outputs = defineOutputs();
  problem_builder->SetOutputs(outputs);
  {
    // Call LineSampler to save values
    std::string gridfunction_name("magnetic_flux_density");
    const int num_pts = 100;
    // Mesh bounding box (for the full serial mesh).
    mfem::Vector pos_min(3), pos_max(3);
    pos_min(0) = 1.6 / 2.0 / 1000.0; // 0<x<1.6 (mm) -> halfway point
    pos_max(0) = pos_min(0);
    pos_min(1) = 0.0;
    pos_max(1) = 0.0;
    pos_min(2) = 0.0;
    pos_max(2) = 60.0 / 1000.0;
    std::string csv_name_ab("fluxDensityProfileAB.csv");
    std::shared_ptr<hephaestus::LineSamplerAux> linesamplerwriter_ab =
        std::make_shared<hephaestus::LineSamplerAux>(gridfunction_name,
                                                     pos_min,
                                                     pos_max,
                                                     num_pts,
                                                     csv_name_ab,
                                                     "t_s, x_m, y_m, z_m, B_x_T, B_y_T, B_z_T");
    linesamplerwriter_ab->SetPriority(5);
    problem_builder->AddPostprocessor("LineSamplerWriterAB", linesamplerwriter_ab);

    pos_min(0) = 2.1 / 1000.0; // 0<x<1.6 (mm) -> halfway point
    pos_max(0) = 122.1 / 1000.0;
    pos_min(1) = (65.0 + 15.0) / 2.0 / 1000.0;
    pos_max(1) = (65.0 + 15.0) / 2.0 / 1000.0;
    pos_min(2) = (63.2 + 60.0) / 2.0 / 1000.0;
    pos_max(2) = (63.2 + 60.0) / 2.0 / 1000.0;
    std::string csv_name_cd("fluxDensityProfileCD.csv");
    std::shared_ptr<hephaestus::LineSamplerAux> linesamplerwriter_cd =
        std::make_shared<hephaestus::LineSamplerAux>(gridfunction_name,
                                                     pos_min,
                                                     pos_max,
                                                     num_pts,
                                                     csv_name_cd,
                                                     "t_s, x_m, y_m, z_m, B_x_T, B_y_T, B_z_T");
    linesamplerwriter_cd->SetPriority(5);
    problem_builder->AddPostprocessor("LineSamplerWriterCD", linesamplerwriter_cd);

    pos_min(0) = (125.3 + 122.1) / 2.0 / 1000.0;
    pos_max(0) = (125.3 + 122.1) / 2.0 / 1000.0;
    pos_min(1) = (65.0 + 15.0) / 2.0 / 1000.0;
    pos_max(1) = (65.0 + 15.0) / 2.0 / 1000.0;
    pos_min(2) = 0.0;
    pos_max(2) = (60.0) / 1000.0;
    std::string csv_name_ef("fluxDensityProfileEF.csv");
    std::shared_ptr<hephaestus::LineSamplerAux> linesamplerwriter_ef =
        std::make_shared<hephaestus::LineSamplerAux>(gridfunction_name,
                                                     pos_min,
                                                     pos_max,
                                                     num_pts,
                                                     csv_name_ef,
                                                     "t_s, x_m, y_m, z_m, B_x_T, B_y_T, B_z_T");
    linesamplerwriter_ef->SetPriority(5);
    problem_builder->AddPostprocessor("LineSamplerWriterEF", linesamplerwriter_ef);

    hephaestus::InputParameters solver_options;
    solver_options.SetParam("Tolerance", float(1.0e-13));
    solver_options.SetParam("AbsTolerance", float(1.0e-16));
    solver_options.SetParam("MaxIter", (unsigned int)500);
    problem_builder->SetSolverOptions(solver_options);

    problem_builder->FinalizeProblem();

    auto problem = problem_builder->ReturnProblem();
    hephaestus::InputParameters exec_params;
    exec_params.SetParam("VisualisationSteps", int(1));
    exec_params.SetParam("UseGLVis", true);
    exec_params.SetParam("Problem", static_cast<hephaestus::SteadyStateProblem *>(problem.get()));

    auto executioner = std::make_unique<hephaestus::SteadyExecutioner>(exec_params);

    hephaestus::logger.info("Created executioner");
    executioner->Execute();
  }

  MPI_Finalize();
}