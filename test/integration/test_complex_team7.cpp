#include "auxsolvers.hpp"
#include "steady_executioner.hpp"

#include "hephaestus.hpp"

#include <catch2/catch_test_macros.hpp>

extern const char * DATA_DIR;

class TestComplexTeam7
{
protected:
  static void SourceCurrent(const mfem::Vector & xv, mfem::Vector & J)
  {
    double x0(194e-3);  // Coil centre x coordinate
    double y0(100e-3);  // Coil centre y coordinate
    double a(50e-3);    // Coil thickness
    double i0(2742);    // Coil current in Ampere-turns
    double s(2.5e-3);   // Coil cross sectional area
    double freq(200.0); // Frequency in Hz

    double x = xv(0);
    double y = xv(1);

    // Current density magnitude
    double jmag = (i0 / s);

    // Calculate x component of current density unit vector
    if (abs(x - x0) < a)
    {
      J(0) = -(y - y0) / abs(y - y0);
    }
    else if (abs(y - y0) < a)
    {
      J(0) = 0.0;
    }
    else
    {
      J(0) =
          -(y - (y0 + a * ((y - y0) / abs(y - y0)))) /
          hypot(x - (x0 + a * ((x - x0) / abs(x - x0))), y - (y0 + a * ((y - y0) / abs(y - y0))));
    }

    // Calculate y component of current density unit vector
    if (abs(y - y0) < a)
    {
      J(1) = (x - x0) / abs(x - x0);
    }
    else if (abs(x - x0) < a)
    {
      J(1) = 0.0;
    }
    else
    {
      J(1) =
          (x - (x0 + a * ((x - x0) / abs(x - x0)))) /
          hypot(x - (x0 + a * ((x - x0) / abs(x - x0))), y - (y0 + a * ((y - y0) / abs(y - y0))));
    }

    // Calculate z component of current density unit vector
    J(2) = 0.0;

    // Scale by current density magnitude
    J *= jmag;
  }

  hephaestus::InputParameters TestParams()
  {
    hephaestus::Subdomain air("air", 1);
    air._scalar_coefficients.Register(
        "electrical_conductivity", new mfem::ConstantCoefficient(1.0), true);
    hephaestus::Subdomain plate("plate", 2);
    plate._scalar_coefficients.Register(
        "electrical_conductivity", new mfem::ConstantCoefficient(3.526e7), true);
    hephaestus::Subdomain coil1("coil1", 3);
    coil1._scalar_coefficients.Register(
        "electrical_conductivity", new mfem::ConstantCoefficient(1.0), true);
    hephaestus::Subdomain coil2("coil2", 4);
    coil2._scalar_coefficients.Register(
        "electrical_conductivity", new mfem::ConstantCoefficient(1.0), true);
    hephaestus::Subdomain coil3("coil3", 5);
    coil3._scalar_coefficients.Register(
        "electrical_conductivity", new mfem::ConstantCoefficient(1.0), true);
    hephaestus::Subdomain coil4("coil4", 6);
    coil4._scalar_coefficients.Register(
        "electrical_conductivity", new mfem::ConstantCoefficient(1.0), true);

    hephaestus::Coefficients coefficients(
        std::vector<hephaestus::Subdomain>({air, plate, coil1, coil2, coil3, coil4}));

    coefficients._scalars.Register("frequency", new mfem::ConstantCoefficient(200.0), true);
    coefficients._scalars.Register(
        "magnetic_permeability", new mfem::ConstantCoefficient(M_PI * 4.0e-7), true);
    coefficients._scalars.Register(
        "dielectric_permittivity", new mfem::ConstantCoefficient(0.0), true);

    hephaestus::BCMap bc_map;

    mfem::Mesh mesh((std::string(DATA_DIR) + std::string("./team7.g")).c_str(), 1, 1);

    hephaestus::Outputs outputs;
    outputs.Register(
        "VisItDataCollection", new mfem::VisItDataCollection("ComplexMaxwellTeam7VisIt"), true);
    outputs.Register("ParaViewDataCollection",
                     new mfem::ParaViewDataCollection("ComplexMaxwellTeam7ParaView"),
                     true);

    hephaestus::AuxSolvers postprocessors;
    hephaestus::AuxSolvers preprocessors;

    hephaestus::Sources sources;
    auto * j_src_coef = new mfem::VectorFunctionCoefficient(3, SourceCurrent);
    mfem::Array<mfem::VectorCoefficient *> sourcecoefs(4);
    sourcecoefs[0] = j_src_coef;
    sourcecoefs[1] = j_src_coef;
    sourcecoefs[2] = j_src_coef;
    sourcecoefs[3] = j_src_coef;
    mfem::Array<int> coilsegments(4);
    coilsegments[0] = 3;
    coilsegments[1] = 4;
    coilsegments[2] = 5;
    coilsegments[3] = 6;
    auto * j_src_restricted = new mfem::PWVectorCoefficient(3, coilsegments, sourcecoefs);
    coefficients._vectors.Register("source", j_src_restricted, true);

    hephaestus::InputParameters div_free_source_params;
    div_free_source_params.SetParam("SourceName", std::string("source"));
    div_free_source_params.SetParam("PotentialName", std::string("electric_potential"));
    div_free_source_params.SetParam("ConductivityCoefName", std::string("electrical_conductivity"));
    div_free_source_params.SetParam("HCurlFESpaceName", std::string("HCurl"));
    div_free_source_params.SetParam("H1FESpaceName", std::string("H1"));
    hephaestus::InputParameters current_solver_options;
    current_solver_options.SetParam("Tolerance", float(1.0e-12));
    current_solver_options.SetParam("MaxIter", (unsigned int)200);
    current_solver_options.SetParam("PrintLevel", 1);
    div_free_source_params.SetParam("SolverOptions", current_solver_options);
    sources.Register("source", new hephaestus::DivFreeSource(div_free_source_params), true);

    hephaestus::InputParameters solver_options;
    solver_options.SetParam("Tolerance", float(1.0e-16));
    solver_options.SetParam("MaxIter", (unsigned int)1000);
    solver_options.SetParam("PrintLevel", 0);

    hephaestus::InputParameters params;
    params.SetParam("UseGLVis", true);

    params.SetParam("Mesh", mfem::ParMesh(MPI_COMM_WORLD, mesh));
    params.SetParam("BoundaryConditions", bc_map);
    params.SetParam("Coefficients", coefficients);
    params.SetParam("PreProcessors", preprocessors);
    params.SetParam("PostProcessors", postprocessors);
    params.SetParam("Outputs", outputs);
    params.SetParam("Sources", sources);
    params.SetParam("SolverOptions", solver_options);
    std::cout << "Created params ";
    return params;
  }
};

TEST_CASE_METHOD(TestComplexTeam7, "TestComplexTeam7", "[CheckRun]")
{
  hephaestus::InputParameters params(TestParams());
  auto pmesh = std::make_shared<mfem::ParMesh>(params.GetParam<mfem::ParMesh>("Mesh"));

  auto bc_map(params.GetParam<hephaestus::BCMap>("BoundaryConditions"));
  auto coefficients(params.GetParam<hephaestus::Coefficients>("Coefficients"));
  auto preprocessors(params.GetParam<hephaestus::AuxSolvers>("PreProcessors"));
  auto postprocessors(params.GetParam<hephaestus::AuxSolvers>("PostProcessors"));
  auto sources(params.GetParam<hephaestus::Sources>("Sources"));
  auto outputs(params.GetParam<hephaestus::Outputs>("Outputs"));
  auto solver_options(params.GetOptionalParam<hephaestus::InputParameters>(
      "SolverOptions", hephaestus::InputParameters()));

  auto problem_builder =
      std::make_unique<hephaestus::ComplexAFormulation>("magnetic_reluctivity",
                                                        "electrical_conductivity",
                                                        "dielectric_permittivity",
                                                        "frequency",
                                                        "magnetic_vector_potential",
                                                        "magnetic_vector_potential_real",
                                                        "magnetic_vector_potential_imag");
  problem_builder->SetMesh(pmesh);
  problem_builder->AddFESpace(std::string("HCurl"), std::string("ND_3D_P1"));
  problem_builder->AddFESpace(std::string("HDiv"), std::string("RT_3D_P0"));
  problem_builder->AddFESpace(std::string("H1"), std::string("H1_3D_P1"));
  problem_builder->AddGridFunction(std::string("magnetic_vector_potential_real"),
                                   std::string("HCurl"));
  problem_builder->AddGridFunction(std::string("magnetic_vector_potential_imag"),
                                   std::string("HCurl"));
  problem_builder->AddGridFunction(std::string("magnetic_flux_density_real"), std::string("HDiv"));
  problem_builder->AddGridFunction(std::string("magnetic_flux_density_imag"), std::string("HDiv"));
  problem_builder->SetBoundaryConditions(bc_map);
  problem_builder->SetAuxSolvers(preprocessors);
  problem_builder->SetCoefficients(coefficients);
  problem_builder->SetPostprocessors(postprocessors);
  problem_builder->SetSources(sources);
  problem_builder->SetOutputs(outputs);
  problem_builder->SetSolverOptions(solver_options);

  hephaestus::ProblemBuildSequencer sequencer(problem_builder.get());
  sequencer.ConstructOperatorProblem();
  auto problem = problem_builder->ReturnProblem();

  hephaestus::InputParameters exec_params;
  exec_params.SetParam("Problem", problem.get());

  auto executioner = std::make_unique<hephaestus::SteadyExecutioner>(exec_params);

  std::cout << "Created exec ";

  executioner->Execute();
}
