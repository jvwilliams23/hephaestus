#pragma once
#include <fstream>
#include <iostream>
#include <memory>

#include "../common/pfem_extras.hpp"
#include "gridfunctions.hpp"
#include "mesh_extras.hpp"

namespace hephaestus {

class Outputs : public mfem::NamedFieldsMap<mfem::DataCollection> {
  friend class ProblemBuilder;

public:
  Outputs();
  Outputs(hephaestus::GridFunctions &gridfunctions);

  void SetGridFunctions(hephaestus::GridFunctions &gridfunctions) {
    _gridfunctions = &gridfunctions;
  }

  // Enable GLVis streams for visualisation
  void EnableGLVis(bool &use_glvis) { _use_glvis = use_glvis; }

  // Reset Outputs and re-register output fields from GridFunctions
  void Reset() {
    // Reset cycle counter
    _cycle = 0;
    // Set up DataCollections to track fields of interest.
    RegisterOutputFields();
    // Write initial fields to disk
    WriteOutputFields(0.0);

    // Initialize GLVis _use_glvis and send the initial condition
    // by socket to a GLVis server.
    if (_use_glvis) {
      InitializeGLVis(_my_rank);
      DisplayToGLVis();
    }
  }

  // Write outputs out to requested streams
  void Write(double t = 1.0) {
    // Wait for all ranks to finish updating their solution before output
    MPI_Barrier(_my_comm);
    // Update cycle counter
    _cycle++;
    // Output timestep summary to console
    WriteConsoleSummary(_my_rank, t);
    if (_use_glvis) {
      DisplayToGLVis();
    }
    // Save output fields at timestep to DataCollections
    WriteOutputFields(t);
  }

private:
  std::map<std::string, mfem::socketstream *> socks_;
  hephaestus::GridFunctions *_gridfunctions;
  int _cycle{0};
  bool _use_glvis{false};
  MPI_Comm _my_comm{MPI_COMM_WORLD};
  int _n_ranks, _my_rank;

  // Initialize Outputs with Gridfunctions; used in ProblemBuilder
  void Init(hephaestus::GridFunctions &gridfunctions) {
    SetGridFunctions(gridfunctions);
    Reset();
  }

  // Register fields (gridfunctions) to write to DataCollections
  void RegisterOutputFields() {
    for (auto output = begin(); output != end(); ++output) {
      auto const &dc_(output->second);
      mfem::ParMesh *pmesh_(
          _gridfunctions->begin()->second->ParFESpace()->GetParMesh());
      dc_->SetMesh(pmesh_);
      for (auto var = _gridfunctions->begin(); var != _gridfunctions->end();
           ++var) {
        dc_->RegisterField(var->first, var->second);
      }
    }
  }

  // Write out fields (gridfunctions) to DataCollections
  void WriteOutputFields(double t) {
    // Write fields to disk
    for (auto output = begin(); output != end(); ++output) {
      auto const &dc_(output->second);
      dc_->SetCycle(_cycle);
      dc_->SetTime(t);
      dc_->Save();
    }
  }

  // Write out summary of last timestep to console
  void WriteConsoleSummary(int _my_rank, double t) {
    if (_my_rank == 0) {
      std::cout << std::fixed;
      std::cout << "step " << std::setw(6) << _cycle
                << ",\tt = " << std::setw(6) << std::setprecision(3) << t
                << std::endl;
    }
  }

  // Initialize GLVis sockets and fields
  void InitializeGLVis(int _my_rank) {
    if (_my_rank == 0) {
      std::cout << "Opening GLVis sockets." << std::endl;
    }

    for (auto var = _gridfunctions->begin(); var != _gridfunctions->end();
         ++var) {
      socks_[var->first] = new mfem::socketstream;
      socks_[var->first]->precision(8);
    }

    if (_my_rank == 0) {
      std::cout << "GLVis sockets open." << std::endl;
    }
  }

  // Update GLVis display of output fields (gridfunctions)
  void DisplayToGLVis() {
    char vishost[] = "localhost";
    int visport = 19916;

    int Wx = 0, Wy = 0;                 // window position
    int Ww = 350, Wh = 350;             // window size
    int offx = Ww + 10, offy = Wh + 45; // window offsets

    for (auto var = _gridfunctions->begin(); var != _gridfunctions->end();
         ++var) {
      mfem::common::VisualizeField(*socks_[var->first], vishost, visport,
                                   *(var->second), (var->first).c_str(), Wx, Wy,
                                   Ww, Wh);
      Wx += offx;
    }
  }
};

} // namespace hephaestus
