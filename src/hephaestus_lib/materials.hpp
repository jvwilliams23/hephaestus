#pragma once
#include "mesh_extras.hpp"
#include <fstream>
#include <iostream>
#include <memory>
#include <unordered_set>

namespace hephaestus {

double prodFunc(double a, double b);
double fracFunc(double a, double b);

class Subdomain {
public:
  Subdomain(const std::string &name_, int id_);

  std::string name;
  int id;
  mfem::NamedFieldsMap<mfem::Coefficient> property_map;
};

// Coefficients - stores all scalar and vector coefficients
//--SetTime
//--scalar_property_map
//--vector_property_map

// Stores all coefficients defined over
class Coefficients {
  double t; // Time at which time-dependent coefficients are evaluated
public:
  Coefficients();
  Coefficients(std::vector<Subdomain> subdomains_);
  void SetTime(double t);
  void AddGlobalCoefficientsFromSubdomains();

  mfem::NamedFieldsMap<mfem::Coefficient> scalar_property_map;
  mfem::NamedFieldsMap<mfem::VectorCoefficient> vector_property_map;
  std::vector<Subdomain> subdomains;
};

} // namespace hephaestus
