#pragma once
#include <CAM/building_units.hxx>
#include <CAM/cellular_automaton.hxx>
#include <CAM/domain.hxx>
#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
namespace CAM
{
template <auto nx, typename fields_array_t = std::array<CAM::fieldNumbers_t, CAM::n_fields<nx>()>>
class CAMInterface
{
  CAM::Domain<nx, fields_array_t> domain;

 public:
  CAMInterface() {}
  ~CAMInterface() {}

  void placeBURandomly(double _porosity = 0.5,
                       double _jump_parameter = 1,
                       unsigned int _random_seed = 0)
  {
    domain.placeBURandomly(_porosity, _jump_parameter, _random_seed);
  }
  void placeSphere(double _jump_parameter = 1, unsigned int _random_seed = 0)
  {
    domain.placeSphere(_jump_parameter, _random_seed);
  }

  void print_array() { domain.print_array(); }

  void doCAM() { CAM::CellularAutomaton<nx, fields_array_t>::apply(domain); }

  const fields_array_t& fields() const { return domain.domainFields; }
  std::array<double, 12> eval_measures() { return domain.eval_measures(); }
};
}  // namespace CAM