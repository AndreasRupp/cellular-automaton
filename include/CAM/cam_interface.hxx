/*!*********************************************************************************************
 * \file cam_interface.hxx
 * \brief Interface for Python
 **********************************************************************************************/
#pragma once
#include <CAM/building_units.hxx>
#include <CAM/cellular_automaton.hxx>
#include <CAM/composite.hxx>
#include <CAM/domain.hxx>
#include <CAM/evaluation.hxx>
#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
namespace CAM
{
template <auto nx, typename fields_array_t = std::array<unsigned int, CAM::n_fields<nx>()>>
class CAMInterface
{
  CAM::Domain<nx, fields_array_t> domain;
  CAM::Evaluation<nx, fields_array_t> evaluation;

 public:
  CAMInterface(double _jump_parameter_composites = 5)
  {
    CAM::jump_parameter_composite = _jump_parameter_composites;
    evaluation = CAM::Evaluation<nx, fields_array_t>(&domain);
  }
  ~CAMInterface() {}

  void place_single_cell_bu_randomly(double _porosity = 0.5,
                                     double _jump_parameter = 1,
                                     unsigned int _random_seed = 0)
  {
    domain.place_single_cell_bu_randomly(_porosity, _jump_parameter, _random_seed);
  }
  bool place_sphere(int _position = -1, double _radius = 1, double _jump_parameter = 1)
  {
    return domain.place_sphere(_position, _radius, _jump_parameter);
  }
  bool place_plane(int _position = -1,
                   double _jump_parameter = 1,
                   std::vector<unsigned int> _extent_v = std::vector<unsigned int>(nx.size(), 0))
  {
    std::array<unsigned int, nx.size()> extent_a;
    for (unsigned int i = 0; i < _extent_v.size(); i++)
      extent_a[i] = _extent_v[i];
    return domain.place_plane(_position, _jump_parameter, extent_a);
  }
  void place_particles()
  {
    // TODO parameterize this function
    unsigned int rand_seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::srand(rand_seed);
    unsigned int randomPoint, centerpoint = 0;

    // std::vector<unsigned int> stencil1 = CAM::ParticleBU<nx>::get_stencil(0, 30, 20, 100);
    // std::vector<unsigned int> stencil2 = CAM::ParticleBU<nx>::get_stencil(0, 20, 10, 200);

    std::vector<unsigned int> stencil1 = CAM::BuildingUnit<nx>::get_stencil(0, 10, 5, 10);
    std::vector<unsigned int> stencil2 = CAM::BuildingUnit<nx>::get_stencil(0, 5, 1, 1);
    unsigned int index = 0;
    while (index < 10)
    {
      randomPoint = 0;
      for (unsigned int i = 0; i < nx.size(); ++i)
      {
        centerpoint += nx[i] / 2 * direct_neigh<nx>(2 * i + 1);
        randomPoint += std::rand() % nx[i] * direct_neigh<nx>(2 * i + 1);
      }
      index = index + (unsigned int)domain.place_bu(
                        CAM::BuildingUnit<nx>(index + 1, 5, randomPoint, stencil1));
    }
    while (index < 20)
    {
      randomPoint = 0;
      for (unsigned int i = 0; i < nx.size(); ++i)
      {
        centerpoint += nx[i] / 2 * direct_neigh<nx>(2 * i + 1);
        randomPoint += std::rand() % nx[i] * direct_neigh<nx>(2 * i + 1);
      }
      index = index + (unsigned int)domain.place_bu(
                        CAM::BuildingUnit<nx>(index + 1, 5, randomPoint, stencil2));
    }
  }

  void print_array() { domain.print_array(); }

  void do_cam() { CAM::CellularAutomaton<nx, fields_array_t>::apply(domain); }

  const fields_array_t& fields() const { return domain.domain_fields; }
  std::array<double, 12> eval_measures() { return evaluation.eval_measures(); }

  double average_particle_size() { return evaluation.average_particle_size(); }
  std::vector<unsigned int> particle_size_distribution()
  {
    return evaluation.particle_size_distribution();
  }
  unsigned int bulk_distance(const fields_array_t& domain_a, const fields_array_t& domain_b)
  {
    return CAM::Evaluation<nx, fields_array_t>::bulk_distance(domain_a, domain_b);
  }
};
}  // namespace CAM
