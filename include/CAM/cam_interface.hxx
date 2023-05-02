/*!*********************************************************************************************
 * \file cam_interface.hxx
 * \brief Interface for Python
 **********************************************************************************************/
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
template <auto nx, typename fields_array_t = std::array<unsigned int, CAM::n_fields<nx>()>>
class CAMInterface
{
  CAM::Domain<nx, fields_array_t> domain;

 public:
  CAMInterface(double _jump_parameter_composites = 5)
  {
    domain.jump_parameter_composites = _jump_parameter_composites;
  }
  ~CAMInterface() {}

  void place_singleCellBU_randomly(double _porosity = 0.5,
                                   double _jump_parameter = 1,
                                   unsigned int _random_seed = 0)
  {
    domain.place_singleCellBU_randomly(_porosity, _jump_parameter, _random_seed);
  }
  bool place_sphere(int _position = -1, double _radius = 1, double _jump_parameter = 1)
  {
    return domain.place_sphere(_position, _radius, _jump_parameter);
  }
  bool place_plane(int _position = -1,
                   std::vector<unsigned int> _extent = std::vector<unsigned int>(nx.size(), 0),
                   double _jump_parameter = 1)
  {
    return domain.place_plane(_position, _extent, _jump_parameter);
  }
  void place_particles()
  {
    // TODO parameterize this function
    unsigned int rand_seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::srand(rand_seed);
    unsigned int randomPoint, centerpoint = 0;

    // std::vector<unsigned int> stencil1 = CAM::ParticleBU<nx>::get_stencil(0, 30, 20, 100);
    // std::vector<unsigned int> stencil2 = CAM::ParticleBU<nx>::get_stencil(0, 20, 10, 200);

    std::vector<unsigned int> stencil1 = CAM::ParticleBU<nx>::get_stencil(0, 10, 5, 10);
    std::vector<unsigned int> stencil2 = CAM::ParticleBU<nx>::get_stencil(0, 5, 1, 1);
    for (int a = 0; a < 10; a++)
    {
      randomPoint = 0;
      for (unsigned int i = 0; i < nx.size(); ++i)
      {
        centerpoint += nx[i] / 2 * direct_neigh<nx>(2 * i + 1);
        randomPoint += std::rand() % nx[i] * direct_neigh<nx>(2 * i + 1);
      }
      CAM::ParticleBU<nx>* particle = new CAM::ParticleBU<nx>(a + 1, 5, randomPoint, stencil1);
      domain.place_BU(particle);
    }

    for (int a = 10; a < 20; a++)
    {
      randomPoint = 0;
      for (unsigned int i = 0; i < nx.size(); ++i)
      {
        centerpoint += nx[i] / 2 * direct_neigh<nx>(2 * i + 1);
        randomPoint += std::rand() % nx[i] * direct_neigh<nx>(2 * i + 1);
      }
      CAM::ParticleBU<nx>* particle = new CAM::ParticleBU<nx>(a + 1, 5, randomPoint, stencil2);
      domain.place_BU(particle);
      ;
    }
  }

  void print_array() { domain.print_array(); }

  void do_CAM() { CAM::CellularAutomaton<nx, fields_array_t>::apply(domain); }

  const fields_array_t& fields() const { return domain.domainFields; }
  std::array<double, 12> eval_measures() { return domain.eval_measures(); }

  double average_particle_size() { return domain.average_particle_size(); }
  std::vector<unsigned int> particle_size_distribution()
  {
    return domain.particle_size_distribution();
  }
  unsigned int bulk_distance(const fields_array_t& domain_a, const fields_array_t& domain_b)
  {
    return CAM::Domain<nx, fields_array_t>::bulk_distance(domain_a, domain_b);
  }
};
}  // namespace CAM
