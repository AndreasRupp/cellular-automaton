/**
 * @file cam_interface.hxx
 * @brief Interface for Python
 *
 *
 */
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
  CAMInterface(double _jump_parameter_composites = 5)
  {
    domain.jump_parameter_composites = _jump_parameter_composites;
  }
  ~CAMInterface() {}

  void placeSingleCellBURandomly(double _porosity = 0.5,
                                 double _jump_parameter = 1,
                                 unsigned int _random_seed = 0)
  {
    domain.placeSingleCellBURandomly(_porosity, _jump_parameter, _random_seed);
  }
  bool placeSphere(int _position = -1, double _radius = 1, double _jump_parameter = 1)
  {
    return domain.placeSphere(_position, _radius, _jump_parameter);
  }
  bool placePlane(int _position = -1,
                  std::vector<unsigned int> _extent = std::vector<unsigned int>(nx.size(), 0),
                  double _jump_parameter = 1)
  {
    return domain.placePlane(_position, _extent, _jump_parameter);
  }
  void placeParticles()
  {
    // TODO parameterize this function
    unsigned int rand_seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::srand(rand_seed);
    unsigned int centerpoint = 0;

    std::vector<unsigned int> stencil1 = CAM::ParticleBU<nx>::getStencil(0, 30, 20, 100);
    std::vector<unsigned int> stencil2 = CAM::ParticleBU<nx>::getStencil(0, 20, 10, 200);
    for (int a = 0; a < 10; a++)
    {
      unsigned int randomPoint = 0;
      for (unsigned int i = 0; i < nx.size(); ++i)
      {
        centerpoint += nx[i] / 2 * direct_neigh<nx>(2 * i + 1);
        randomPoint += std::rand() % nx[i] * direct_neigh<nx>(2 * i + 1);
      }
      CAM::ParticleBU<nx>* particle = new CAM::ParticleBU<nx>(a + 1, 5, randomPoint, stencil1);
      domain.placeBU(particle);
    }

    for (int a = 10; a < 20; a++)
    {
      unsigned int randomPoint = 0;
      for (unsigned int i = 0; i < nx.size(); ++i)
      {
        centerpoint += nx[i] / 2 * direct_neigh<nx>(2 * i + 1);
        randomPoint += std::rand() % nx[i] * direct_neigh<nx>(2 * i + 1);
      }
      CAM::ParticleBU<nx>* particle = new CAM::ParticleBU<nx>(a + 1, 5, randomPoint, stencil2);
      domain.placeBU(particle);
      ;
    }
  }

  void print_array() { domain.print_array(); }

  void doCAM() { CAM::CellularAutomaton<nx, fields_array_t>::apply(domain); }

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