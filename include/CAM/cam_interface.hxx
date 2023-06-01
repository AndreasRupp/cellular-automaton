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
  unsigned int rand_seed;

 public:
  CAMInterface(double _jump_parameter_composites = 5)
  {
    CAM::jump_parameter_composite = _jump_parameter_composites;
  }
  /*!*********************************************************************************************
   * \brief Creates domain with building units containing only one cell
   *
   * \param _porosity The percentage of void space, not occupied by solid/bu.
   * \param _jump_parameter How far individual particles are allowed to jump.
   * \param random_seed If given, sets random seed to given seed.
   ************************************************************************************************/
  void place_single_cell_bu_randomly(double _porosity = 0.5,
                                     double _jump_parameter = 1,
                                     unsigned int _random_seed = 0)
  {
    if (_random_seed == 0)
    {
      rand_seed = std::chrono::system_clock::now().time_since_epoch().count();
      std::srand(rand_seed);
    }

    else
      std::srand(_random_seed);

    unsigned int n_particles = (1. - _porosity) * domain.n_fields_;
    unsigned int position = std::rand() % (domain.n_fields_);
    for (unsigned int i = 0; i < n_particles; ++i)
    {
      while (domain.domain_fields[position] != 0)
        position = std::rand() % (domain.n_fields_);
      domain.place_bu(CAM::BuildingUnit<nx>::create_single_cell_bu(_jump_parameter), position);
    }
  }
  /*!*********************************************************************************************
   * \brief Places HyperSphere (2D: Circle, 3D: Sphere) into domain
   *
   * \param _position field index of center point; no position is given (-1) -> random position
   * \param _radius all cells are completely within the radius
   * \param _jump_parameter How far sphere is allowed to jump.
   * \return true: sphere is placed
   * \return false: sphere could not be placed. all its cells are removed
   ************************************************************************************************/
  bool place_sphere(int _position = -1, double _jump_parameter = 1, double _radius = 1)
  {
    const CAM::BuildingUnit<nx> bu =
      CAM::BuildingUnit<nx>::create_hyper_sphere(_jump_parameter, _radius);
    return domain.place_bu(bu, _position);
  }
  /*!*********************************************************************************************
   * \brief Place a limited hyper plane (2D: Rectangle, 3D: cuboid) into domain
   *
   * \param _position field index of center point; no position is given (-1) -> random position
   * \param _extent size of plane in each dimension
   * \param _jump_parameter How far hyper plane is allowed to jump.
   * \return true: plane is placed
   * \return false: plane could not be placed. all its cells are removed
   ************************************************************************************************/
  bool place_plane(int _position = -1,
                   double _jump_parameter = 1,
                   std::vector<unsigned int> _extent_v = std::vector<unsigned int>(nx.size(), 0))
  {
    std::array<unsigned int, nx.size()> extent_a;
    for (unsigned int i = 0; i < _extent_v.size(); i++)
      extent_a[i] = _extent_v[i];

    const BuildingUnit<nx> bu = BuildingUnit<nx>::create_hyper_plane(_jump_parameter, extent_a);
    return domain.place_bu(bu, _position);
  }
  // TODO parameterize this function
  void place_particles()
  {
    rand_seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::srand(rand_seed);
    unsigned int random_point;

    // std::vector<unsigned int> stencil1 = CAM::ParticleBU<nx>::get_custom_particle_stencil(0, 30,
    // 20, 100); std::vector<unsigned int> stencil2 =
    // CAM::ParticleBU<nx>::get_custom_particle_stencil(0, 20, 10, 200);

    std::vector<unsigned int> stencil1 =
      CAM::BuildingUnit<nx>::get_custom_particle_stencil(0, 10, 5, 10);
    std::vector<unsigned int> stencil2 =
      CAM::BuildingUnit<nx>::get_custom_particle_stencil(0, 5, 1, 1);
    unsigned int index = 0;
    while (index < 10)
    {
      random_point = 0;
      for (unsigned int i = 0; i < nx.size(); ++i)
      {
        random_point += std::rand() % nx[i] * direct_neigh<nx>(2 * i + 1);
      }
      index = index + (unsigned int)domain.place_bu(
                        CAM::BuildingUnit<nx>::create_custom_bu(5, stencil1), random_point);
    }
    while (index < 20)
    {
      random_point = 0;
      for (unsigned int i = 0; i < nx.size(); ++i)
      {
        random_point += std::rand() % nx[i] * direct_neigh<nx>(2 * i + 1);
      }
      index = index + (unsigned int)domain.place_bu(
                        CAM::BuildingUnit<nx>::create_custom_bu(5, stencil2), random_point);
    }
  }

  void print_array() { domain.print_array(); }

  void do_cam() { CAM::CellularAutomaton<nx, fields_array_t>::apply(domain); }

  const fields_array_t& fields() const { return domain.domain_fields; }

  std::array<double, 12> eval_measures()
  {
    return CAM::Evaluation<nx, fields_array_t>::eval_measures(&domain);
  }

  double average_particle_size()
  {
    return CAM::Evaluation<nx, fields_array_t>::average_particle_size(&domain);
  }
  std::vector<unsigned int> particle_size_distribution()
  {
    return CAM::Evaluation<nx, fields_array_t>::particle_size_distribution(&domain);
  }
  unsigned int bulk_distance(const fields_array_t& domain_a, const fields_array_t& domain_b)
  {
    return CAM::Evaluation<nx, fields_array_t>::bulk_distance(domain_a, domain_b);
  }
};
}  // namespace CAM
