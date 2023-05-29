/*!*********************************************************************************************
 * \file domain.hxx
 *
 * \brief This class implements a domain and holds all its properties
 * TODO Mayber outsource all evaluation functions to own class
 *
 * \tparam  nx              The size of a row for each dimension of the matrix
 * \tparam  n_fields        Size of the domain
 *
 ************************************************************************************************/
#pragma once

#include <CAM/building_units.hxx>
#include <CAM/composite.hxx>
#include <CAM/utils.hxx>
#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
namespace CAM
{
/*!***********************************************************************************************
 * \brief   Describes one connected set of bulk cells/fields.
 * Information gathered only by analysing the domain independ of defintion of bu and composites
 ************************************************************************************************/
struct Particle
{
  Particle(const std::vector<unsigned int>& _field_indices,
           const std::vector<unsigned int>& _numbers)
  {
    field_indices = _field_indices;
    numbers = _numbers;
  }
  /*!*********************************************************************************************
   * \brief   Indices of building units contained in a particle
   * min_size = 1, max_size = number of building units in domain
   **********************************************************************************************/
  std::vector<unsigned int> numbers;
  /*!*********************************************************************************************
   * \brief   Location of the particle.
   **********************************************************************************************/
  std::vector<unsigned int> field_indices;
};
template <auto nx, typename fields_array_t>
class Domain
{
  static const unsigned int dim = nx.size();

 public:
  static constexpr unsigned int n_fields_ = n_fields<nx>();
  // storing index for new particle
  unsigned int index_bu_max;

  ~Domain() {}
  /*!*********************************************************************************************
   * \brief Construct a new Domain object
   *
   * \param _jump_parameter_composites How far composites particles are allowed to jump.
   ************************************************************************************************/
  Domain(const double _jump_parameter_composites = -1)
  {
    if (_jump_parameter_composites != -1.)
      CAM::jump_parameter_composite = _jump_parameter_composites;
    else
      CAM::jump_parameter_composite = 5;  // double _jump_parameter_composites = 1.0
    if constexpr (std::is_same<fields_array_t,
                               std::vector<typename fields_array_t::value_type>>::value)
      domain_fields.resize(n_fields_, 0);
    else
    {
      static_assert(
        std::is_same<fields_array_t,
                     std::array<typename fields_array_t::value_type, n_fields<nx>()>>::value,
        "The fields array has incorrect size");
      domain_fields.fill(0);
    }
    index_bu_max = 0;
  }
  /*!*********************************************************************************************
   * \brief Place any building unit (bu) into the domain
   * Basis function for placement for all special BUs
   * \param _unit building unit
   * \return true: cells pf bu are marked with individual index, bu is stored in bu vector
   * \return false: bu could not be placed. all its cells are removed
   ************************************************************************************************/
  bool constexpr place_bu(CAM::BuildingUnit<nx> _unit)
  {
    bool success = true;

    std::vector<unsigned int> fields = _unit.get_field_indices();
    std::for_each(fields.begin(), fields.end(),
                  [&](unsigned int field)
                  {
                    if (domain_fields[field] == 0)
                    {
                      success = success && true;
                      domain_fields[field] = _unit.number;
                    }
                    else
                    {
                      success = 0;
                    }
                  });

    if (success == 0)
    {
      std::for_each(fields.begin(), fields.end(),
                    [&](unsigned int field)
                    {
                      if (domain_fields[field] == _unit.number)
                        domain_fields[field] = 0;
                    });
    }
    else
    {
      building_units.push_back(_unit);
      index_bu_max++;
    }
    return success;
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
  bool constexpr place_sphere(const int _position = -1,
                              const double _radius = 1,
                              const double _jump_parameter = 1)
  {
    int position = _position;
    if (_position == -1)
    {
      rand_seed = std::chrono::system_clock::now().time_since_epoch().count();
      std::srand(rand_seed);
      position = std::rand() % (n_fields_);
    }
    bool success = place_bu(BuildingUnit<nx>(index_bu_max + 1, _jump_parameter, position, _radius));
    return success;
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
  bool constexpr place_plane(
    const int _position = -1,
    const double _jump_parameter = 1,
    const std::array<unsigned int, nx.size()>& _extent = std::array<unsigned int, nx.size()>{0})
  {
    int position = _position;
    if (_position == -1)
    {
      rand_seed = std::chrono::system_clock::now().time_since_epoch().count();
      std::srand(rand_seed);
      position = std::rand() % (n_fields_);
    }
    bool success = place_bu(BuildingUnit<nx>(index_bu_max + 1, _jump_parameter, position, _extent));
    return success;
  }
  /*!*********************************************************************************************
   * \brief Creates domain with building units containing only one cell
   *
   * \param _porosity The percentage of void space, not occupied by solid/bu.
   * \param _jump_parameter How far individual particles are allowed to jump.
   * \param random_seed If given, sets random seed to given seed.
   ************************************************************************************************/
  void constexpr place_single_cell_bu_randomly(const double _porosity = 0.90,
                                               const double _jump_parameter = 1,
                                               const unsigned int random_seed = 0)
  {
    if (random_seed == 0)
    {
      rand_seed = std::chrono::system_clock::now().time_since_epoch().count();
    }
    else
      rand_seed = random_seed;

    std::srand(rand_seed);

    unsigned int n_particles = (1. - _porosity) * n_fields_;
    unsigned int position = std::rand() % (n_fields_);
    for (unsigned int i = 0; i < n_particles; ++i)
    {
      while (domain_fields[position] != 0)
        position = std::rand() % (n_fields_);
      building_units.push_back(BuildingUnit<nx>(i + 1, _jump_parameter, position));
      domain_fields[position] = i + 1;
      index_bu_max = i + 1;
    }
  }
  /*!*********************************************************************************************
   * \brief  Finds composites (particles containing more then one bu) in domain and stores
   * information in std::vector<CAM::Composite<nx>*> composites; std::vector<Particle> particles;
   * \deprecated
   ************************************************************************************************/
  void find_composites()
  {
    fields_array_t fields = domain_fields;
    constexpr unsigned int dim = nx.size();
    unsigned int solids_size, field, neigh_field, number;
    std::vector<unsigned int> found_solids;
    particles.clear();
    composites.clear();
    std::for_each(fields.begin(), fields.end(), [](unsigned int& field) { field = (field == 0); });

    for (auto first_solid = std::find(fields.begin(), fields.end(), 0); first_solid != fields.end();
         first_solid = std::find(first_solid, fields.end(), 0))
    {
      std::vector<unsigned int> composite_components;

      found_solids = std::vector<unsigned int>(1, std::distance(fields.begin(), first_solid));
      fields[found_solids[0]] = uint_max;

      composite_components.push_back(domain_fields[found_solids[0]]);
      solids_size = 1;
      for (unsigned int k = 0; k < solids_size; ++k, solids_size = found_solids.size())
      {
        field = found_solids[k];
        for (unsigned int i = 0; i < 2 * dim; ++i)
        {
          neigh_field = aim<nx>(field, direct_neigh<nx>(i));
          if (fields[neigh_field] == 0)
          {
            fields[neigh_field] = uint_max;
            found_solids.push_back(neigh_field);
            number = domain_fields[neigh_field];
            if (std::find(composite_components.begin(), composite_components.end(), number) ==
                composite_components.end())
            {
              composite_components.push_back(number);
            }
          }
        }
      }
      if (composite_components.size() > 1)
      {
        CAM::Composite<nx> new_composite;
        for (unsigned int i = 0; i < composite_components.size(); i++)
        {
          typename std::vector<CAM::BuildingUnit<nx>>::iterator it =
            std::find_if(building_units.begin(), building_units.end(),
                         [&](CAM::BuildingUnit<nx> unit) -> bool
                         { return unit.number == composite_components[i]; });
          new_composite.building_units.push_back(&(*it));
        }
        new_composite.field_indices = found_solids;
        new_composite.jump_parameter =
          CAM::get_jump_range_composite<nx>(new_composite.field_indices.size());
        composites.push_back(new_composite);
      }
      particles.push_back(Particle(found_solids, composite_components));
    }
  }
  /*!*********************************************************************************************
   * \brief  Finds composites (particles containing more then one bu) in domain and stores
   * information in std::vector<CAM::Composite<nx>*> composites; std::vector<Particle> particles;
   * Faster version using the connected particles and their borders
   ************************************************************************************************/
  void find_composites_via_bu_border()
  {
    fields_array_t fields = domain_fields;
    constexpr unsigned int dim = nx.size();
    unsigned int neigh_field, borders_size;
    unsigned int field_number;
    std::vector<unsigned int> borders, found_solids, helper;
    particles.clear();
    composites.clear();
    std::vector<bool> is_bu_visited(index_bu_max, false);

    for (unsigned int i = 0; i < building_units.size(); i++)
    {
      if (is_bu_visited[building_units[i].number] == true)
        continue;
      CAM::Composite<nx> new_composite;
      new_composite.building_units.push_back(&building_units[i]);

      borders.clear();
      borders = building_units[i].get_border_indices();

      found_solids.clear();
      helper = building_units[i].get_field_indices();
      found_solids.insert(found_solids.end(), helper.begin(), helper.end());

      field_number = building_units[i].number;
      is_bu_visited[field_number] = true;
      borders_size = borders.size();

      std::vector<unsigned int> composite_components;
      composite_components.push_back(field_number);

      for (unsigned int j = 0; j < borders_size; j++, borders_size = borders.size())
      {
        for (unsigned int k = 0; k < 2 * dim; ++k)
        {
          neigh_field = aim<nx>(borders[j], direct_neigh<nx>(k));
          field_number = fields[neigh_field];
          if (field_number != 0 && is_bu_visited[field_number] != true)
          {
            is_bu_visited[field_number] = true;

            typename std::vector<CAM::BuildingUnit<nx>>::iterator it = std::find_if(
              building_units.begin(), building_units.end(),
              [&](CAM::BuildingUnit<nx> unit) -> bool { return unit.number == field_number; });

            new_composite.building_units.push_back(&(*it));
            composite_components.push_back(field_number);

            helper = (*it).get_border_indices();
            borders.insert(borders.end(), helper.begin(), helper.end());
            helper = (*it).get_field_indices();
            found_solids.insert(found_solids.end(), helper.begin(), helper.end());
          }
        }
      }
      if (composite_components.size() > 1)
      {
        new_composite.field_indices = found_solids;
        new_composite.jump_parameter =
          CAM::get_jump_range_composite<nx>(new_composite.field_indices.size());
        composites.push_back(new_composite);
      }
      particles.push_back(Particle(found_solids, composite_components));
    }
  }
  /*!***********************************************************************************************
   * \brief   Returns Domain
   *
   * \retval  domain_fields      Information about  status of each cell in domain
   ************************************************************************************************/
  const fields_array_t& fields() const { return domain_fields; }
  /*!***********************************************************************************************
   * \brief   Array of particle locations.
   ************************************************************************************************/
  fields_array_t domain_fields;
  /*!***********************************************************************************************
   * \brief   Vector of particles.
   ************************************************************************************************/
  std::vector<CAM::BuildingUnit<nx>> building_units;
  std::vector<CAM::Composite<nx>> composites;

  /*!*********************************************************************************************
   * \brief vector of particles (connected components)
   * contains
   ************************************************************************************************/
  std::vector<Particle> particles;
  /*!***********************************************************************************************
   * \brief   Random seed.
   ************************************************************************************************/
  unsigned int rand_seed;

  // -------------------------------------------------------------------------------------------------
  // PRINTING SECTION STARTS HERE
  // -------------------------------------------------------------------------------------------------

  /*!*************************************************************************************************
   * \brief   Prints array of domain. Only used in print_n_dim().
   *
   * \param   fields      Particle index and size
   * \param   nx             Nx dimension and size
   * \param   init_index     Temporary index in current slice
   **************************************************************************************************/
  void print_1d(const fields_array_t& _fields, unsigned int init_index = 0)
  {
    const unsigned int len_numbers = std::log10(_fields.size() - 1) + 1;
    unsigned int index;

    for (unsigned int x = 0; x < nx[0]; ++x)
    {
      index = init_index + x;
      if (_fields[index] == 0)
        std::cout << std::setw(len_numbers) << std::setfill('0') << _fields[index] << "  ";
      else
        std::cout << "\033[0;31m" << std::setw(len_numbers) << std::setfill('0') << _fields[index]
                  << "\033[0m  ";
    }
    std::cout << std::endl;
  }
  /*!*************************************************************************************************
   * \brief   Prints array of fields in n dimensions. Only used in print_array().
   *
   * \tparam  n_dim          Temporary dimension of current slice
   * \param   fields      Particle index and size
   * \param   nx             Nx dimension and size
   * \param   init_index     Temporary index in current slice
   **************************************************************************************************/
  template <unsigned int n_dim>
  void print_n_dim(const fields_array_t& _fields, unsigned int init_index = 0)
  {
    if constexpr (n_dim == 1)
      print_1d(_fields, init_index);
    else
    {
      if (n_dim == 2 && nx.size() > 2)
      {
        unsigned int coord_i_helper = init_index;
        std::cout << std::endl;
        for (unsigned int i = 3; i < nx.size() + 1; ++i)
        {
          coord_i_helper = coord_i_helper / nx[i - 2];
          std::cout << i;
          if (i % 10 == 1 && i != 11)
            std::cout << "st";
          else if (i % 10 == 2 && i != 12)
            std::cout << "nd";
          else if (i % 10 == 3 && i != 13)
            std::cout << "rd";
          else
            std::cout << "th";
          std::cout << " coord: " << coord_i_helper % nx[i - 1] << "  ";
        }
        std::cout << std::endl;
      }
      for (unsigned int i = 0; i < nx[n_dim - 1]; ++i)
        print_n_dim<n_dim - 1>(_fields, (init_index + i) * nx[n_dim - 2]);
    }
  }
  /*!*************************************************************************************************
   * \brief   Runs print_n_dim() function
   *
   * \param   fields      Particle index and size
   * \param   nx             Nx dimension and size
   **************************************************************************************************/

  void print_array()
  {
    unsigned int sum = 0;
    std::for_each(domain_fields.begin(), domain_fields.end(), [&](unsigned int t) { sum += t; });
    print_n_dim<nx.size()>(domain_fields);
  }
};
}  // namespace CAM
