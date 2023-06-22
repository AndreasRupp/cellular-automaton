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
 public:
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
      domain_fields.resize(n_fields<nx>(), 0);
    else
    {
      static_assert(
        std::is_same<fields_array_t,
                     std::array<typename fields_array_t::value_type, n_fields<nx>()>>::value,
        "The fields array has incorrect size");
      domain_fields.fill(0);
    }
    max_field_number = 0;
  }
  /*!*********************************************************************************************
   * \brief Place any building unit (bu) into the domain
   * Basis function for placement for all special BUs
   * \param _unit building unit
   * \return true: cells of bu are marked with individual index, bu is stored in bu vector
   * \return false: bu could not be placed. all its cells are removed
   ************************************************************************************************/
  bool constexpr place_bu(const CAM::BuildingUnit<nx>& _unit)
  {
    // check if the number is already in the domain
    // auto iter = std::find_if(building_units.begin(), building_units.end(),
    //                          [&](const CAM::BuildingUnit<nx>& bu)
    //                          { return bu.get_number() == _unit.get_number(); });
    // if (_unit.get_number() == 0 || iter != building_units.end())
    //   return false;
    if (_unit.get_number() > max_field_number)
      max_field_number = _unit.get_number();

    unsigned int field;
    for (unsigned int i = 0; i < _unit.get_shape().size(); i++)
    {
      field = CAM::bu_in_world<nx>(_unit.get_reference_field(), _unit.get_shape()[i],
                                   _unit.get_rotation());
      if (domain_fields[field] != 0)
        return false;
    }

    for (unsigned int i = 0; i < _unit.get_shape().size(); i++)
    {
      field = CAM::bu_in_world<nx>(_unit.get_reference_field(), _unit.get_shape()[i],
                                   _unit.get_rotation());
      domain_fields[field] = _unit.get_number();
    }
    building_units.push_back(_unit);
    return true;
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
                         { return unit.get_number() == composite_components[i]; });
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
   * Faster version using the connected particles and their boundaries
   ************************************************************************************************/
  void find_composites_via_bu_boundary()
  {
    fields_array_t fields = domain_fields;
    constexpr unsigned int dim = nx.size();
    unsigned int neigh_field, boundaries_size;
    unsigned int field_number;
    std::vector<unsigned int> boundaries, found_solids, helper;
    particles.clear();
    composites.clear();
    std::vector<bool> is_bu_visited(building_units.size(), false);

    for (unsigned int i = 0; i < building_units.size(); i++)
    {
      if (is_bu_visited[building_units[i].get_number()] == true)
        continue;
      CAM::Composite<nx> new_composite;
      new_composite.building_units.push_back(&building_units[i]);

      boundaries.clear();
      for (unsigned int boundary_field : building_units[i].get_boundary())
        boundaries.push_back(CAM::bu_in_world<nx>(building_units[i].get_reference_field(),
                                                  boundary_field,
                                                  building_units[i].get_rotation()));

      found_solids.clear();
      for (unsigned int shape_field : building_units[i].get_shape())
        found_solids.push_back(CAM::bu_in_world<nx>(building_units[i].get_reference_field(),
                                                    shape_field, building_units[i].get_rotation()));

      field_number = building_units[i].get_number();
      is_bu_visited[field_number] = true;
      boundaries_size = boundaries.size();

      std::vector<unsigned int> composite_components;
      composite_components.push_back(field_number);

      for (unsigned int j = 0; j < boundaries_size; j++, boundaries_size = boundaries.size())
      {
        for (unsigned int k = 0; k < 2 * dim; ++k)
        {
          neigh_field = aim<nx>(boundaries[j], direct_neigh<nx>(k));
          field_number = fields[neigh_field];
          if (field_number != 0 && is_bu_visited[field_number] != true)
          {
            is_bu_visited[field_number] = true;

            unsigned int index = field_number_2_index[field_number];

            new_composite.building_units.push_back(&(building_units[index]));
            composite_components.push_back(field_number);

            for (unsigned int boundary_field : (building_units[index]).get_boundary())
              boundaries.push_back(
                CAM::aim<nx>((building_units[index]).get_reference_field(), boundary_field));

            for (unsigned int shape_field : (building_units[index]).get_shape())
              found_solids.push_back(
                CAM::bu_in_world<nx>((building_units[index]).get_reference_field(), shape_field,
                                     building_units[i].get_rotation()));
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

  unsigned int max_field_number;
  std::vector<unsigned int> field_number_2_index;

  /*!*********************************************************************************************
   * \brief vector of particles (connected components)
   * contains
   ************************************************************************************************/
  std::vector<Particle> particles;

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
