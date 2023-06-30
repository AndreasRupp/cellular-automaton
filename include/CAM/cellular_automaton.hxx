/*!*********************************************************************************************
 * \file cellular_automaton.hxx
 * \brief This class implements the rules of a cellular automaton
 * Can be applied on a domain
 *TODO calculate attraction only on basis of bu or composites boundary
 **********************************************************************************************/
#pragma once

#include <CAM/domain.hxx>
#include <CAM/utils.hxx>
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <utility>
#include <vector>

#define FACE_ATTRACTIVITY 1
#define ROTATION 1
#define STENCIL_4_ALL_BUS 0
namespace CAM
{

template <auto nx, typename fields_array_t>  //, unsigned int default_jump_parameter
class CellularAutomaton
{
 public:
  /*!*********************************************************************************************
   * \brief Applies rules of cellular automaton on domain
   * Generates next generation of CA
   * Single building units and composites are moved
   *
   * \param _domain domain object
   **********************************************************************************************/
  static void apply(Domain<nx, fields_array_t>& _domain)
  {
    std::shuffle(_domain.building_units.begin(), _domain.building_units.end(),
                 std::default_random_engine(std::rand()));
    // this way easier to handle deprecated bu and random bu numbers.
    _domain.field_number_2_index.resize(_domain.max_field_number + 1, _domain.max_field_number + 1);
    for (unsigned int i = 0; i < _domain.building_units.size(); i++)
    {
      _domain.field_number_2_index[_domain.building_units[i].get_number()] = i;
    }

    std::cout << "Number of bu: " << _domain.building_units.size() << std::endl;
    std::for_each(_domain.building_units.begin(), _domain.building_units.end(),
                  [&](CAM::BuildingUnit<nx>& unit) { move_bu(unit, _domain); });
    std::cout << "move_bu done: " << std::endl;
    _domain.find_composites_via_bu_boundary();
    std::cout << "Number of composites: " << _domain.composites.size() << std::endl;
    std::shuffle(_domain.composites.begin(), _domain.composites.end(),
                 std::default_random_engine(std::rand()));
    std::for_each(_domain.composites.begin(), _domain.composites.end(),
                  [&](CAM::Composite<nx> composite) { move_composites(composite, _domain); });
    std::cout << "move_comp done: " << std::endl;
  }

 private:
  static constexpr std::array<std::array<int, CAM::n_DoF_basis_rotation<nx>()>,
                              CAM::n_DoF_basis_rotation<nx>() * 2 + 1>
  get_90_degree_rotations()
  {
    std::array<std::array<int, CAM::n_DoF_basis_rotation<nx>()>,
               CAM::n_DoF_basis_rotation<nx>() * 2 + 1>
      possible_rotations;
    for (unsigned int i = 0; i < n_DoF_basis_rotation<nx>() * 2 + 1; i++)
    {
      std::fill(possible_rotations[i].begin(), possible_rotations[i].end(), 0);
      if (i < n_DoF_basis_rotation<nx>() * 2)
      {
        possible_rotations[i][i / 2] = (i % 2 == 0) ? 1 : -1;
      }
    }
    return possible_rotations;
  }
  static constexpr std::array<std::array<int, CAM::n_DoF_basis_rotation<nx>()>,
                              CAM::n_DoF_basis_rotation<nx>() + 1>
  get_90_degree_rotations_symmetric()
  {
    std::array<std::array<int, CAM::n_DoF_basis_rotation<nx>()>,
               CAM::n_DoF_basis_rotation<nx>() + 1>
      possible_rotations;
    for (unsigned int i = 0; i < n_DoF_basis_rotation<nx>() + 1; i++)
    {
      std::fill(possible_rotations[i].begin(), possible_rotations[i].end(), 0);
      if (i < n_DoF_basis_rotation<nx>())
      {
        possible_rotations[i][i] = 1;
      }
    }
    return possible_rotations;
  }
  static constexpr std::array<std::array<int, CAM::n_DoF_basis_rotation<nx>()>,
                              CAM::n_DoF_basis_rotation<nx>()* 2 + 1>
    possible_rotations = get_90_degree_rotations();
  static constexpr std::array<int, CAM::n_DoF_basis_rotation<nx>()> no_rotation =
    get_90_degree_rotations()[CAM::n_DoF_basis_rotation<nx>() * 2];

  /*!*********************************************************************************************
   * \brief Finds a new position for an individual building unit
   *
   * \param _unit building unit
   * \param _domain Domain object
   **********************************************************************************************/
  static void move_bu(CAM::BuildingUnit<nx>& _unit, Domain<nx, fields_array_t>& _domain)
  {
#if STENCIL_4_ALL_BUS >= 1
    constexpr std::array<unsigned int, get_stencil_size<nx, STENCIL_4_ALL_BUS>()> possible_moves =
      CAM::get_stencil_c<nx, STENCIL_4_ALL_BUS>();
#else
    const std::vector<unsigned int> possible_moves =
      CAM::get_stencil<nx>(_unit.get_jump_parameter());
#endif
    std::vector<std::pair<unsigned int, std::array<int, CAM::n_DoF_basis_rotation<nx>()>>>
      best_move_rotation;
    best_move_rotation.push_back(std::make_pair(0, no_rotation));
    double current_attraction, attraction = 0.;

#if ROTATION == 1
    for (const std::array<int, CAM::n_DoF_basis_rotation<nx>()>& rotation : possible_rotations)
    {
      BuildingUnit<nx> rotated_unit(_unit);
      rotated_unit.rotate(rotation,
                          CAM::aim<nx>(_unit.get_reference_field(), _unit.get_center_field()));
#else
    std::array<int, CAM::n_DoF_basis_rotation<nx>()> rotation;
    std::fill(rotation.begin(), rotation.end(), 0);
#endif
      std::for_each(possible_moves.begin(), possible_moves.end(),
                    [&](int move)
                    {
#if ROTATION == 1
                      current_attraction = get_attraction_bu(move, rotated_unit, _domain);
#else
current_attraction = get_attraction_bu(move,_unit, _domain);
#endif
                      if (current_attraction > attraction)
                      {
                        best_move_rotation.clear();
                        best_move_rotation.push_back(std::make_pair(move, rotation));
                        attraction = current_attraction;
                      }
                      else if (current_attraction == attraction)
                        best_move_rotation.push_back(std::make_pair(move, rotation));
                    });
#if ROTATION == 1
    }
#endif
    do_move_bu(best_move_rotation[std::rand() % best_move_rotation.size()], _unit,
               _domain.domain_fields);
  }

  /*!*********************************************************************************************
   * \brief  Executes the actual movement for an individual building unit
   *
   * \param {move, rotation} {Index shift induced by possible move, rotation around rotation_point}
   * \param _unit building unit
   * \param _domain_fields cells of domain
   * \param _rotation_point point bu is rotated around
   **********************************************************************************************/
  static void do_move_bu(
    const std::pair<unsigned int, std::array<int, CAM::n_DoF_basis_rotation<nx>()>>& move_rotation,
    CAM::BuildingUnit<nx>& _unit,
    fields_array_t& _domain_fields,
    unsigned int _rotation_point = 0)
  {
    const unsigned int reference_field_old = _unit.get_reference_field();
    const std::vector<unsigned int> shape_old = _unit.get_shape();
#if ROTATION == 1
    if (_rotation_point == 0)
      _rotation_point = CAM::aim<nx>(_unit.get_reference_field(), _unit.get_center_field());
    _unit.rotate(move_rotation.second, _rotation_point);  //
#endif
    unsigned int field_new, field_old;

    _unit.set_reference_field(CAM::aim<nx>(_unit.get_reference_field(), move_rotation.first));
    const unsigned int reference_field_new = _unit.get_reference_field();
    for (unsigned int i = 0; i < _unit.get_shape().size(); i++)
    {
      field_old = CAM::aim<nx>(reference_field_old, shape_old[i]);
      field_new = CAM::aim<nx>(reference_field_new, _unit.get_shape()[i]);
      _domain_fields[field_new] += _unit.get_number();
      _domain_fields[field_old] -= _unit.get_number();
    }
  }
  /*!*********************************************************************************************
   * \brief Finds a new position for merged building units (composites)
   *
   * \param _composite composites
   * \param _domain Domain object
   **********************************************************************************************/
  static void move_composites(CAM::Composite<nx>& _composite, Domain<nx, fields_array_t>& _domain)
  {
    fields_array_t& _domain_fields = _domain.domain_fields;
    std::vector<unsigned int> indices(_composite.field_indices.size(), 0);
    for (unsigned int i = 0; i < _composite.field_indices.size(); i++)
    {
      indices[i] = _domain_fields[_composite.field_indices[i]];
      _domain_fields[_composite.field_indices[i]] = 0;
    }

    const std::vector<unsigned int> possible_moves =
      CAM::get_stencil<nx>(_composite.jump_parameter);

    std::vector<std::pair<unsigned int, std::array<int, CAM::n_DoF_basis_rotation<nx>()>>>
      best_move_rotation;
    best_move_rotation.push_back(std::make_pair(0, no_rotation));

    double current_attraction, attraction = 0.;
    std::vector<CAM::BuildingUnit<nx>> bus;
    for (const CAM::BuildingUnit<nx>* bu : _composite.building_units)
      bus.push_back(*bu);

#if ROTATION == 1
    unsigned int rotation_center = CAM::get_center_field<nx>(_composite.field_indices);
    for (const std::array<int, CAM::n_DoF_basis_rotation<nx>()>& rotation : possible_rotations)
    {
      for (unsigned int i = 0; i < _composite.building_units.size(); i++)
      {
        BuildingUnit<nx> rotated_unit(*(_composite.building_units[i]));
        rotated_unit.rotate(rotation, rotation_center);
        bus[i] = rotated_unit;
      }

#else
    std::array<int, CAM::n_DoF_basis_rotation<nx>()> rotation;
    std::fill(rotation.begin(), rotation.end(), 0);
#endif
      std::for_each(possible_moves.begin(), possible_moves.end(),
                    [&](int move)
                    {
                      current_attraction = get_attraction_composite(move, bus, _domain);

                      if (current_attraction > attraction)
                      {
                        best_move_rotation.clear();
                        best_move_rotation.push_back(std::make_pair(move, rotation));
                        attraction = current_attraction;
                      }
                      else if (current_attraction == attraction)
                        best_move_rotation.push_back(std::make_pair(move, rotation));
                    });
#if ROTATION == 1
    }
#endif
    for (unsigned int i = 0; i < _composite.field_indices.size(); i++)
    {
      _domain_fields[_composite.field_indices[i]] = indices[i];
    }

    std::pair<unsigned int, std::array<int, CAM::n_DoF_basis_rotation<nx>()>> chosen_move_rotation =
      best_move_rotation[std::rand() % best_move_rotation.size()];

    for (unsigned int i = 0; i < _composite.building_units.size(); i++)
    {
#if ROTATION == 1
      do_move_bu(chosen_move_rotation, *(_composite.building_units[i]), _domain_fields,
                 rotation_center);
#else
      do_move_bu(chosen_move_rotation, *(_composite.building_units[i]), _domain_fields);
#endif
    }
  }
  /*!*********************************************************************************************
   * \brief Definition of face attractivity
   **********************************************************************************************/
  static double get_attractivity_between_two_faces(const double charge_face1,
                                                   const double charge_face2)
  {
    // return charge_face1 && charge_face2;
    return -(charge_face1 * charge_face2);
  }
  /*!*********************************************************************************************
   * \brief   Check attraction of the possible moves for a single bu.
   *
   * \param   move            Index shift induced by possible move.
   * \param _domain_fields     cells of domain
   * \param _unit               building unit
   * \retval  attraction      Amount of neighbours.
   **********************************************************************************************/
  static double get_attraction_bu(const unsigned int move,
                                  const CAM::BuildingUnit<nx>& _unit,
                                  const Domain<nx, fields_array_t>& _domain)
  {
    double attraction = 0.;
    const fields_array_t& _domain_fields = _domain.domain_fields;
    unsigned int aiming, neigh_index, field_index;
    unsigned int index_begin_bound = (_unit.get_shape().size() - _unit.get_boundary().size());

    for (unsigned int i = 0; i < _unit.get_shape().size(); ++i)
    {
      field_index = CAM::aim<nx>(_unit.get_reference_field(), _unit.get_shape()[i]);
      aiming = aim<nx>(field_index, move);
      if (_domain_fields[aiming] != _unit.get_number() &&
          _domain_fields[aiming] != 0)  // occupied by cell
        return double_min;
      // boundary cell
      if (i >= index_begin_bound)
      {
        for (unsigned int j = 0; j < 2 * nx.size(); j++)
        {
          neigh_index = aim<nx>(aiming, direct_neigh<nx>(j));
          if (_domain_fields[neigh_index] != _unit.get_number() &&
              _domain_fields[neigh_index] != 0)  // boundary cell with neighbor
          {
#if FACE_ATTRACTIVITY == 1

            // neigbor bu
            const CAM::BuildingUnit<nx>& neigh_bu =
              _domain.building_units[_domain.field_number_2_index[_domain_fields[neigh_index]]];
            // faces
            unsigned int neigh_boundary_cell =
              CAM::aim<nx>(neigh_index, -neigh_bu.get_reference_field());

            const std::array<double, nx.size()* 2>& faces_neigh =
              neigh_bu.get_face_charges_of_boundary_cell(neigh_boundary_cell);

            const std::array<double, nx.size()* 2>& faces_unit =
              _unit.get_face_charges()[i - index_begin_bound];

            // opposite face
            unsigned int opposite_face = (j % 2 == 0) ? j + 1 : j - 1;
            // attraction

            attraction +=
              get_attractivity_between_two_faces(faces_neigh[opposite_face], faces_unit[j]);
#else
            attraction += 1;
#endif
          }
        }
      }
    }
    return attraction;
  }
  /*!*********************************************************************************************
   * \brief   Check attraction of the possible moves for merged bus.
   *
   * \param   move            Index shift induced by possible move.
   * \param   _bus            bu of composite
   * \param _domain           Domain object
   * \retval  attraction      Amount of neighbours.
   **********************************************************************************************/
  static double get_attraction_composite(const unsigned int move,
                                         const std::vector<CAM::BuildingUnit<nx>>& bus,
                                         const Domain<nx, fields_array_t>& _domain)
  {
    double attraction_bu, attraction = 0.;
    for (const CAM::BuildingUnit<nx>& bu : bus)
    {
      attraction_bu = get_attraction_bu(move, (bu), _domain);
      if (attraction_bu == double_min)
        return double_min;
      else
        attraction += attraction_bu;
    }
    return attraction;
  }
};
}  // namespace CAM
