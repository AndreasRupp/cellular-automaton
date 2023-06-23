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
#define EVALUATE_ATTRACTIVITY 1
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
  /*!*********************************************************************************************
   * \brief Finds a new position for an individual building unit
   *
   * \param _unit building unit
   * \param _domain_fields cells of domain
   **********************************************************************************************/
  static void move_bu(CAM::BuildingUnit<nx>& _unit, Domain<nx, fields_array_t>& _domain)
  {
    // if((_unit.get_jump_parameter() % 1) == 0 && )
    //  if(_default_jump_parameter == 0 )
    //  {
    const std::vector<unsigned int> possible_moves =
      CAM::get_stencil<nx>(_unit.get_jump_parameter());
    // }
    // else{
    //   constexpr std::array<unsigned int, get_stencil_size<nx,_default_jump_parameter>()>
    //   possible_moves = CAM::get_stencil_c<nx,_default_jump_parameter()>();
    // }

    // std::cout<<get_stencil_size<nx,1>()<<"nr_cells "<<std::endl;
    // constexpr std::array<unsigned int, get_stencil_size<nx,(unsigned int)_unit.jump_parameter>()>
    // stencil = CAM::get_stencil_c<nx,_unit.get_jump_parameter()>();
    // static_assert(std::cout<<get_stencil_c<nx,5>()[0] == 0);
    //  std::cout<<"start"<<std::endl;
    //  for(unsigned int i = 0; i < get_stencil_c<nx,1>().size(); i++)
    //    std::cout<<get_stencil_c<nx,5>()[i]<<" "<<possible_moves[i]<<std::endl;
    //    std::cout<<"ende"<<std::endl;
    //::vector<unsigned int> best_moves(1, 0);
    std::vector<std::pair<unsigned int, std::array<int, CAM::n_DoF_basis_rotation<nx>()>>>
      best_trans_rot;
    double current_attraction, attraction = 0.;

    std::array<int, CAM::n_DoF_basis_rotation<nx>()> rotation;
    for (unsigned int i = 0; i < n_DoF_basis_rotation<nx>() * 2 + 1; i++)
    {
      std::fill(rotation.begin(), rotation.end(), 0);
      if (i < n_DoF_basis_rotation<nx>() * 2)
      {
        rotation[i / 2] = (i % 2 == 0) ? 1 : -1;
        ;
      }
      BuildingUnit<nx> current_rotation_bu(_unit);
      current_rotation_bu.rotate(rotation);

      std::for_each(possible_moves.begin(), possible_moves.end(),
                    [&](int move)
                    {
                      current_attraction = get_attraction_bu(move, current_rotation_bu, _domain);
                      if (current_attraction > attraction)
                      {
                        best_trans_rot.clear();
                        best_trans_rot.push_back(std::make_pair(move, rotation));
                        attraction = current_attraction;
                      }
                      else if (current_attraction == attraction)
                        best_trans_rot.push_back(std::make_pair(move, rotation));
                    });
    }
    do_move_bu(best_trans_rot[std::rand() % best_trans_rot.size()], _unit, _domain.domain_fields);
  }

  /*!*********************************************************************************************
   * \brief  Executes the actual movement for an individual building unit
   *
   * \param move Index shift induced by possible move.
   * \param _unit building unit
   * \param _domain_fields cells of domain
   **********************************************************************************************/
  static void do_move_bu(
    const std::pair<unsigned int, std::array<int, CAM::n_DoF_basis_rotation<nx>()>>& move,
    CAM::BuildingUnit<nx>& _unit,
    fields_array_t& _domain_fields)
  {
    BuildingUnit bu_old = _unit;
    _unit.rotate(move.second);

    unsigned int field_new, field_old;
    const unsigned int reference_field_old = _unit.get_reference_field();
    _unit.set_reference_field(CAM::aim<nx>(_unit.get_reference_field(), move.first));
    const unsigned int reference_field_new = _unit.get_reference_field();
    for (unsigned int i = 0; i < _unit.get_shape().size(); i++)
    {
      // field_old =
      //   CAM::bu_in_world<nx>(reference_field_old, _unit.get_shape()[i], _unit.get_rotation());
      field_old = CAM::aim<nx>(reference_field_old, bu_old.get_shape()[i]);
      // field_new =
      //   CAM::bu_in_world<nx>(reference_field_new, _unit.get_shape()[i], _unit.get_rotation());
      field_new = CAM::aim<nx>(reference_field_new, _unit.get_shape()[i]);
      _domain_fields[field_new] += _unit.get_number();
      _domain_fields[field_old] -= _unit.get_number();
    }
  }
  /*!*********************************************************************************************
   * \brief Finds a new position for merged building units (composites)
   *
   * \param _unit composites
   * \param _domain_fields cells of domain
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
    std::vector<unsigned int> best_moves(1, 0);
    double current_attraction, attraction = 0.;
    std::for_each(possible_moves.begin(), possible_moves.end(),
                  [&](int move)
                  {
                    current_attraction = get_attraction_composite(move, _composite, _domain);
                    if (current_attraction > attraction)
                    {
                      best_moves = std::vector<unsigned int>(1, move);
                      attraction = current_attraction;
                    }
                    else if (current_attraction == attraction)
                      best_moves.push_back(move);
                  });
    for (unsigned int i = 0; i < _composite.field_indices.size(); i++)
    {
      _domain_fields[_composite.field_indices[i]] = indices[i];
    }

    int best_move = best_moves[std::rand() % best_moves.size()];
    std::array<int, CAM::n_DoF_basis_rotation<nx>()> rotation;
    std::fill(rotation.begin(), rotation.end(), 0);
    for (unsigned int i = 0; i < _composite.building_units.size(); i++)
    {
      do_move_bu(std::make_pair(best_move, rotation), (*_composite.building_units[i]),
                 _domain_fields);
    }
  }
  static double get_attractivity_between_two_faces(const double charge_face1,
                                                   const double charge_face2)
  {
    // return charge_face1 && charge_face2;
    return -(charge_face1 * charge_face2);
  }
  static double get_attraction_bu_(const unsigned int move,
                                   const CAM::BuildingUnit<nx>& _unit,
                                   const fields_array_t& _domain_fields)
  {
    double attraction = 0.;
    unsigned int aiming, field_index;
    for (unsigned int i = 0; i < _unit.get_shape().size(); ++i)
    {
      // field_index = CAM::bu_in_world<nx>(_unit.get_reference_field(), _unit.get_shape()[i],
      //                                    _unit.get_rotation());
      field_index = CAM::aim<nx>(_unit.get_reference_field(), _unit.get_shape()[i]);
      aiming = aim<nx>(field_index, move);
      if (_domain_fields[aiming] != _unit.get_number() && _domain_fields[aiming] != 0)
        return double_min;
      for (unsigned int i = 0; i < 2 * nx.size(); ++i)
        attraction += _domain_fields[aim<nx>(aiming, direct_neigh<nx>(i))] != _unit.get_number() &&
                      _domain_fields[aim<nx>(aiming, direct_neigh<nx>(i))] != 0;
    }
    return attraction;
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
      // field_index = CAM::bu_in_world<nx>(_unit.get_reference_field(), _unit.get_shape()[i],
      //                                    _unit.get_rotation());
      field_index = CAM::aim<nx>(_unit.get_reference_field(), _unit.get_shape()[i]);
      aiming = aim<nx>(field_index, move);
      if (_domain_fields[aiming] != _unit.get_number() &&
          _domain_fields[aiming] != 0)  // occupied by cell
        return double_min;
      // boundary cell
      if (i >= index_begin_bound)
      {
        for (unsigned int j = 0; j < 2 * nx.size(); ++j)
        {
          neigh_index = aim<nx>(aiming, direct_neigh<nx>(j));
          if (_domain_fields[neigh_index] != _unit.get_number() &&
              _domain_fields[neigh_index] != 0)  // boundary cell with neighbor
          {
#if EVALUATE_ATTRACTIVITY == 1

            // neigbor bu
            const CAM::BuildingUnit<nx>& neigh_bu =
              _domain.building_units[_domain.field_number_2_index[_domain_fields[neigh_index]]];
            // faces

            // unsigned int neigh_boundary_cell =
            //   CAM::bu_in_world<nx>(0, CAM::aim<nx>(neigh_index, -neigh_bu.get_reference_field()),
            //                        neigh_bu.get_rotation());
            unsigned int neigh_boundary_cell =
              CAM::aim<nx>(neigh_index, -neigh_bu.get_reference_field());
            // std::cout<<neigh_boundary_cell<<std::endl;
            const std::array<double, nx.size()* 2>& faces_neigh =
              neigh_bu.get_face_charges_of_boundary_cell(neigh_boundary_cell);

            const std::array<double, nx.size()* 2>& faces_unit =
              _unit.get_face_charges()[i - index_begin_bound];

            // opposite face
            unsigned int opposite_face = (j % 2 == 0) ? j + 1 : j - 1;
            // attraction
            attraction +=
              get_attractivity_between_two_faces(faces_neigh[opposite_face], faces_unit[j]);
            // std::cout<<attraction<<std::endl;

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
   * \param   _composite      composites
   * \retval  attraction      Amount of neighbours.
   **********************************************************************************************/
  static double get_attraction_composite(const unsigned int move,
                                         const CAM::Composite<nx>& _composite,
                                         const Domain<nx, fields_array_t>& _domain)
  {
    double attraction_bu, attraction = 0.;
    for (const CAM::BuildingUnit<nx>* bu : _composite.building_units)
    {
      attraction_bu = get_attraction_bu(move, (*bu), _domain);
      if (attraction_bu == double_min)
        return double_min;
      else
        attraction += attraction_bu;
    }
    return attraction;
  }
};
}  // namespace CAM
