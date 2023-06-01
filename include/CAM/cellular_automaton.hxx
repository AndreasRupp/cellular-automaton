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
#include <vector>

namespace CAM
{

template <auto nx, typename fields_array_t>
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
    //std::cout << "Number of bu: " << _domain.building_units.size() << std::endl;
    std::for_each(_domain.building_units.begin(), _domain.building_units.end(),
                  [&](CAM::BuildingUnit<nx>& unit) { move_bu(unit, _domain.domain_fields); });
    _domain.find_composites_via_bu_boundary();
    //std::cout << "Number of composites: " << _domain.composites.size() << std::endl;
    std::shuffle(_domain.composites.begin(), _domain.composites.end(),
                 std::default_random_engine(std::rand()));
    std::for_each(_domain.composites.begin(), _domain.composites.end(),
                  [&](CAM::Composite<nx> composite)
                  { move_composites(composite, _domain.domain_fields); });
  }
  private:
  static const unsigned int dim = nx.size();
  /*!*********************************************************************************************
   * \brief Finds a new position for an individual building unit
   *
   * \param _unit building unit
   * \param _domain_fields cells of domain
   **********************************************************************************************/
  static void move_bu(CAM::BuildingUnit<nx>& _unit, fields_array_t& _domain_fields)
  {
    const std::vector<unsigned int> possible_moves =
      CAM::get_stencil<nx>(_unit.get_jump_parameter());
    std::vector<unsigned int> best_moves(1, 0);
    double current_attraction, attraction = 0.;
    std::for_each(possible_moves.begin(), possible_moves.end(),
                  [&](int move)
                  {
                    current_attraction = get_attraction_bu(move, _unit, _domain_fields);
                    if (current_attraction > attraction)
                    {
                      best_moves = std::vector<unsigned int>(1, move);
                      attraction = current_attraction;
                    }
                    else if (current_attraction == attraction)
                      best_moves.push_back(move);
                  });
    do_move_bu(best_moves[std::rand() % best_moves.size()], _unit, _domain_fields);
  }
  /*!*********************************************************************************************
   * \brief  Executes the actual movement for an individual building unit
   *
   * \param move Index shift induced by possible move.
   * \param _unit building unit
   * \param _domain_fields cells of domain
   **********************************************************************************************/
  static void do_move_bu(const int move,
                         CAM::BuildingUnit<nx>& _unit,
                         fields_array_t& _domain_fields)
  {
    unsigned int field_new, field_old;
    const unsigned int reference_field_old = _unit.get_reference_field();
    _unit.set_reference_field(CAM::aim<nx>(_unit.get_reference_field(), move));
    const unsigned int reference_field_new = _unit.get_reference_field();
    for (unsigned int i = 0; i < _unit.get_shape().size(); i++)
    {
      field_old = CAM::aim<nx>(reference_field_old, _unit.get_shape()[i]);
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
  static void move_composites(CAM::Composite<nx>& _composite, fields_array_t& _domain_fields)
  {
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
                    current_attraction = get_attraction_composite(move, _composite, _domain_fields);
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
    for (unsigned int i = 0; i < _composite.building_units.size(); i++)
    {
      do_move_bu(best_move, (*_composite.building_units[i]), _domain_fields);
    }
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
                                  CAM::BuildingUnit<nx>& _unit,
                                  const fields_array_t& _domain_fields)
  {
    double attraction = 0.;
    unsigned int aiming, field_index;
    for (unsigned int i = 0; i < _unit.get_shape().size(); ++i)
    {
      field_index = CAM::aim<nx>(_unit.get_reference_field(), _unit.get_shape()[i]);
      aiming = aim<nx>(field_index, move);
      if (_domain_fields[aiming] != _unit.get_number() && _domain_fields[aiming] != 0)
        return double_min;
      for (unsigned int i = 0; i < 2 * dim; ++i)
        attraction += _domain_fields[aim<nx>(aiming, direct_neigh<nx>(i))] != _unit.get_number() &&
                      _domain_fields[aim<nx>(aiming, direct_neigh<nx>(i))] != 0;
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
                                         const fields_array_t& _domain_fields)
  {
    double attraction = 0.;
    unsigned int aiming;
    for (unsigned int i = 0; i < _composite.field_indices.size(); i++)
    {
      aiming = aim<nx>(_composite.field_indices[i], move);
      if (_domain_fields[aiming] != 0 && move != 0)
        return double_min;
      for (unsigned int i = 0; i < 2 * dim; ++i)
        attraction += _domain_fields[aim<nx>(aiming, direct_neigh<nx>(i))] != 0;
    }
    return attraction;
  }
};
}  // namespace CAM
