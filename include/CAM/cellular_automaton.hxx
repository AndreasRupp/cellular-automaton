/*!*********************************************************************************************
 * \file cellular_automaton.hxx
 * \brief This class implements the rules of a cellular automaton
 * Can be applied on a domain
 *TODO calculate attraction only on basis of BU or composites border
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
  static const unsigned int dim = nx.size();

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
    std::shuffle(_domain.buildingUnits.begin(), _domain.buildingUnits.end(),
                 std::default_random_engine(std::rand()));
    // std::cout << "Number of BU: " << _domain.buildingUnits.size() << std::endl;
    std::for_each(_domain.buildingUnits.begin(), _domain.buildingUnits.end(),
                  [&](CAM::BuildingUnit<nx>* unit) { move_BU(unit, _domain.domainFields); });
    _domain.find_composites_via_BU_border();
    // std::cout << "Number of composites: " << _domain.composites.size() << std::endl;
    // _domain.find_composites();
    std::shuffle(_domain.composites.begin(), _domain.composites.end(),
                 std::default_random_engine(std::rand()));
    std::for_each(_domain.composites.begin(), _domain.composites.end(),
                  [&](CAM::Composite<nx>* composite)
                  { move_composites(composite, _domain.domainFields); });
  }
  /*!*********************************************************************************************
   * \brief Finds a new position for an individual building unit
   *
   * \param _unit building unit
   * \param _domainFields cells of domain
   **********************************************************************************************/
  static void move_BU(CAM::BuildingUnit<nx>* _unit, fields_array_t& _domainFields)
  {
    std::vector<unsigned int> possible_moves = CAM::get_stencil<nx>(_unit->jump_parameter);
    std::vector<unsigned int> best_moves(1, 0);
    double current_attraction, attraction = 0.;
    std::for_each(possible_moves.begin(), possible_moves.end(),
                  [&](int move)
                  {
                    current_attraction = get_attraction_BU(move, _unit, _domainFields);
                    if (current_attraction > attraction)
                    {
                      best_moves = std::vector<unsigned int>(1, move);
                      attraction = current_attraction;
                    }
                    else if (current_attraction == attraction)
                      best_moves.push_back(move);
                  });
    do_move_BU(best_moves[std::rand() % best_moves.size()], _unit, _domainFields);
  }
  /*!*********************************************************************************************
   * \brief  Executes the actual movement for an individual building unit
   *
   * \param move Index shift induced by possible move.
   * \param _unit building unit
   * \param _domainFields cells of domain
   **********************************************************************************************/
  static void do_move_BU(const int move,
                         CAM::BuildingUnit<nx>* _unit,
                         fields_array_t& _domainFields)
  {
    std::vector<unsigned int> fields_old = _unit->get_field_indices();
    std::for_each(_unit->referenceFields.begin(), _unit->referenceFields.end(),
                  [&](unsigned int& field) { field = aim<nx>(field, move); });
    std::vector<unsigned int> fields_new = _unit->get_field_indices();
    for (unsigned int i = 0; i < fields_new.size(); i++)
    {
      _domainFields[fields_new[i]] += _unit->number;
      _domainFields[fields_old[i]] -= _unit->number;
    }
  }
  /*!*********************************************************************************************
   * \brief Finds a new position for merged building units (composites)
   *
   * \param _unit composites
   * \param _domainFields cells of domain
   **********************************************************************************************/
  static void move_composites(CAM::Composite<nx>* _composite, fields_array_t& _domainFields)
  {
    std::vector<unsigned int> indices(_composite->fieldIndices.size(), 0);
    for (unsigned int i = 0; i < _composite->fieldIndices.size(); i++)
    {
      indices[i] = _domainFields[_composite->fieldIndices[i]];
      _domainFields[_composite->fieldIndices[i]] = 0;
    }

    std::vector<unsigned int> possible_moves = CAM::get_stencil<nx>(_composite->jump_parameter);
    std::vector<unsigned int> best_moves(1, 0);
    double current_attraction, attraction = 0.;
    std::for_each(possible_moves.begin(), possible_moves.end(),
                  [&](int move)
                  {
                    current_attraction = get_attraction_composite(move, _composite, _domainFields);
                    if (current_attraction > attraction)
                    {
                      best_moves = std::vector<unsigned int>(1, move);
                      attraction = current_attraction;
                    }
                    else if (current_attraction == attraction)
                      best_moves.push_back(move);
                  });
    for (unsigned int i = 0; i < _composite->fieldIndices.size(); i++)
    {
      _domainFields[_composite->fieldIndices[i]] = indices[i];
    }

    int best_move = best_moves[std::rand() % best_moves.size()];
    for (unsigned int i = 0; i < _composite->buildingUnits.size(); i++)
    {
      do_move_BU(best_move, _composite->buildingUnits[i], _domainFields);
    }
  }
  /*!*********************************************************************************************
   * \brief   Check attraction of the possible moves for a single BU.
   *
   * \param   move            Index shift induced by possible move.
   * \param _domainFields     cells of domain
   * \param _BU               building unit
   * \retval  attraction      Amount of neighbours.
   **********************************************************************************************/
  static double get_attraction_BU(const unsigned int move,
                                  CAM::BuildingUnit<nx>* _BU,
                                  fields_array_t& _domainFields)
  {
    double attraction = 0.;
    unsigned int aiming;
    // TODO evaluate Attraction only by borderpoints of BU
    std::vector<unsigned int> fieldIndices = _BU->get_field_indices();
    for (unsigned int i = 0; i < fieldIndices.size(); ++i)
    {
      aiming = aim<nx>(fieldIndices[i], move);
      if (_domainFields[aiming] != _BU->number && _domainFields[aiming] != 0)
        return double_min;
      for (unsigned int i = 0; i < 2 * dim; ++i)
        attraction += _domainFields[aim<nx>(aiming, direct_neigh<nx>(i))] != _BU->number &&
                      _domainFields[aim<nx>(aiming, direct_neigh<nx>(i))] != 0;
    }
    return attraction;
  }
  /*!*********************************************************************************************
   * \brief   Check attraction of the possible moves for merged BUs.
   *
   * \param   move            Index shift induced by possible move.
   * \param   _composite      composites
   * \retval  attraction      Amount of neighbours.
   **********************************************************************************************/
  static double get_attraction_composite(const unsigned int move,
                                         Composite<nx>* _composite,
                                         fields_array_t& _domainFields)
  {
    double attraction = 0.;
    unsigned int aiming;
    for (unsigned int i = 0; i < _composite->fieldIndices.size(); i++)
    {
      aiming = aim<nx>(_composite->fieldIndices[i], move);
      if (_domainFields[aiming] != 0 && move != 0)
        return double_min;
      for (unsigned int i = 0; i < 2 * dim; ++i)
        attraction += _domainFields[aim<nx>(aiming, direct_neigh<nx>(i))] != 0;
    }
    return attraction;
  }
};
}  // namespace CAM
