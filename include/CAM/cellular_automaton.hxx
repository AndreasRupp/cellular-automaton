/**
 * @file cellular_automaton.hxx
 * \brief This class implements the rules of a cellular automaton
 * Can be applied on a domain
 *
 */
#pragma once
#ifndef CELLULAR_AUTOMATON_HXX
#define CELLULAR_AUTOMATON_HXX

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
  /**
   * @brief Applies rules of cellular automaton on domain
   * Generates next generation of CA
   * Single building units and composites are moved
   *
   * @param _domain domain object
   */
  static void apply(Domain<nx, fields_array_t>& _domain)
  {
    _domain.aggregates.clear();
    std::shuffle(_domain.buildingUnits.begin(), _domain.buildingUnits.end(),
                 std::default_random_engine(std::rand()));
    std::cout << "Number of BU: " << _domain.buildingUnits.size() << std::endl;
    std::for_each(_domain.buildingUnits.begin(), _domain.buildingUnits.end(),
                  [&](CAM::BuildingUnit<nx>* unit) { moveBU(unit, _domain.domainFields); });
    _domain.findAggregates();
    std::cout << "Number of aggregates: " << _domain.aggregates.size() << std::endl;
    std::shuffle(_domain.aggregates.begin(), _domain.aggregates.end(),
                 std::default_random_engine(std::rand()));
    std::for_each(_domain.aggregates.begin(), _domain.aggregates.end(),
                  [&](CAM::Aggregate<nx>* aggregate)
                  { moveAggregate(aggregate, _domain.domainFields); });
  }

  /**
   * @brief Finds a new position for an individual building unit
   *
   * @param _unit building unit
   * @param _domainFields cells of domain
   */
  static void moveBU(CAM::BuildingUnit<nx>* _unit, fields_array_t& _domainFields)
  {
    std::vector<unsigned int> possible_moves = CAM::getStencil<nx>(_unit->jump_parameter);
    std::vector<unsigned int> best_moves(1, 0);
    double attraction = 0.;
    std::for_each(possible_moves.begin(), possible_moves.end(),
                  [&](int move)
                  {
                    double current_attraction = getAttractionBU(move, _unit, _domainFields);
                    if (current_attraction > attraction)
                    {
                      best_moves = std::vector<unsigned int>(1, move);
                      attraction = current_attraction;
                    }
                    else if (current_attraction == attraction)
                      best_moves.push_back(move);
                  });
    doMoveBU(best_moves[std::rand() % best_moves.size()], _unit, _domainFields);
  }
  /**
   * @brief  Executes the actual movement for an individual building unit
   *
   * @param move Index shift induced by possible move.
   * @param _unit building unit
   * @param _domainFields cells of domain
   */
  static void doMoveBU(const int move, CAM::BuildingUnit<nx>* _unit, fields_array_t& _domainFields)
  {
    std::vector<unsigned int> fields_old = _unit->getFieldIndices();
    std::for_each(_unit->referencePoints.begin(), _unit->referencePoints.end(),
                  [&](unsigned int& field) { field = aim<nx>(field, move); });
    std::vector<unsigned int> fields_new = _unit->getFieldIndices();
    for (unsigned int i = 0; i < fields_new.size(); i++)
    {
      _domainFields[fields_new[i]] += _unit->number;
      _domainFields[fields_old[i]] -= _unit->number;
    }
  }
  /**
   * @brief Finds a new position for merged building units (composites)
   *
   * @param _unit composites
   * @param _domainFields cells of domain
   */
  static void moveAggregate(CAM::Aggregate<nx>* _aggregate, fields_array_t& _domainFields)
  {
    std::vector<unsigned int> indices(_aggregate->fieldIndices.size(), 0);
    for (unsigned int i = 0; i < _aggregate->fieldIndices.size(); i++)
    {
      indices[i] = _domainFields[_aggregate->fieldIndices[i]];
      _domainFields[_aggregate->fieldIndices[i]] = 0;
    }

    std::vector<unsigned int> possible_moves = CAM::getStencil<nx>(_aggregate->jump_parameter);
    std::vector<unsigned int> best_moves(1, 0);
    double attraction = 0.;
    std::for_each(possible_moves.begin(), possible_moves.end(),
                  [&](int move)
                  {
                    double current_attraction =
                      getAttractionAggregate(move, _aggregate, _domainFields);
                    if (current_attraction > attraction)
                    {
                      best_moves = std::vector<unsigned int>(1, move);
                      attraction = current_attraction;
                    }
                    else if (current_attraction == attraction)
                      best_moves.push_back(move);
                  });
    for (unsigned int i = 0; i < _aggregate->fieldIndices.size(); i++)
    {
      _domainFields[_aggregate->fieldIndices[i]] = indices[i];
    }

    int best_move = best_moves[std::rand() % best_moves.size()];
    for (unsigned int i = 0; i < _aggregate->buildingUnits.size(); i++)
    {
      doMoveBU(best_move, _aggregate->buildingUnits[i], _domainFields);
    }
  }
  /*!*********************************************************************************************
   * \brief   Check attraction of the possible moves for a single BU.
   *
   * \param   move            Index shift induced by possible move.
   * @param _domainFields     cells of domain
   * @param _BU               building unit
   * \retval  attraction      Amount of neighbours.
   **********************************************************************************************/
  static double getAttractionBU(const unsigned int move,
                                CAM::BuildingUnit<nx>* _BU,
                                fields_array_t& _domainFields)
  {
    double attraction = 0.;
    // TODO evaluate Attraction only by borderpoints of BU
    std::vector<unsigned int> fieldIndices = _BU->getFieldIndices();
    for (unsigned int i = 0; i < fieldIndices.size(); ++i)
    {
      unsigned int aiming =
        aim<nx>(fieldIndices[i], move);  // TODO hier könnten man auch einfach stencil addieren
      // aim<nx>(refPoint, move + stencil[i])
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
   * \param   _aggregate      composites
   * \retval  attraction      Amount of neighbours.
   **********************************************************************************************/
  static double getAttractionAggregate(const unsigned int move,
                                       Aggregate<nx>* _aggregate,
                                       fields_array_t& _domainFields)
  {
    double attraction = 0.;
    for (unsigned int i = 0; i < _aggregate->fieldIndices.size(); i++)
    {
      unsigned int aiming = aim<nx>(_aggregate->fieldIndices[i], move);
      if (_domainFields[aiming] != 0 && move != 0)
        return double_min;
      for (unsigned int i = 0; i < 2 * dim; ++i)
        attraction += _domainFields[aim<nx>(aiming, direct_neigh<nx>(i))] != 0;
    }
    return attraction;
  }
};
}  // namespace CAM
#endif  // CELLULAR_AUTOMATON_HXX