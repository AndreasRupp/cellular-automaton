#pragma once
#ifndef CELLULAR_AUTOMATON_HXX
#define CELLULAR_AUTOMATON_HXX

#include <CAM/domain.hxx>
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

template <auto nx, typename fields_array_t>// = std::array<unsigned int, n_fields<nx>()>
class CellularAutomaton
{
  static const unsigned int dim = nx.size();

 public:
  // CellularAutomaton();
  // ~CellularAutomaton();
  static void moveBU(CAM::BuildingUnit* _unit, fields_array_t& _domainFields)
  {
    std::vector<int> possible_moves = getStencil(_unit->jump_parameter);
    std::vector<int> best_moves(1, 0);
    double attraction = 0.;
    std::for_each(possible_moves.begin(), possible_moves.end(),
                  [&](int move)
                  {
                    double current_attraction = getAttractionBU(move, _unit, _domainFields);
                    if (current_attraction > attraction)
                    {
                      best_moves = std::vector<int>(1, move);
                      attraction = current_attraction;
                    }
                    else if (current_attraction == attraction)
                      best_moves.push_back(move);
                  });
    doMoveBU(best_moves[std::rand() % best_moves.size()], _unit, _domainFields);
  }
  static void doMoveBU(const int move, CAM::BuildingUnit* _unit, fields_array_t& _domainFields)
  {
    std::vector<unsigned int> fields_old = _unit->getFieldIndices();
    std::for_each(_unit->referencePoints.begin(), _unit->referencePoints.end(),
                  [&](unsigned int& field) { field = aim<nx>(field, move); });
    std::vector<unsigned int> fields_new = _unit->getFieldIndices();
    for (unsigned int i = 0; i < fields_new.size(); i++)
    {
      if (_domainFields[fields_old[i]] == _unit->number)
        _domainFields[fields_old[i]] = 0;
      _domainFields[fields_new[i]] = _unit->number;
    }
  }

  static void apply(Domain<nx, fields_array_t>& _domain)
  {
    _domain.aggregates.clear();
    _domain.findAggregates();
    std::for_each(_domain.aggregates.begin(), _domain.aggregates.end(),
                  [&](CAM::Aggregate* aggregate)
                  { moveAggregate(aggregate, _domain.domainFields); });
    std::shuffle(_domain.buildingUnits.begin(), _domain.buildingUnits.end(),
                 std::default_random_engine(std::rand()));
    std::for_each(_domain.buildingUnits.begin(), _domain.buildingUnits.end(),
                  [&](CAM::BuildingUnit* unit) { moveBU(unit, _domain.domainFields); });
  }

  static void moveAggregate(CAM::Aggregate* _aggregate, fields_array_t& _domainFields)
  {
    std::vector<int> possible_moves = getStencil(_aggregate->jump_parameter);  //
    std::vector<int> best_moves(1, 0);
    double attraction = 0.;
    for (unsigned int i = 0; i < _aggregate->fieldIndices.size(); i++)
    {
      _domainFields[_aggregate->fieldIndices[i]] = 0;
    }
    std::for_each(possible_moves.begin(), possible_moves.end(),
                  [&](int move)
                  {
                    double current_attraction =
                      getAttractionAggregate(move, _aggregate, _domainFields);
                    if (current_attraction > attraction)
                    {
                      best_moves = std::vector<int>(1, move);
                      attraction = current_attraction;
                    }
                    else if (current_attraction == attraction)
                      best_moves.push_back(move);
                  });
    int best_move = best_moves[std::rand() % best_moves.size()];
    for (unsigned int i = 0; i < _aggregate->buildingUnits.size(); i++)
    {
      doMoveBU(best_move, _aggregate->buildingUnits[i], _domainFields);
    }
  }

  static std::vector<int> getStencil(double _jump_parameter)
  {
    unsigned int layers = std::max(1., _jump_parameter);
    std::vector<int> stencil(1, 0);
    unsigned int index = 0;
    unsigned int old_size = stencil.size();

    for (unsigned int lay = 0; lay < layers; ++lay)
    {
      for (; index < old_size; ++index)
        for (unsigned int i = 0; i < 2 * dim; ++i)
          if (std::find(stencil.begin(), stencil.end(), stencil[index] + direct_neigh<nx>(i)) ==
              stencil.end())
            stencil.push_back(stencil[index] + direct_neigh<nx>(i));
      old_size = stencil.size();
    }
    return stencil;
  }
  static double getAttractionBU(const int move,
                                CAM::BuildingUnit* _BU,
                                fields_array_t& _domainFields)
  {
    double attraction = 0.;

    for (unsigned int i = 0; i < _BU->getFieldIndices().size(); ++i)
    {
      unsigned int aiming = aim<nx>(_BU->getFieldIndices()[i], move);
      if (_domainFields[aiming] != _BU->number && _domainFields[aiming] != 0)
        return double_min;
      for (unsigned int i = 0; i < 2 * dim; ++i)
        attraction += _domainFields[aim<nx>(aiming, direct_neigh<nx>(i))] != _BU->number &&
                      _domainFields[aim<nx>(aiming, direct_neigh<nx>(i))] != 0;
    }
    return attraction;
  }

  static double getAttractionAggregate(const int move,
                                       Aggregate* _aggregate,
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
#endif //CELLULAR_AUTOMATON_HXX