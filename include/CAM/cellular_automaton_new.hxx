#pragma once

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <limits>
#include <random>
#include <vector>
#include<iostream>
#include <CAM/domain_new.hxx>

namespace CAM
{
template <auto nx, typename fields_array_t = std::array<unsigned int, n_fields<nx>()>>
class CellularAutomaton
{
    static constexpr unsigned int uint_max = std::numeric_limits<unsigned int>::max();
    /*!***********************************************************************************************
    * \brief   Smallest (negative) double.
    ************************************************************************************************/
    static constexpr double double_min = std::numeric_limits<double>::lowest();
    static const unsigned int dim = nx.size();
    public:
    CellularAutomaton();
    ~CellularAutomaton();
    static void moveBU(CAM::BuildingUnit* _unit, fields_array_t& _domainFields)
    {
      std::vector<int> possible_moves = getStencil(4);
      std::vector<int> best_moves(1, 0);
      double attraction = 0.;
      std::for_each(possible_moves.begin(), possible_moves.end(),
                    [&](int move)
                    {
                      double current_attraction = getAttraction(move, _unit, _domainFields);
                      if (current_attraction > attraction)
                      {
                        best_moves = std::vector<int>(1, move);
                        attraction = current_attraction;
                      }
                      else if (current_attraction == attraction)
                        best_moves.push_back(move);
                    });
     // doMove(best_moves[std::rand() % best_moves.size()], _unit, _domainFields);
      //_unit->move<nx>(best_moves[std::rand() % best_moves.size()]);
    }
    // void doMove(const int move, CAM::BuildingUnit* _unit, fields_array_t& _domainFields)
    // {
    //     std::vector<unsigned int> fields = _unit->getFieldIndices();
    //     std::for_each(fields_.begin(), fields_.end(),
    //                 [&](unsigned int& field)
    //                 {
    //                   _domain_fields[field] -= number_;
    //                   field = aim<nx>(field, move);
    //                   _domain_fields[field] += number_;
    //                 });
    // }
    static void apply(Domain<nx>& _domain)
    {
        //std::vector<CAM::BuildingUnit*> buildingUnits = _domain.buildingUnits;
        std::shuffle(_domain.buildingUnits.begin(), _domain.buildingUnits.end(), std::default_random_engine(std::rand()));
        std::for_each(_domain.buildingUnits.begin(), _domain.buildingUnits.end(), [&](CAM::BuildingUnit* unit) { 
        moveBU(unit, _domain.domainFields);
        std::cout<<"dsf"<<unit->number<<std::endl; });
    //    for( unit : buildingUnits) 
    //    {
        
    //     std::cout<<unit->_number<<std::endl;
    //    }
        
        /* std::shuffle(particles_.begin(), particles_.end(), std::default_random_engine(std::rand()));
        std::for_each(particles_.begin(), particles_.end(),
                  [&](particle& part) { part.move_singles(); });
        unsigned int particles_size = particles_.size();
        for (unsigned int k = 0; k < particles_size; ++k, particles_size = particles_.size())
         particles_[k].update_particle();
        particles_.erase(
        std::remove_if(particles_.begin(), particles_.end(),
                     [&](const particle& particle) -> bool { return particle.is_deprecated(); }),
        particles_.end()); */

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
    static double getAttraction(const int move, CAM::BuildingUnit* _BU, fields_array_t& _domainFields)
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
};
}