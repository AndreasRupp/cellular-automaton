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
    static const unsigned int dim = nx.size();
    public:
    CellularAutomaton();
    ~CellularAutomaton();
    static void moveAggregate()
    {
      //get Attraction (eigene Funktion)
        //getFieldIndices 
        //in domain auf null stzten
        //dann wieder mit particl number befüllen
      //beim Bewegen jede BU bewegen (sind als pointer drin) (doMoveBU)
    }
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
      //_unit->move<nx>(best_moves[std::rand() % best_moves.size()]);
    }
    static void doMoveBU(const int move, CAM::BuildingUnit* _unit, fields_array_t& _domainFields)
    {
        std::vector<unsigned int> fields_old = _unit->getFieldIndices();
        std::for_each(_unit->referencePoints.begin(), _unit->referencePoints.end(),
                    [&](unsigned int& field)
                    {
                      field = aim<nx>(field, move);
                    });
        std::vector<unsigned int> fields_new = _unit->getFieldIndices();
        for(unsigned int i = 0; i < fields_new.size(); i++)
        {
            // _domainFields[fields_old[i]] -= _unit->number;
            // _domainFields[fields_new[i]] += _unit->number;
            // _domainFields[fields_new[i]] = _domainFields[fields_old[i]];
            // _domainFields[fields_old[i]] = 0;
            _domainFields[fields_old[i]] = 0;
            _domainFields[fields_new[i]] = _unit->number;
            
        }

    }
    static void findAggregates(Domain<nx>& _domain)
    {
      fields_array_t fields = _domain.domainFields;
      constexpr unsigned int dim = nx.size();
      unsigned int solids_size, field, neigh_field;
      std::vector<unsigned int> found_solids, distribution;

      std::for_each(fields.begin(), fields.end(), [](unsigned int& field) { field = (field == 0); });

      for (auto first_solid = std::find(fields.begin(), fields.end(), 0); first_solid != fields.end();
            first_solid = std::find(first_solid, fields.end(), 0))
      {
        std::vector<unsigned int> aggregateComponents;
        
        found_solids = std::vector<unsigned int>(1, std::distance(fields.begin(), first_solid));
        fields[found_solids[0]] = uint_max;

        aggregateComponents.push_back(_domain.domainFields[found_solids[0]]);
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
              unsigned int number = _domain.domainFields[neigh_field];
              if(std::find(aggregateComponents.begin(), aggregateComponents.end(), number) == aggregateComponents.end())
              {
                aggregateComponents.push_back(number);
              }
            }
          }
        }
        if(aggregateComponents.size() > 1)
        {
          
          CAM::Aggregate* newAggregate = new CAM::Aggregate();
          for(unsigned int i = 0; i<aggregateComponents.size();i++)
          {
            std::cout<<aggregateComponents[i] <<" ";

            std::vector<CAM::BuildingUnit*>::iterator it =  std::find_if(_domain.buildingUnits.begin(), _domain.buildingUnits.end(),
                                      [&](CAM::BuildingUnit* unit) -> bool
                                      { return unit->number == aggregateComponents[i]; });
            newAggregate->buildingUnits.push_back(*it);
            newAggregate->fieldIndices = found_solids;
            newAggregate->jump_parameter =  std::max(1., 1 / std::pow(newAggregate->fieldIndices.size(), 1.0 / (double)dim));
          }

          _domain.aggregates.push_back(newAggregate);
        }

        std::cout<<std::endl;
        distribution.push_back(found_solids.size());
      }
      std::sort(distribution.begin(), distribution.end());
      std::cout<<"Distribution"<<std::endl;
      for(unsigned int i = 0; i<distribution.size();i++)
        std::cout<<distribution[i] <<" ";
      std::cout<<std::endl;
    }
    static void apply(Domain<nx>& _domain)
    {
        //std::vector<CAM::BuildingUnit*> buildingUnits = _domain.buildingUnits;
        std::shuffle(_domain.buildingUnits.begin(), _domain.buildingUnits.end(), std::default_random_engine(std::rand()));
        std::for_each(_domain.buildingUnits.begin(), _domain.buildingUnits.end(), [&](CAM::BuildingUnit* unit) { 
        moveBU(unit, _domain.domainFields);});
        std::cout<<"numberofBU" << _domain.buildingUnits.size()<<std::endl;
        _domain.aggregates.clear();
        std::cout<<"numberofBU " << _domain.buildingUnits.size()<<std::endl;
        findAggregates(_domain);
        std::for_each(_domain.aggregates.begin(), _domain.aggregates.end(), [&](CAM::Aggregate* aggregate) { 
        moveAggregate(aggregate, _domain.domainFields);});
        std::cout<<"numberofAgg "<<_domain.aggregates.size()<<" "<<std::endl;
        for(unsigned int i = 0; i < _domain.aggregates.size(); i++)
        {
            std::cout<<"das waren nummern "<<std::endl;;
            std::for_each(_domain.aggregates[i]->buildingUnits.begin(),_domain.aggregates[i]->buildingUnits.end(), [&](CAM::BuildingUnit* unit){
              //std::cout<<unit<<" "<<unit->number<<" "<<unit->getFieldIndices()[0]<< std::endl;
              //  std::cout<<unit->number<<" "<<unit->getFieldIndices()[0]<< std::endl;
             //doMoveBU(1,unit,_domain.domainFields);
              // std::cout<<unit->number<<" "<<unit->getFieldIndices()[0]<< std::endl;
              // std::cout<<std::endl;
              });
          	std::cout<<std::endl;
            //std::cout<<"sizeField "<<_domain.aggregates[i]->fieldIndices.size()<<std::endl;
        }
        std::for_each(_domain.buildingUnits.begin(),_domain.buildingUnits.end(), [&](CAM::BuildingUnit* unit){
             // std::cout<<unit<<" "<<unit->number<<" "<<unit->getFieldIndices()[0]<< std::endl;

              doMoveBU(0,unit,_domain.domainFields);
               });
        

    //     for (unsigned int k = 0; k < particles_size; ++k, particles_size = particles_.size())
    //   particles_[k].update_particle();
    // particles_.erase(
    //   std::remove_if(particles_.begin(), particles_.end(),
    //                  [&](const particle& particle) -> bool { return particle.is_deprecated(); }),
    //   particles_.end());
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
    static void moveAggregate(CAM::Aggregate* _aggregate, fields_array_t& _domainFields)
    {
      std::vector<int> possible_moves = getStencil(_aggregate->jump_parameter);
      std::vector<int> best_moves(1, 0);
      double attraction = 0.;
      for (unsigned int i = 0; i < _aggregate->fieldIndices.size(); i++)
      {
          _domainFields[_aggregate->fieldIndices[i]] = 0;
      }
      std::for_each(possible_moves.begin(), possible_moves.end(),
                    [&](int move)
                    {
                      double current_attraction = getAttractionAggregate(move, _aggregate, _domainFields);
                      if (current_attraction > attraction)
                      {
                        best_moves = std::vector<int>(1, move);
                        attraction = current_attraction;
                      }
                      else if (current_attraction == attraction)
                        best_moves.push_back(move);
                    });
      int best_move = 1;//best_moves[std::rand() % best_moves.size()];
      for(unsigned int i = 0; i < _aggregate->buildingUnits.size(); i++)
      {
          doMoveBU(best_move, _aggregate->buildingUnits[i], _domainFields);
      }
      
      // for (unsigned int i = 0; i < _aggregate->fieldIndices.size(); ++i)
      // {
      //     _domainFields[_aggregate->fieldIndices[i]] = 0;
      // }
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
    static double getAttractionBU(const int move, CAM::BuildingUnit* _BU, fields_array_t& _domainFields)
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

    static double getAttractionAggregate(const int move, Aggregate* _aggregate, fields_array_t& _domainFields)
    {
       double attraction = 0.;
      for (unsigned int i = 0; i < _aggregate->fieldIndices.size(); ++i)
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
}