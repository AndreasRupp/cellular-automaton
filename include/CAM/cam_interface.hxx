#pragma once 
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <array>
#include <CAM/domain.hxx>
#include <CAM/building_units.hxx>
#include <CAM/cellular_automaton.hxx>
namespace CAM
{
template <auto nx, typename fields_array_t = std::array<unsigned int, CAM::n_fields<nx>()>>
class CAMInterface
{
     CAM::Domain<nx, fields_array_t> domain;
    //CAM::Domain<nx, fields_array_t>& ref = domain;
    public:
   // static constexpr std::array<unsigned int, 2> nx = {10, 10};
//    CAMInterface(double _jump_parameter_composites = 1.0){
//         domain.jump_parameter_composites = _jump_parameter_composites;
//    }
    CAMInterface(){}//domain = new CAM::Domain<nx, fields_array_t>();
   ~CAMInterface(){}//delete domain;

    void placeBURandomly(double _porosity = 0.5, double _jump_parameter = 1, unsigned int _random_seed = 0)
    {domain.placeBURandomly(_porosity, _jump_parameter, _random_seed);}

    void print_array(){domain.print_array();}

    void doCAM(){CAM::CellularAutomaton<nx, fields_array_t>::apply(domain);}
    
    const fields_array_t& fields() const { //(domain.domainFields.begin(), domain.domainFields.end());
        return domain.domainFields; }
    std::array<double, 12> eval_measures(){return domain.eval_measures();}
  
};
}