#pragma once
#include <CAM/building_units.hxx>
#include <CAM/cellular_automaton.hxx>
#include <CAM/domain.hxx>
#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
namespace CAM
{
template <auto nx, typename fields_array_t = std::array<CAM::fieldNumbers_t, CAM::n_fields<nx>()>>
class CAMInterface
{
  CAM::Domain<nx, fields_array_t> domain;

 public:
  CAMInterface() {}
  ~CAMInterface() {}

  void placeSingleCellBURandomly(double _porosity = 0.5,
                       double _jump_parameter = 1,
                       unsigned int _random_seed = 0)
  {
    domain.placeSingleCellBURandomly(_porosity, _jump_parameter, _random_seed);
  }
  void placeSphere(double _jump_parameter = 1, unsigned int _random_seed = 0)
  {
    domain.placeSphere(_jump_parameter, _random_seed);
  }
  void placeBU()
  {
    unsigned int rand_seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::srand(rand_seed);
    unsigned int centerpoint = 0;
    
    std::vector<int> stencil1 = CAM::ParticleBU<nx>::getStencil(0,10,8);
    std::vector<int> stencil2 = CAM::ParticleBU<nx>::getStencil(0,6,3);
   // CAM::ParticleBU<nx>* particle = new CAM::ParticleBU<nx>(1, 1, centerpoint, stencil);
    for(int a = 0; a<  10; a++)
    {
    unsigned int randomPoint = 0;
    for (unsigned int i = 0; i < nx.size(); ++i)
    {
      centerpoint +=  nx[i]/2 * direct_neigh<nx>(2 * i + 1);
      randomPoint += std::rand() % nx[i] * direct_neigh<nx>(2 * i + 1);
    }
      CAM::ParticleBU<nx>* particle = new CAM::ParticleBU<nx>(a+1, 5, randomPoint, stencil1);
       bool success = domain.placeBU(particle);
       std::cout<<"sucss "<<success<<std::endl;
    }
    
    for(int a= 10; a<  20; a++)
    {
    unsigned int randomPoint = 0;
    for (unsigned int i = 0; i < nx.size(); ++i)
    {
      centerpoint +=  nx[i]/2 * direct_neigh<nx>(2 * i + 1);
      randomPoint += std::rand() % nx[i] * direct_neigh<nx>(2 * i + 1);
    }
      CAM::ParticleBU<nx>* particle = new CAM::ParticleBU<nx>(a+1, 5, randomPoint, stencil2);
       bool success = domain.placeBU(particle);
       std::cout<<"sucss "<<success<<std::endl;
    }
    //domain.placeBU(particle);
  }

  void print_array() { domain.print_array(); }

  void doCAM() { CAM::CellularAutomaton<nx, fields_array_t>::apply(domain); }

  const fields_array_t& fields() const { return domain.domainFields; }
  std::array<double, 12> eval_measures() { return domain.eval_measures(); }
};
}  // namespace CAM