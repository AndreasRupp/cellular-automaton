//#include <CAM/cellular_automaton.hxx>
#include <CAM/cellular_automaton_new.hxx>
#include <CAM/domain_new.hxx>
#include <CAM/building_units.hxx>
/*!*************************************************************************************************
 * \brief   Main function.
 *
 * Runs CAM, updates and prints the matrix.
 *
 * Parameters:
 *    nx              the size of a row for each dimension of the domain
 *    n_moves         the number of iterations of the CAM
 *    porosity        the percentage of void space, not occupied by solid
 *    jump_param      how far individual particles are allowed to jump
 **************************************************************************************************/
int main()
{
  constexpr std::array<unsigned int, 2> nx = {10, 10};
  //const unsigned int n_moves = 1;
  const double porosity = 0.5;
  const double jump_param = 1.;


  CAM::ParticleBU p;
  //CAM::BuildingUnit p1;
  // p.move(3);
  // p1.move(3);
  p.number = 2;
  std::vector<CAM::BuildingUnit*> buildingUnits;
  buildingUnits.push_back(&p);
  std::for_each(buildingUnits.begin(), buildingUnits.end(), [&](CAM::BuildingUnit* unit) { std::cout<<"dsf"<<unit->number<<std::endl; });
  for(auto unit : buildingUnits) 
  {
        
    std::cout<<unit->number<<std::endl;
    //std::cout<<unit->getFieldIndices().size()<<std::endl;
  }



  CAM::Domain<nx> domain;
  domain.placeBU();
  domain.print_array();

  CAM::CellularAutomaton<nx>::apply(domain);
  domain.print_array();



  // CAM::cellular_automaton<nx> domain(porosity, jump_param);
  // std::cout << "Seed: " << domain.random_seed() << std::endl;

  // CAM::print_array<nx>(domain.fields());

  // for (unsigned int i = 0; i < n_moves; ++i)
  // {
  //   std::cout << std::endl;
  //   CAM::print_array<nx>(domain.move_particles());
  // }

  // std::cout << std::endl << "Characteristics / Measures:" << std::endl;
  // const std::array<double, 12> meas = domain.eval_measures();
  // for (unsigned int k = 0; k < 12; ++k)
  //   std::cout << "Meas[" << k << "] = " << meas[k] << std::endl;
}
