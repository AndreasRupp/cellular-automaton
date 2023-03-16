
#include <CAM/cellular_automaton.hxx>
#include <CAM/domain.hxx>
#include <CAM/building_units.hxx>
#include <iostream>
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
  const unsigned int n_moves = 10;
  const double porosity = 0.5;
  const double jump_param = 1.;

  CAM::Domain<nx> domain;
  domain.placeBURandomly(porosity, jump_param);
  domain.print_array();

  for (unsigned int i = 0; i < n_moves; ++i)
  {
    CAM::CellularAutomaton<nx>::apply(domain);
    domain.print_array();
    std::cout << std::endl;
  }
  std::cout << std::endl << "Characteristics / Measures:" << std::endl;
  const std::array<double, 12> meas = domain.eval_measures();
  for (unsigned int k = 0; k < 12; ++k)
    std::cout << "Meas[" << k << "] = " << meas[k] << std::endl;
}
