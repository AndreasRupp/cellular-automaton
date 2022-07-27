#include <CAM/cellular_automaton.hxx>
#include <CAM/print.hxx>

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
  constexpr std::array<unsigned int, 3> nx = {2, 3, 4};
  const unsigned int n_moves = 1;
  const double porosity = 0.9;
  const double jump_param = 1.;
  cellular_automaton<nx> domain(porosity, jump_param);
  std::cout << "Seed: " << domain.random_seed() << std::endl;

  print_array(domain.fields(), nx);

  for (unsigned int i = 0; i < n_moves; ++i)
  {
    std::cout << std::endl;
    print_array(domain.move_particles(), nx);
  }
}
