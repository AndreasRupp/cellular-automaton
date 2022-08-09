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
  constexpr std::array<unsigned int, 2> nx = {5, 5};
  const unsigned int n_moves = 1;
  const double porosity = 0.5;
  const double jump_param = 1.;
  cellular_automaton<nx> domain(porosity, jump_param, 123456789);
  std::cout << "Seed: " << domain.random_seed() << std::endl;

  print_array(domain.fields(), nx);

  for (unsigned int i = 0; i < n_moves; ++i)
  {
    std::cout << std::endl;
    print_array(domain.move_particles(), nx);
  }
  const std::array<double, 12> meas = domain.eval_measures();
  for (unsigned int k = 0; k < 12; ++k)
  {
    std::cout << "Meas[" << k << "] = " << meas[k] << std::endl;
  }
  auto helper = domain.fields();
  for (unsigned int i = 0; i < helper.size(); ++i)
    helper[i] = i;
  print_array(helper, nx);
  std::cout << (domain.fields())[9];
}
