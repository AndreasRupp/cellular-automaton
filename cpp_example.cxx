#include <cmath>
#include <iomanip>
#include <iostream>

#include <CAM/cellular_automaton.hxx>

template <std::size_t size>
/*!*************************************************************************************************
 * \brief   Prints matrix of particles.
 *
 * \param   particles      Particle number and size
 **************************************************************************************************/
void print_array(const std::array<unsigned int, size> particles)
{
  const unsigned int len_numbers = std::log10(size - 1) + 1;
  const unsigned int nx = std::sqrt(size);
  for (unsigned int y = 0; y < nx; ++y)
  {
    for (unsigned int x = 0; x < nx; ++x)
      if (particles[y * nx + x] == 0)
        std::cout << std::setw(len_numbers) << std::setfill('0') << particles[y * nx + x] << "  ";
      else
        std::cout << "\033[0;31m" << std::setw(len_numbers) << std::setfill('0')
                  << particles[y * nx + x] << "\033[0m"
                  << "  ";
    std::cout << std::endl;
  }
}
/*!*************************************************************************************************
 * \brief   Main function.
 *
 * Runs CAM, updates and prints the matrix.
 *
 * Parameters:
 *    nx              (the number of rows of the domain)
 *    ny              (the number of columns of the domain)
 *    n_moves         (the number of iterations of the CAM)
 *    porosity        (the percentage of void space, not occupied by solid)
 *    jump_param      (how far individual particles are allowed to jump)
 **************************************************************************************************/
int main()
{
  constexpr std::array<unsigned int, 2> nx = {10, 10};
  const unsigned int n_moves = 3;
  const double porosity = 0.3;
  const double jump_param = 1.;
  cellular_automaton<decltype(nx), nx> domain(porosity, jump_param);

  print_array(domain.fields());

  for (unsigned int i = 0; i < n_moves; ++i)
  {
    std::cout << std::endl;
    print_array(domain.move_particles());
  }
}
