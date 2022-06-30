#include <cmath>
#include <iomanip>
#include <iostream>

#include <CAM/cellular_automaton.hxx>

// TODO: Do everything that you did in the cellular_automaton.hxx in this file, too.

template <std::size_t size>
/*!***********************************************************************************************
 * \brief   Prints matrix of particles.
 *
 * \param   particles      Particle number and size
 ************************************************************************************************/
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
/*!***********************************************************************************************
 * \brief   Main function.
 *
 * TODO: Detailed description goes here.
 ************************************************************************************************/
int main()
{
  constexpr unsigned int nx = 10;
  constexpr unsigned int ny = 10;
  const unsigned int n_moves = 3;
  const double porosity = 0.3;
  const double jump_param = 1.;
  cellular_automaton<nx, ny> domain(porosity, jump_param);

  print_array(domain.fields());

  for (unsigned int i = 0; i < n_moves; ++i)
  {
    std::cout << std::endl;
    print_array(domain.move_particles());
  }
}
