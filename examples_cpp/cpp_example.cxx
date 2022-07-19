#include <cmath>
#include <iomanip>
#include <iostream>

#include <CAM/cellular_automaton.hxx>

/*!*************************************************************************************************
 * \brief   Prints matrix of particles in 2D and 3D.
 *
 * \param   particles      Particle index and size
 * \param   nx             Nx dimension and size
 **************************************************************************************************/
template <std::size_t size_p, std::size_t dim>
void print_array(const std::array<unsigned int, size_p> particles, std::array<unsigned int, dim> nx)
{
  const unsigned int len_numbers = std::log10(size_p - 1) + 1;
  unsigned int index;
  if (dim == 3)
  {
    for (unsigned int z = 0; z < nx[2]; ++z)
    {
      std::cout << "z = " << z << std::endl;
      for (unsigned int y = 0; y < nx[1]; ++y)
      {
        for (unsigned int x = 0; x < nx[0]; ++x)
        {
          index = (z * nx[1] + y) * nx[0] + x;
          if (particles[index] == 0)
          {
            std::cout << std::setw(len_numbers) << std::setfill('0') << particles[index] << "  ";
          }
          else
          {
            std::cout << "\033[0;31m" << std::setw(len_numbers) << std::setfill('0')
                      << particles[index] << "\033[0m"
                      << "  ";
          }
        }
        std::cout << std::endl;
      }
      std::cout << std::endl;
    }
  }
  else if (dim == 2)
  {
    for (unsigned int y = 0; y < nx[1]; ++y)
    {
      for (unsigned int x = 0; x < nx[0]; ++x)
      {
        if (particles[y * nx[0] + x] == 0)
        {
          std::cout << std::setw(len_numbers) << std::setfill('0') << particles[y * nx[0] + x]
                    << "  ";
        }
        else
        {
          std::cout << "\033[0;31m" << std::setw(len_numbers) << std::setfill('0')
                    << particles[y * nx[0] + x] << "\033[0m"
                    << "  ";
        }
      }
      std::cout << std::endl;
    }
  }
  else
  {
    std::cout << "Printing works only for 2D and 3D." << std::endl;
    return;
  }
}
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
  const unsigned int n_moves = 3;
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
