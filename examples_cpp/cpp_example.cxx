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
// void print_array(const std::array<unsigned int, size> particles)
// {
//   const unsigned int len_numbers = std::log10(size - 1) + 1;
//   const unsigned int nx = std::sqrt(size);
//   for (unsigned int y = 0; y < nx; ++y)
//   {
//     for (unsigned int x = 0; x < nx; ++x)
//       if (particles[y * nx + x] == 0)
//         std::cout << std::setw(len_numbers) << std::setfill('0') << particles[y * nx + x] << "  ";
//       else
//         std::cout << "\033[0;31m" << std::setw(len_numbers) << std::setfill('0')
//                   << particles[y * nx + x] << "\033[0m"
//                   << "  ";
//     std::cout << std::endl;
//   }
// }
void print_array(const std::array<unsigned int, size> particles, std::array<unsigned int, 2> nx)
{
	const unsigned int len_numbers = std::log10(size - 1) + 1;
	const unsigned int dim = nx.size();
	const unsigned int dim_2_size = nx[0] * nx[1];
	const unsigned int n = std::sqrt(dim_2_size);
	if (dim == 3)
	{
		for (unsigned int z = 0; z < nx[2]; ++z)
		{
      std::cout << "z = " << z << std::endl;
			for (unsigned int y = 0; y < nx[1]; ++y)
  		{
				for (unsigned int x = 0; x < nx[0]; ++x)
				{
			    if (particles[y * n + x] == 0)
			    {
						std::cout << std::setw(len_numbers) << std::setfill('0') << particles[y * n + x] << "  ";
					}
					else
					{
				    std::cout << "\033[0;31m" << std::setw(len_numbers) << std::setfill('0')
				              << particles[y * n + x] << "\033[0m"
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
		    if (particles[y * n + x] == 0)
		  	{
					std::cout << std::setw(len_numbers) << std::setfill('0') << particles[y * n + x] << "  ";
				}
				else
				{
			    std::cout << "\033[0;31m" << std::setw(len_numbers) << std::setfill('0')
			              << particles[y * n + x] << "\033[0m"
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
 *    nx              (the number of rows of the domain)
 *    ny              (the number of columns of the domain)
 *    n_moves         (the number of iterations of the CAM)
 *    porosity        (the percentage of void space, not occupied by solid)
 *    jump_param      (how far individual particles are allowed to jump)
 **************************************************************************************************/
int main()
{
  constexpr std::array<unsigned int, 2> nx = {10, 10};
  const unsigned int n_moves = 5;
  const double porosity = 0.9;
  const double jump_param = 1.;
  cellular_automaton<nx> domain(porosity, jump_param);
  std::cout<< "Seed: " << domain.random_seed() <<std::endl;

  print_array(domain.fields(), nx);

  for (unsigned int i = 0; i < n_moves; ++i)
  {
    std::cout << std::endl;
    print_array(domain.move_particles(), nx);
  }
}
