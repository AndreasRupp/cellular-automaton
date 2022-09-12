#pragma once

#include <cmath>
#include <iomanip>
#include <iostream>

namespace CAM
{

/*!*************************************************************************************************
 * \brief   Prints array of particles. Only used in print_n_dim().
 *
 * \param   particles      Particle index and size
 * \param   nx             Nx dimension and size
 * \param   init_index     Temporary index in current slice
 **************************************************************************************************/
template <typename fields_array_t, std::size_t dim>
void print_1d(const fields_array_t& particles,
              const std::array<unsigned int, dim>& nx,
              unsigned int init_index = 0)
{
  const unsigned int len_numbers = std::log10(particles.size() - 1) + 1;
  unsigned int index;

  for (unsigned int x = 0; x < nx[0]; ++x)
  {
    index = init_index + x;
    if (particles[index] == 0)
    {
      std::cout << std::setw(len_numbers) << std::setfill('0') << particles[index] << "  ";
    }
    else
    {
      std::cout << "\033[0;31m" << std::setw(len_numbers) << std::setfill('0') << particles[index]
                << "\033[0m"
                << "  ";
    }
  }
  std::cout << std::endl;
}
/*!*************************************************************************************************
 * \brief   Prints array of particles in n dimensions. Only used in print_array().
 *
 * \tparam  n_dim          Temporary dimension of current slice
 * \param   particles      Particle index and size
 * \param   nx             Nx dimension and size
 * \param   init_index     Temporary index in current slice
 **************************************************************************************************/
template <unsigned int n_dim, typename fields_array_t, std::size_t dim>
void print_n_dim(const fields_array_t& particles,
                 const std::array<unsigned int, dim>& nx,
                 unsigned int init_index = 0)
{
  if constexpr (n_dim == 1)
  {
    print_1d(particles, nx, init_index);
  }
  else
  {
    if (n_dim == 2 && dim > 2)
    {
      unsigned int coord_i_helper = init_index;
      std::cout << std::endl;
      for (unsigned int i = 3; i < dim + 1; ++i)
      {
        coord_i_helper = coord_i_helper / nx[i - 2];
        std::cout << i;
        if (i % 10 == 1 && i != 11)
        {
          std::cout << "st";
        }
        else if (i % 10 == 2 && i != 12)
        {
          std::cout << "nd";
        }
        else if (i % 10 == 3 && i != 13)
        {
          std::cout << "rd";
        }
        else
        {
          std::cout << "th";
        }
        std::cout << " coord: " << coord_i_helper % nx[i - 1] << "  ";
      }
      std::cout << std::endl;
    }
    for (unsigned int i = 0; i < nx[n_dim - 1]; ++i)
    {
      unsigned int init_index_loc = (init_index + i) * nx[n_dim - 2];
      print_n_dim<n_dim - 1>(particles, nx, init_index_loc);
    }
  }
}
/*!*************************************************************************************************
 * \brief   Runs print_n_dim() function
 *
 * \param   particles      Particle index and size
 * \param   nx             Nx dimension and size
 **************************************************************************************************/
template <typename fields_array_t, std::size_t dim>
void print_array(const fields_array_t& particles, const std::array<unsigned int, dim>& nx)
{
  print_n_dim<dim>(particles, nx);
}

}  // end of namespace CAM
