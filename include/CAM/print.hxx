#pragma once

#include <cmath>
#include <iomanip>
#include <iostream>

template <std::size_t size_p, std::size_t dim>
void print_1d(const std::array<unsigned int, size_p>& particles,
              const std::array<unsigned int, dim>& nx,
              unsigned int init_index = 0)
{
  const unsigned int len_numbers = std::log10(size_p - 1) + 1;
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

template <unsigned int n_dim, std::size_t size_p, std::size_t dim>
void print_n_dim(const std::array<unsigned int, size_p>& particles,
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
      for (unsigned int i = dim; i > 2; --i)
      {
        coord_i_helper = coord_i_helper / nx[i - 2];
        std::cout << i << "th coord: " << coord_i_helper % nx[i - 1] << "  ";
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

template <std::size_t size_p, std::size_t dim>
void print_array(const std::array<unsigned int, size_p>& particles,
                 const std::array<unsigned int, dim>& nx)
{
  print_n_dim<dim>(particles, nx);
}
