#pragma once

#include <cmath>
#include <iomanip>
#include <iostream>

namespace CAM
{

/*!*************************************************************************************************
 * \brief   Calculates the size of the domain.
 *
 * \retval n_field        size of the domain
 **************************************************************************************************/
template <auto nx>
static constexpr unsigned int n_fields()
{
  unsigned int n_field = 1;
  for (unsigned int i = 0; i < nx.size(); ++i)
    n_field *= nx[i];
  return n_field;
}

template <auto nx>
static constexpr int direct_neigh(const unsigned int index)
{
  static_assert(nx.size() != 0, "Dimension of zero does not make sense.");
  int direct_neigh = (index % 2 == 0) ? -1 : 1;
  for (unsigned int i = 0; i < index / 2; ++i)
    direct_neigh *= nx[i];
  return direct_neigh;
}

/*!*************************************************************************************************
 * \brief   Find field if one moves from position to move.
 *
 * \param   position  Current position of field that may move.
 * \param   move      Index shift induced by possible move.
 * \retval  index     Index of move target.
 **************************************************************************************************/
template <auto nx>
static constexpr unsigned int aim(const unsigned int position, const int move)
{
  unsigned int coord, new_pos = 0;
  for (unsigned int i = 0; i < nx.size(); ++i)
  {
    coord = (position / direct_neigh<nx>(2 * i + 1) + move / (int)direct_neigh<nx>(2 * i + 1) +
             n_fields<nx>()) %
            nx[i];
    new_pos += coord * direct_neigh<nx>(2 * i + 1);
  }
  return new_pos;
}

template <typename fields_array_t>
static constexpr unsigned int bulk_distance(const fields_array_t& domain_a,
                                            const fields_array_t& domain_b)
{
  unsigned int distance = 0;
  for (unsigned int field = 0; field < domain_a.size(); ++field)
    distance += (domain_a[field] == 0) != (domain_b[field] == 0);
  return distance;
}

template <auto nx, typename fields_array_t>
static constexpr unsigned int skeleton_distance(const fields_array_t& domain_a,
                                                const fields_array_t& domain_b)
{
  unsigned int neigh_field, distance = 0;
  for (unsigned int field = 0; field < domain_a.size(); ++field)
    for (unsigned int j = 0; j < 2 * nx.size(); j += 2)
    {
      neigh_field = aim<nx>(field, direct_neigh<nx>(j));
      distance += std::abs((domain_a[field] == 0) - (domain_a[neigh_field] == 0) +
                           (domain_b[neigh_field] == 0) - (domain_b[field] == 0));
    }
  return distance;
}

// -------------------------------------------------------------------------------------------------
// PRINTING SECTION STARTS HERE
// -------------------------------------------------------------------------------------------------

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
      std::cout << std::setw(len_numbers) << std::setfill('0') << particles[index] << "  ";
    else
      std::cout << "\033[0;31m" << std::setw(len_numbers) << std::setfill('0') << particles[index]
                << "\033[0m  ";
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
          std::cout << "st";
        else if (i % 10 == 2 && i != 12)
          std::cout << "nd";
        else if (i % 10 == 3 && i != 13)
          std::cout << "rd";
        else
          std::cout << "th";
        std::cout << " coord: " << coord_i_helper % nx[i - 1] << "  ";
      }
      std::cout << std::endl;
    }
    for (unsigned int i = 0; i < nx[n_dim - 1]; ++i)
      print_n_dim<n_dim - 1>(particles, nx, (init_index + i) * nx[n_dim - 2]);
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
