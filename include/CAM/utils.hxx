/*!*********************************************************************************************
 * \file utils.hxx
 *
 * \brief Holds static utility functions which can be used in CAM
 *
 **************************************************************************************************/
#pragma once
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iostream>
#include <random>
#include <type_traits>
#include <vector>
namespace CAM
{
/*!*************************************************************************************************
 * \brief   Maximum unsigned integer.
 **************************************************************************************************/
static constexpr unsigned int uint_max = std::numeric_limits<unsigned int>::max();
/*!***********************************************************************************************
 * \brief   Smallest (negative) double.
 ************************************************************************************************/
static constexpr double double_min = std::numeric_limits<double>::lowest();
/*!***********************************************************************************************
 * \brief  Biggest double.
 ************************************************************************************************/
static constexpr double double_max = std::numeric_limits<double>::max();

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
static constexpr unsigned int get_random_field_index()
{
  unsigned int rand_seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::srand(rand_seed);
  return std::rand() % (n_fields<nx>());
}
/*!*********************************************************************************************
 * \brief Finds the movement to a direct neighbor within von Neumann neighborhood
 * each dimension two neighbors (left and right)
 * \param index 0 < index < 2 * nx.size(). index of neighbor
 * \return move
 **************************************************************************************************/
template <auto nx>
static constexpr int direct_neigh(const unsigned int index)
{
  static_assert(nx.size() != 0, "Dimension of zero does not make sense.");
  int direct_neigh = (index % 2 == 0) ? -1 : 1;
  for (unsigned int i = 0; i < index / 2; ++i)
    direct_neigh *= nx[i];
  return direct_neigh;
}
/*!*********************************************************************************************
 * \brief Constexpr power function
 **************************************************************************************************/
template <typename T>
constexpr T ipow(T num, unsigned int pow)
{
  return pow == 0 ? 1 : num * ipow(num, pow - 1);
}
/*!*********************************************************************************************
 * \brief How many corner points of a cell
 *
 * \return number of cornerpoints
 **************************************************************************************************/
template <auto nx>
static constexpr unsigned int n_field_corner_points()
{
  return ipow(2, nx.size());
}
/*!*************************************************************************************************
 * \brief   Find field if one moves from position to move.
 *
 * \param   position  Current position of field that may move.
 * \param   move      Index shift induced by possible move.
 * \retval  index     Index of move target.
 **************************************************************************************************/
template <auto nx>
static constexpr unsigned int aim(const int position, const int move)
{
  unsigned int coord, new_pos = 0;
  for (unsigned int i = 0; i < nx.size(); ++i)
  {
    coord = ((position) / direct_neigh<nx>(2 * i + 1) + (move) / direct_neigh<nx>(2 * i + 1) +
             n_fields<nx>()) %
            nx[i];
    new_pos += coord * direct_neigh<nx>(2 * i + 1);
  }
  return new_pos;
}
/*!*********************************************************************************************
 * \brief Get the all 2^dim corner points of a cell
 *
 * \tparam nx
 * \param _field index of cell
 * \return indices of points
 **************************************************************************************************/
template <auto nx>
static constexpr std::array<unsigned int, n_field_corner_points<nx>()> get_corner_points(
  const unsigned int _field)
{
  std::array<unsigned int, n_field_corner_points<nx>()> points;
  unsigned int leftright;
  for (unsigned int a = 0; a < pow(2, nx.size()); a++)
  {
    points[a] = _field;
    for (unsigned int i = 0; i < nx.size(); ++i)
    {
      leftright = (a & (1 << i)) >> i;
      points[a] = aim<nx>(points[a], leftright * direct_neigh<nx>(2 * i + 1));
    }
  }
  // TODO shorter without using aim
  return points;
}
/*!*********************************************************************************************
 * \brief Distance induced by p-Norm of two points/cells in a periodic domain
 *
 * \tparam nx
 * \tparam p type of norm (p == 2: euclidean norm, p == 0: maximum norm)
 * \param _position1
 * \param _position2
 * \return distance
 **************************************************************************************************/
template <auto nx, unsigned int p>
double constexpr p_norm_distance(const unsigned int _position1, const unsigned int _position2)
{
  unsigned int coord1, coord2, dist;
  double norm = 0;
  for (unsigned int i = 0; i < nx.size(); ++i)
  {
    coord1 = (_position1 / direct_neigh<nx>(2 * i + 1) + n_fields<nx>()) % nx[i];
    coord2 = (_position2 / direct_neigh<nx>(2 * i + 1) + n_fields<nx>()) % nx[i];
    dist = std::abs((int)coord1 - (int)coord2);
    if (dist > nx[i] / 2)
      dist = nx[i] - dist;
    if ((p == 0) && (norm < dist))
      norm = dist;
    else if (p != 0)
      norm += std::pow(dist, p);
  }
  if (p == 0)
    return norm;
  else
    return std::pow(norm, 1.0 / (double)p);
}
/*!*********************************************************************************************
 * \brief shape with cells inside certain radius
 *
 * \tparam nx
 * \tparam p defines type of norm/distance
 * \param _radius maximum distance of a cell from the central point
 * \return std::vector<unsigned int>
 **************************************************************************************************/
template <auto nx, unsigned int p>
static constexpr std::vector<unsigned int> get_p_normed_shape(const double _radius)
{
  double radius = std::max(1., _radius);
  std::vector<unsigned int> shape(1, 0);
  std::array<unsigned int, n_field_corner_points<nx>()> points;
  unsigned int new_move, index = 0, old_size = shape.size();
  bool isInside, do_next_layer = true;
  while (do_next_layer)
  {
    do_next_layer = false;
    for (; index < old_size; ++index)
    {
      for (unsigned int i = 0; i < 2 * nx.size(); ++i)
      {
        new_move = aim<nx>(shape[index], direct_neigh<nx>(i));  // newCell = aim<nx>(0, new_move);
        points = CAM::get_corner_points<nx>(new_move);
        isInside = true;
        for (unsigned int i = 0; i < points.size(); i++)
          isInside = isInside && (p_norm_distance<nx, p>(0, points[i]) < radius);

        if (std::find(shape.begin(), shape.end(), new_move) == shape.end() && isInside)
        {
          shape.push_back(new_move);
          do_next_layer = true;
        }
      }
    }
    old_size = shape.size();
  }
  return shape;
}
/*!*********************************************************************************************
 * \brief   Checks possible moves for a certain jump parameter.
 * \param   jump_parameter Manhattan Distance/L1 Distance
 * Order of possible moves: down, right, up, left.
 *
 * \retval stencil
 **********************************************************************************************/
template <auto nx>
static constexpr std::vector<unsigned int> get_stencil(const double _jump_parameter)
{
  // if(std::is_constant_evaluated())
  // {
  //   std::cout<<"compile time"<<std::endl;
  // }
  // else
  // {
  //   std::cout<<"runtime"<<std::endl;
  // }
  std::vector<unsigned int> stencil(1, 0);
  unsigned int new_neigh, layers = std::max(1., _jump_parameter), index = 0,
                          old_size = stencil.size();
  for (unsigned int lay = 0; lay < layers; ++lay)
  {
    for (; index < old_size; ++index)
    {
      for (unsigned int i = 0; i < 2 * nx.size(); ++i)
      {
        new_neigh = aim<nx>(stencil[index], direct_neigh<nx>(i));
        if (std::find(stencil.begin(), stencil.end(), new_neigh) == stencil.end())
        {
          stencil.push_back(new_neigh);
        }
      }
    }
    old_size = stencil.size();
  }
  return stencil;
}
template <auto nx, unsigned int jump_parameter>
static constexpr unsigned int get_stencil_size()
{
  // if(std::is_constant_evaluated())
  // {
  //   std::cout<<"compile time"<<std::endl;
  // }
  // else
  // {
  //   std::cout<<"runtime"<<std::endl;
  // }
  std::array<unsigned int, ipow((2 * jump_parameter + 1), nx.size())> stencil;
  std::fill(stencil.begin(), stencil.end(), n_fields<nx>() + 1);
  stencil[0] = 0;
  unsigned int new_neigh, layers = jump_parameter, index = 0, nr_of_cells = 1;
  unsigned int old_size = 1;
  for (unsigned int lay = 0; lay < layers; ++lay)
  {
    for (; index < old_size; ++index)
    {
      for (unsigned int i = 0; i < 2 * nx.size(); ++i)
      {
        new_neigh = aim<nx>(stencil[index], direct_neigh<nx>(i));
        if (std::find(stencil.begin(), stencil.end(), new_neigh) == stencil.end())
        {
          stencil[nr_of_cells] = new_neigh;
          nr_of_cells++;
        }
      }
    }
    old_size = nr_of_cells;
  }
  return nr_of_cells;
}
template <auto nx, unsigned int jump_parameter>
static constexpr std::array<unsigned int, get_stencil_size<nx, jump_parameter>()> get_stencil_c()
{
  // if(std::is_constant_evaluated())
  // {
  //   std::cout<<"compile time"<<std::endl;
  // }
  // else
  // {
  //   std::cout<<"runtime"<<std::endl;
  // }
  std::array<unsigned int, get_stencil_size<nx, jump_parameter>()> stencil;
  std::fill(stencil.begin(), stencil.end(), n_fields<nx>() + 1);
  stencil[0] = 0;

  unsigned int new_neigh, layers = jump_parameter, index = 0, nr_of_cells = 1;
  unsigned int old_size = 1;
  for (unsigned int lay = 0; lay < layers; ++lay)
  {
    for (; index < old_size; ++index)
    {
      for (unsigned int i = 0; i < 2 * nx.size(); ++i)
      {
        new_neigh = aim<nx>(stencil[index], direct_neigh<nx>(i));
        if (std::find(stencil.begin(), stencil.end(), new_neigh) == stencil.end())
        {
          stencil[nr_of_cells] = new_neigh;
          nr_of_cells++;
        }
      }
    }
    old_size = nr_of_cells;
  }
  return stencil;
}

/*!*********************************************************************************************
 * \brief Calculate maximal feret Diameter
 * Attentation: coords must be in R^n -> particle is not allowed to cross periodic domain boundary
 * Shift particle in the middle of domain
 * \tparam nx
 * \param _fields field indices of particle
 * \return max feret diamater
 **********************************************************************************************/
template <auto nx>
static constexpr double feret_diameter_max_by_fields(const std::vector<unsigned int>& _fields)
{
  unsigned int coord_a, coord_b, bit_a, bit_b;
  double distance_max = 0;
  for (unsigned int a = 0; a < _fields.size(); a++)
  {
    for (unsigned int b = 0; b < _fields.size(); b++)
    {
      for (unsigned int c = 0; c < pow(2, nx.size()); c++)
      {
        for (unsigned int d = 0; d < pow(2, nx.size()); d++)
        {
          double distance = 0;
          for (unsigned int i = 0; i < nx.size(); ++i)
          {
            bit_a = (c & (1 << i)) >> i;
            bit_b = (d & (1 << i)) >> i;

            coord_a = (_fields[a] / (int)direct_neigh<nx>(2 * i + 1) + n_fields<nx>()) % nx[i];
            coord_a = coord_a + bit_a;
            coord_b = (_fields[b] / (int)direct_neigh<nx>(2 * i + 1) + n_fields<nx>()) % nx[i];
            coord_b = coord_b + bit_b;

            distance += (coord_a - coord_b) * (coord_a - coord_b);
          }
          distance = sqrt(distance);
          if (distance > distance_max)
            distance_max = distance;
        }
      }
    }
  }
  return distance_max;
}
/*!*********************************************************************************************
 * \brief Calculate maximal feret Diameter
 * Attentation: coords must be in R^n -> particle is not allowed to cross periodic domain boundary
 * Shift particle in the middle of domain
 * \tparam nx
 * \param _fields indices of corner points of cells
 * \return max feret diamater
 **********************************************************************************************/
template <auto nx>
static constexpr double feret_diameter_max(const std::vector<unsigned int>& _points)
{
  unsigned int coord_a, coord_b;
  double distance, distance_max = 0;

  for (unsigned int a = 0; a < _points.size(); a++)
  {
    for (unsigned int b = 0; b < _points.size(); b++)
    {
      distance = 0;
      for (unsigned int i = 0; i < nx.size(); ++i)
      {
        coord_a = (_points[a] / (int)direct_neigh<nx>(2 * i + 1) + n_fields<nx>()) % nx[i];
        coord_b = (_points[b] / (int)direct_neigh<nx>(2 * i + 1) + n_fields<nx>()) % nx[i];
        distance += (coord_a - coord_b) * (coord_a - coord_b);
      }
      distance = sqrt(distance);
      if (distance > distance_max)
        distance_max = distance;
    }
  }
  return distance_max;
}

/*!*********************************************************************************************
 * \brief Get the boundary cells and points of particle
 *TODO Implement convexHull
 * \tparam nx
 * \param _fields connected cells
 * \return {boundary_cells, boundary_corner_points}
 *
 **********************************************************************************************/
template <auto nx>
static constexpr std::array<std::vector<unsigned int>, 2> get_boundary(
  const std::vector<unsigned int>& _fields)
{
  std::vector<unsigned int> boundary_fields, boundary_points, amount_neighbors(_fields.size());
  std::fill(amount_neighbors.begin(), amount_neighbors.end(), 0);
  std::array<unsigned int, n_field_corner_points<nx>()> points;
  unsigned int element, neigh, dim, lr;
  for (unsigned int a = 0; a < _fields.size(); a++)
  {
    element = _fields[a];
    points = CAM::get_corner_points<nx>(element);
    for (unsigned int i = 0; i < 2 * nx.size(); ++i)
    {
      neigh = aim<nx>(element, direct_neigh<nx>(i));
      if (std::find(_fields.begin(), _fields.end(), neigh) != _fields.end())
        amount_neighbors[a]++;
      else
      {
        dim = i / 2;
        lr = i % 2;
        // TODO Shorter: iterate over pow(2,nx.size()-1) and at dim bit sandwich lr
        for (unsigned int j = 0; j < points.size(); j++)
        {
          // bit value represent dimension and 0 or 1 defines which face side of hyper cube
          if (((j & (1 << dim)) >> dim) == lr)
            boundary_points.push_back(points[j]);
        }
      }
    }
    if (amount_neighbors[a] < 2 * nx.size())
      boundary_fields.push_back(element);
  }
  sort(boundary_points.begin(), boundary_points.end());
  boundary_points.erase(unique(boundary_points.begin(), boundary_points.end()),
                        boundary_points.end());
  return {boundary_fields, boundary_points};
}
template <auto nx>
static constexpr std::vector<unsigned int> get_boundary_fields(
  const std::vector<unsigned int>& _fields)
{
  return get_boundary<nx>(_fields)[0];
}
template <auto nx>
static constexpr std::vector<unsigned int> get_boundary_corner_points(
  const std::vector<unsigned int>& _fields)
{
  return get_boundary<nx>(_fields)[1];
}
}  // namespace CAM
