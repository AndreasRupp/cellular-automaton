/*!*********************************************************************************************
 * \file utils.hxx
 *
 * \brief Holds static utility functions which can be used in CAM
 *
 **************************************************************************************************/
#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <vector>
namespace CAM
{
typedef unsigned int fieldNumbers_t;
/*!*************************************************************************************************
 * \brief   Maximum unsigned integer.
 **************************************************************************************************/
static constexpr CAM::fieldNumbers_t uint_max = std::numeric_limits<fieldNumbers_t>::max();
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
 * \brief How many corner points of a cell
 *
 * \return number of cornerpoints
 **************************************************************************************************/
template <auto nx>
static constexpr unsigned int n_fieldpoints()
{
  return std::pow(2, nx.size());
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
 * \brief Get the all 2^dim points of a cell
 *
 * \tparam nx
 * \param _field index of cell
 * \return indices of points
 **************************************************************************************************/
template <auto nx>
static constexpr std::array<unsigned int, n_fieldpoints<nx>()> get_points(const unsigned int _field)
{
  std::array<unsigned int, n_fieldpoints<nx>()> points;
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
 * \tparam p type of norm (p = 2 euclidean norm)
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
    norm += std::pow(dist, p);
  }
  return std::pow(norm, 1.0 / (double)p);
}
/*!*********************************************************************************************
 * \brief Stencil with cells inside certain radius
 *
 * \tparam nx
 * \tparam p
 * \param _radius
 * \return std::vector<unsigned int>
 **************************************************************************************************/
template <auto nx, unsigned int p>
static constexpr std::vector<unsigned int> get_p_normed_particle(const double _radius)
{
  const unsigned int dim = nx.size();
  double radius = std::max(1., _radius);
  std::vector<unsigned int> stencil(1, 0);
  std::array<unsigned int, n_fieldpoints<nx>()> points;
  unsigned int newMove, index = 0, old_size = stencil.size();
  bool isInside, doNextLayer = true;
  while (doNextLayer)
  {
    doNextLayer = false;
    for (; index < old_size; ++index)
    {
      for (unsigned int i = 0; i < 2 * dim; ++i)
      {
        newMove = aim<nx>(stencil[index], direct_neigh<nx>(i));  // newCell = aim<nx>(0, newMove);
        points = CAM::get_points<nx>(newMove);
        isInside = true;
        for (unsigned int i = 0; i < points.size(); i++)
          isInside = isInside && (p_norm_distance<nx, p>(0, points[i]) < radius);

        if (std::find(stencil.begin(), stencil.end(), newMove) == stencil.end() && isInside)
        {
          stencil.push_back(newMove);
          doNextLayer = true;
        }
      }
    }
    old_size = stencil.size();
  }
  return stencil;
}
/*!*********************************************************************************************
 * \brief   Checks possible moves for a certain jump parameter.
 * \param   jump_parameter Manhattan Distance/L1 Distance
 * Order of possible moves: down, right, up, left.
 *
 * \retval stencil
 **********************************************************************************************/
template <auto nx>
static std::vector<unsigned int> get_stencil(double _jump_parameter)
{
  static const unsigned int dim = nx.size();
  std::vector<unsigned int> stencil(1, 0);
  unsigned int newNeigh, layers = std::max(1., _jump_parameter), index = 0,
                         old_size = stencil.size();
  for (unsigned int lay = 0; lay < layers; ++lay)
  {
    for (; index < old_size; ++index)
    {
      for (unsigned int i = 0; i < 2 * dim; ++i)
      {
        newNeigh = aim<nx>(stencil[index], direct_neigh<nx>(i));
        if (std::find(stencil.begin(), stencil.end(), newNeigh) == stencil.end())
        {
          stencil.push_back(newNeigh);
        }
      }
    }
    old_size = stencil.size();
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
static constexpr double feretDiameter_max_byFields(const std::vector<unsigned int>& _fields)
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
 * \param _fields indices of points of cells
 * \return max feret diamater
 **********************************************************************************************/
template <auto nx>
static constexpr double feretDiameter_max(const std::vector<unsigned int>& _points)
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
 * \brief Get the border cells and points of particle
 *TODO Implement convexHull
 * \tparam nx
 * \param _set connected cells
 * \return {borderCells, borderPoints}
 *
 **********************************************************************************************/
template <auto nx>
static constexpr std::array<std::vector<unsigned int>, 2> get_border(
  const std::vector<unsigned int>& _set)
{
  std::vector<unsigned int> borderCells, borderPoints, amountNeighbors(_set.size());
  std::fill(amountNeighbors.begin(), amountNeighbors.end(), 0);
  std::array<unsigned int, n_fieldpoints<nx>()> points;
  unsigned int element, neigh, dim, lr;
  for (unsigned int a = 0; a < _set.size(); a++)
  {
    element = _set[a];
    points = CAM::get_points<nx>(element);
    for (unsigned int i = 0; i < 2 * nx.size(); ++i)
    {
      neigh = aim<nx>(element, direct_neigh<nx>(i));
      if (std::find(_set.begin(), _set.end(), neigh) != _set.end())
        amountNeighbors[a]++;
      else
      {
        dim = i / 2;
        lr = i % 2;
        // TODO Shorter: iterate over pow(2,nx.size()-1) and at dim bit sandwich lr
        for (unsigned int j = 0; j < points.size(); j++)
        {
          if (((j & (1 << dim)) >> dim) == lr)
            borderPoints.push_back(points[j]);
        }
      }
    }
    if (amountNeighbors[a] < 2 * nx.size())
      borderCells.push_back(element);
  }
  sort(borderPoints.begin(), borderPoints.end());
  borderPoints.erase(unique(borderPoints.begin(), borderPoints.end()), borderPoints.end());
  return {borderCells, borderPoints};
}
}  // namespace CAM