/**
 * @file utils.hxx
 *
 * @brief Holds static utility functions which can be used in CAM
 *
 */
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
/**
 * @brief Finds the movement to a direct neighbor within von Neumann neighborhood
 * each dimension two neighbors (left and right)
 * @param index 0 < index < 2 * nx.size(). index of neighbor
 * @return move
 */
template <auto nx>
static constexpr int direct_neigh(const unsigned int index)
{
  static_assert(nx.size() != 0, "Dimension of zero does not make sense.");
  int direct_neigh = (index % 2 == 0) ? -1 : 1;
  for (unsigned int i = 0; i < index / 2; ++i)
    direct_neigh *= nx[i];
  return direct_neigh;
}
/**
 * @brief How many corner points of a cell
 *
 * @return number of cornerpoints
 */
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
/**
 * @brief Get the all 2^dim points of a cell
 *
 * @tparam nx
 * @param _field index of cell
 * @return indices of points
 */
template <auto nx>
static std::array<unsigned int, n_fieldpoints<nx>()> getPoints(unsigned int _field)
{
  std::array<unsigned int, n_fieldpoints<nx>()> points;
  for (unsigned int a = 0; a < pow(2, nx.size()); a++)
  {
    points[a] = _field;
    for (unsigned int i = 0; i < nx.size(); ++i)
    {
      unsigned int leftright = (a & (1 << i)) >> i;
      points[a] = aim<nx>(points[a], leftright * direct_neigh<nx>(2 * i + 1));
    }
  }
  // TODO shorter without using aim
  return points;
}
/**
 * @brief Distance induced by p-Norm of two points/cells in a periodic domain
 *
 * @tparam nx
 * @param _position1
 * @param _position2
 * @param _p type of norm (p = 2 euclidean norm)
 * @return distance
 */
template <auto nx>
double pNormDistance(const unsigned int _position1, const unsigned int _position2, unsigned int _p)
{
  unsigned int coord1, coord2;
  unsigned int norm = 0;
  for (unsigned int i = 0; i < nx.size(); ++i)
  {
    coord1 = (_position1 / direct_neigh<nx>(2 * i + 1) + n_fields<nx>()) % nx[i];
    coord2 = (_position2 / direct_neigh<nx>(2 * i + 1) + n_fields<nx>()) % nx[i];
    std::array<int, 3> ds;
    ds[0] = (int)coord1 - (int)coord2;
    ds[1] = ((int)coord1 - (int)coord2 + nx[i]) % nx[i];
    ds[2] = ds[1] - nx[i];
    ds[0] = std::abs(ds[0]);
    ds[1] = std::abs(ds[1]);
    ds[2] = std::abs(ds[2]);
    auto d = *std::min_element(ds.begin(), ds.end());

    norm += std::pow(d, _p);
  }
  return std::pow(norm, 1.0 / (double)_p);
}
/**
 * @brief Stencil with cells inside certain radius
 *
 * @tparam nx
 * @param _radius
 * @param _p
 * @return std::vector<unsigned int>
 */
template <auto nx>
static std::vector<unsigned int> getPNormedStencil(double _radius, unsigned int _p)
{
  // Distance vielleicht mit Punkten machen
  static const unsigned int dim = nx.size();
  double radius = std::max(1., _radius);
  std::vector<unsigned int> stencil(1, 0);
  unsigned int index = 0;
  unsigned int old_size = stencil.size();
  bool doNextLayer = true;
  while (doNextLayer)
  {
    doNextLayer = false;
    for (; index < old_size; ++index)
    {
      for (unsigned int i = 0; i < 2 * dim; ++i)
      {
        unsigned int newMove = aim<nx>(stencil[index], direct_neigh<nx>(i));
        unsigned int newCell = newMove;  // aim<nx>(0, newMove);
        std::array<unsigned int, n_fieldpoints<nx>()> points = CAM::getPoints<nx>(newCell);
        bool isInside = true;
        for (unsigned int i = 0; i < points.size(); i++)
          isInside = isInside && (pNormDistance<nx>(0, points[i], _p) < radius);

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
 * \param   jump_parameter Manhattan Distance
 * Order of possible moves: down, right, up, left.
 *
 * \retval stencil
 **********************************************************************************************/
template <auto nx>
static std::vector<unsigned int> getStencil(double _jump_parameter)
{
  static const unsigned int dim = nx.size();
  unsigned int layers = std::max(1., _jump_parameter);
  std::vector<unsigned int> stencil(1, 0);
  unsigned int index = 0;
  unsigned int old_size = stencil.size();

  for (unsigned int lay = 0; lay < layers; ++lay)
  {
    for (; index < old_size; ++index)
    {
      for (unsigned int i = 0; i < 2 * dim; ++i)
      {
        unsigned int newNeigh = aim<nx>(stencil[index], direct_neigh<nx>(i));
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
template <typename T>
std::vector<int> findItems(std::vector<T> const& v, T target)
{
  std::vector<int> indices;

  for (auto it = v.begin(); it != v.end(); it++)
  {
    if (*it == target)
    {
      indices.push_back(std::distance(v.begin(), it));
    }
  }

  return indices;
}
}  // namespace CAM