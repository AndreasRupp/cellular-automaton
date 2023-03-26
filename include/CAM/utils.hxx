#pragma once
#include <algorithm>
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
    coord = ((position) / direct_neigh<nx>(2 * i + 1) +
             (move) / (int)direct_neigh<nx>(
                        2 * i + 1) +  // + negative werden in die falsche richtug aufgerundet
             n_fields<nx>()) %
            nx[i];
    new_pos += coord * direct_neigh<nx>(2 * i + 1);
  }
  return new_pos;
}
template <auto nx>
static constexpr unsigned int addMoves(int move1, int move2)
{
  unsigned int coord, new_move = 0;
  for (unsigned int i = 0; i < nx.size(); ++i)
  {
    coord = ((move1) / direct_neigh<nx>(2 * i + 1) + (move2) / direct_neigh<nx>(2 * i + 1) +
             n_fields<nx>()) %
            nx[i];
    new_move += coord * direct_neigh<nx>(2 * i + 1);
  }
  return new_move;
}
template <auto nx>
double pNormDistance(const unsigned int _position, const int _move, unsigned int _p)
{
  unsigned int coord, coord_moved;
  unsigned int norm = 0;
  for (unsigned int i = 0; i < nx.size(); ++i)
  {
    coord = (_position / (int)direct_neigh<nx>(2 * i + 1) + n_fields<nx>()) % nx[i];
    coord_moved = (_position / (int)direct_neigh<nx>(2 * i + 1) +
                   _move / (int)direct_neigh<nx>(2 * i + 1) + n_fields<nx>()) %
                  nx[i];
    unsigned int d = std::abs((int)coord - (int)coord_moved);
    d = std::min(d, nx[i] - d);
    norm += std::pow(d, _p);
  }
  return std::pow(norm, 1.0 / (double)_p);
}
template <auto nx>
static std::vector<int> getPNormedStencil(double _radius, unsigned int _p)
{
  static const unsigned int dim = nx.size();
  double radius = std::max(1., _radius);
  std::vector<int> stencil(1, 0);
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
        int newNeigh = addMoves<nx>(stencil[index], direct_neigh<nx>(i));
        if (std::find(stencil.begin(), stencil.end(), newNeigh) == stencil.end() &&
            pNormDistance<nx>(0, newNeigh, _p) < radius)
        {
          stencil.push_back(newNeigh);
          doNextLayer = true;
        }
      }
    }
    old_size = stencil.size();
  }
  return stencil;
}
template <auto nx>
static std::vector<int> getStencil(double _jump_parameter)
{
  static const unsigned int dim = nx.size();
  unsigned int layers = std::max(1., _jump_parameter);
  std::vector<int> stencil(1, 0);
  unsigned int index = 0;
  unsigned int old_size = stencil.size();

  for (unsigned int lay = 0; lay < layers; ++lay)
  {
    for (; index < old_size; ++index)
    {
      for (unsigned int i = 0; i < 2 * dim; ++i)
      {
        int newNeigh = addMoves<nx>(stencil[index], direct_neigh<nx>(i));
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

}  // namespace CAM