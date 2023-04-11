#pragma once

#include <CAM/utils.hxx>
#include <algorithm>
#include <array>
#include <chrono>
#include <iostream>
#include <memory>
#include <vector>
namespace CAM
{
// static std::vector<unsigned int> convexHull(std::vector<unsigned int> _points)
// {
//   std::vector<unsigned int> convexHull;

//   for(unsigned int a = 0; a < _points.size(); a++)
//     {
//       unsigned int distance = 0;
//       for (unsigned int i = 0; i < nx.size(); ++i)
//       {
//       }
//     }

// }
// template<auto nx>
// std::vector<std::array<unsigned int,nx.size()>> getConvexHull(std::vector<std::array<unsigned
// int, nx.size()>> _points)
// {
//   std::vector<array<unsigned int, nx.size()>> convexHull;
//   convexHull.push_back(_points[0]);
//   for(int i = 1; i < _points.size(); i++)
//   {
//     for(int j = 0; j < nx.size(); j++)
//     {

//     }
//   }

// }

template <auto nx>
double feretDiameter_max(std::vector<unsigned int> _convexHull)
{
  // man geht durch array, verlgleicht immer akutellstes elemt mit elementen mögicher convex hull
  //  wenn das element in einer richtung extremer ist, wird mit dem nicht extremen ausgetauscht.
  // es müssen n-1 gleich sein und 1 unterschiedlich
  unsigned int coord_a, coord_b;
  double distance_max = 0;
  for (unsigned int a = 0; a < _convexHull.size(); a++)
  {
    for (unsigned int b = 0; b < _convexHull.size(); b++)
    {
      for (unsigned int c = 0; c < pow(2, nx.size()); c++)
      {
        for (unsigned int d = 0; d < pow(2, nx.size()); d++)
        {
          double distance = 0;
          for (unsigned int i = 0; i < nx.size(); ++i)
          {
            int bit_a = (c & (1 << i)) >> i;
            int bit_b = (d & (1 << i)) >> i;

            coord_a = (_convexHull[a] / (int)direct_neigh<nx>(2 * i + 1) + n_fields<nx>()) % nx[i];
            coord_a = coord_a + bit_a;
            coord_b = (_convexHull[b] / (int)direct_neigh<nx>(2 * i + 1) + n_fields<nx>()) % nx[i];
            coord_b = coord_b + bit_b;
            // std::cout<<coord_a << " "<<coord_b<<std::endl;
            distance += (coord_a - coord_b) * (coord_a - coord_b);
          }
          distance = sqrt(distance);
          if (distance > distance_max)
            distance_max = distance;
          // std::cout<<"-----"<<std::endl;
        }
      }
    }
  }
  return distance_max;
}
template <auto nx>
// Points, Cells
std::array<std::vector<unsigned int>, 2> getBorder(std::vector<int> _set)
{
  std::vector<unsigned int> amountNeighbors(_set.size());
  std::fill(amountNeighbors.begin(), amountNeighbors.end(), 0);
  std::vector<unsigned int> borderCells;
  std::vector<unsigned int> borderPoints;
  for (unsigned int a = 0; a < _set.size(); a++)
  {
    int element = _set[a];
    for (unsigned int i = 0; i < 2 * nx.size(); ++i)
    {
      int neigh = aim<nx>(element, direct_neigh<nx>(i));
      if (std::find(_set.begin(), _set.end(), neigh) != _set.end())
        amountNeighbors[a]++;
      else  // alle nachbarn in positiver richtung außer die wo man schon hingegangen ist
        if (i % 2 == 1)
        {
          borderPoints.push_back(neigh);
          for (unsigned int j = 0; j < nx.size(); j++)
          {
            if ((2 * j + 1) != i)
            {
              int neigh2 = aim<nx>(neigh, direct_neigh<nx>(2 * j + 1));
              borderPoints.push_back(neigh2);
            }
          }
        }
    }

    if (amountNeighbors[a] < 2 * nx.size())
    {
      borderCells.push_back(_set[a]);
      borderPoints.push_back(_set[a]);
    }
  }
  sort(borderPoints.begin(), borderPoints.end());
  borderPoints.erase(unique(borderPoints.begin(), borderPoints.end()), borderPoints.end());
  // for(int i = 0; i<borderPoints.size(); i++)
  //   std::cout<<borderPoints[i] <<std::endl;
  return {borderCells, borderPoints};
}
// alle nachbarn mit 1er aus border raus, alle nx.size* 2 raus
template <auto nx>
struct BuildingUnit
{
 public:
  virtual ~BuildingUnit(){};
  CAM::fieldNumbers_t number;
  unsigned int jump_parameter;
  std::vector<unsigned int> referencePoints;

  virtual std::vector<unsigned int> getFieldIndices() = 0;

  virtual bool isMember(unsigned int _index) = 0;
};
template <auto nx>
struct CustomBU : public BuildingUnit<nx>
{
  CustomBU(unsigned int _number,
           unsigned int _jump_parameter,
           std::vector<unsigned int> _fieldIndices)
  {
    this->number = _number;
    this->referencePoints = _fieldIndices;
    this->jump_parameter = _jump_parameter;
  }
  ~CustomBU() override {}

  std::vector<unsigned int> getFieldIndices() override { return this->referencePoints; }
  bool isMember(unsigned int _index) override
  {
    return (std::find(this->referencePoints.begin(), this->referencePoints.end(), _index) !=
            this->referencePoints.end());
  }
};
template <auto nx>
struct HyperSphereBU : public BuildingUnit<nx>
{
  double radius;
  std::vector<int> stencilBU;

  HyperSphereBU(unsigned int _number,
                unsigned int _jump_parameter,
                unsigned int _centerPoint,
                double _radius)
  {
    this->number = _number;
    this->referencePoints = {_centerPoint};
    radius = _radius;
    this->jump_parameter = _jump_parameter;
    stencilBU = CAM::getPNormedStencil<nx>(radius, 2);
  }
  ~HyperSphereBU() override {}

  std::vector<unsigned int> getFieldIndices() override
  {
    std::vector<unsigned int> fields;
    for (unsigned int i = 0; i < stencilBU.size(); i++)
    {
      fields.push_back(CAM::aim<nx>(this->referencePoints[0], stencilBU[i]));
    }
    return fields;
  }
  bool isMember(unsigned int _index) override
  {
    return (CAM::pNormDistance<nx>(this->referencePoints[0], (int)_index, 2) < radius);
  }
};
template <auto nx>
struct HyperPlaneBU : public BuildingUnit<nx>
{
  std::vector<unsigned int> extent;
  std::vector<int> stencilBU;
  HyperPlaneBU(unsigned int _number,
               unsigned int _jump_parameter,
               unsigned int _referencePoint,
               std::vector<unsigned int> _extent)
  {
    this->number = _number;
    this->referencePoints = {_referencePoint};
    this->jump_parameter = _jump_parameter;
    extent = _extent;

    stencilBU.push_back(0);
    // static_assert(_extent.size() == nx.size(), "Size");

    for (unsigned int dim = 0; dim < nx.size(); dim++)
    {
      unsigned int size = stencilBU.size();
      for (unsigned int d = 1; d < extent[dim]; d++)
      {
        for (unsigned int i = 0; i < size; i++)
        {
          unsigned int newCell = CAM::addMoves<nx>(stencilBU[i], d * direct_neigh<nx>(2 * dim + 1));
          stencilBU.push_back(newCell);
        }
      }
    }
  }
  ~HyperPlaneBU() override {}
  std::vector<unsigned int> getFieldIndices() override
  {
    std::vector<unsigned int> fields;
    for (unsigned int i = 0; i < stencilBU.size(); i++)
    {
      fields.push_back(CAM::aim<nx>(this->referencePoints[0], stencilBU[i]));
    }
    return fields;
  }
  bool isMember(unsigned int _index) override
  {
    // int coord;
    // int coord_max;
    // for (unsigned int i = 0; i < nx.size(); ++i)
    // {
    //   CAM::pNormDistance<nx>(this->referencePoints[0],_index) (int)direct_neigh<nx>(2 * i + 1)
    // }
    return true;
  }
};
// Vielleicht brauche man diese struct auch nicht, sonder  alles nur CustomBU
//+ create Stencil
template <auto nx>
struct ParticleBU : public BuildingUnit<nx>
{
  std::vector<int> stencilBU;
  ParticleBU(unsigned int _number,
             unsigned int _jump_parameter,
             unsigned int _seedPoint,
             std::vector<int> _stencil)
  {
    this->number = _number;
    this->referencePoints = {_seedPoint};  // vector
    this->jump_parameter = _jump_parameter;
    stencilBU = _stencil;
  }
  ~ParticleBU() override {}
  static const std::vector<int> getStencil(const int _random_seed,
                                           const double _feretDiameter_max,
                                           const double _feretDiameter_min)
  {
    // mimum feret durch kreis mit 2 * radius = feret
    // dann voxel dran machen, sodass max feret diameter erreicht wird.
    //  dann voxel dran machen (bis area erfüllt, oder verhältnis umfang- area),  mit bedingung
    //  max_feret nicht überboten wird.

    // schwerpunkt d
    // punkt finden dass mit max feret diamter passt
    unsigned int rand_seed;
    if (_random_seed == 0)
    {
      rand_seed = std::chrono::system_clock::now().time_since_epoch().count();
    }
    else
      rand_seed = _random_seed;
    std::srand(rand_seed);

    unsigned int centerpoint = 0;
    for (unsigned int i = 0; i < nx.size(); ++i)
    {
      centerpoint += nx[i] / 2 * direct_neigh<nx>(2 * i + 1);
    }
    std::vector<int> stencil = CAM::getPNormedStencil<nx>(_feretDiameter_min / 2.0, 2);

    static std::array<std::vector<unsigned int>, 2> border_ = CAM::getBorder<nx>(stencil);
    std::vector<unsigned int> border = border_[0];
    std::vector<unsigned int> fields(border.size());
    for (unsigned int i = 0; i < fields.size(); i++)
      fields[i] = aim<nx>(centerpoint, border[i]);
    double feretDiameter;
    feretDiameter = CAM::feretDiameter_max<nx>(fields);
    std::cout << "max_feret1 " << feretDiameter << std::endl;

    // border = border_[1];
    // std::vector<unsigned int> fields1(border.size());
    // for(unsigned int i = 0; i < fields1.size(); i++)
    //     fields1[i] = aim<nx>(centerpoint,border[i]);
    // feretDiameter = CAM::feretDiameter_max<nx>(fields1);
    // std::cout<<"max_feret2 "<<feretDiameter<<std::endl;

    int axis_min_feret = std::rand() % nx.size();
    std::vector<unsigned int> min_feret_coord;
    for (unsigned int i = 0; i < stencil.size(); i++)
    {
      min_feret_coord.push_back(
        (stencil[i] / (int)direct_neigh<nx>(2 * axis_min_feret + 1) + n_fields<nx>()) %
        nx[axis_min_feret]);
    }
    sort(min_feret_coord.begin(), min_feret_coord.end());
    min_feret_coord.erase(unique(min_feret_coord.begin(), min_feret_coord.end()),
                          min_feret_coord.end());
    // std::for_each(min_feret_coord.begin(),min_feret_coord.end(),[] (unsigned int& c)
    // {
    //   std::cout<<c<<""<<std::endl;
    // });

    // der abstand muss gleich bleiben
    // min max rausfinden: _feretDiameter_min/2 * (int)direct_neigh<nx>(2 * axis_min_feret + 1)
    // coord_a = (_convexHull[a] / (int)direct_neigh<nx>(2 * axis_min_feret + 1) + n_fields<nx>()) %
    // nx[axis_min_feret];

    int loops = 0;
    do
    {
      loops++;
      std::cout << "looop " << loops << std::endl;
      std::array<std::vector<unsigned int>, 2> border_ = CAM::getBorder<nx>(stencil);
      std::vector<unsigned int> border = border_[0];
      std::vector<unsigned int> fields(border.size());
      for (unsigned int i = 0; i < border.size(); i++)
        fields[i] = aim<nx>(centerpoint, border[i]);

      std::vector<int> newStencilPoints;
      for (unsigned int a = 0; a < border.size(); a++)
      {
        for (unsigned int i = 0; i < 2 * nx.size(); ++i)
        {
          unsigned int neigh = aim<nx>(border[a], direct_neigh<nx>(i));
          if (std::find(stencil.begin(), stencil.end(), neigh) == stencil.end())
          {
            fields.push_back(aim<nx>(centerpoint, neigh));
            feretDiameter = CAM::feretDiameter_max<nx>(fields);
            // std::cout<<"max_feret "<<feretDiameter<<std::endl;
            fields.pop_back();
            //
            unsigned int coord_min =
              (neigh / (int)direct_neigh<nx>(2 * axis_min_feret + 1) + n_fields<nx>()) %
              nx[axis_min_feret];
            if (feretDiameter < _feretDiameter_max &&
                find(min_feret_coord.begin(), min_feret_coord.end(), coord_min) !=
                  min_feret_coord.end())
              newStencilPoints.push_back(neigh);
          }
        }
      }
      int p = newStencilPoints[std::rand() % newStencilPoints.size()];
      stencil.push_back(p);
    } while (loops < 25);  // feretDiameter < _feretDiameter_max &&

    return stencil;
  }
  std::vector<unsigned int> getFieldIndices() override
  {
    std::vector<unsigned int> fields;
    for (unsigned int i = 0; i < stencilBU.size(); i++)
    {
      fields.push_back(CAM::aim<nx>(this->referencePoints[0], stencilBU[i]));
    }
    return fields;
  }
  bool isMember(unsigned int _index) override { return _index == 42; }
};
}  // namespace CAM
