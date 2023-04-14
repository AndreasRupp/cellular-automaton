/**
 * @file building_units.hxx
 * @brief Implementation of different building units
 * CustomBU
 * HyperPlane
 * HyperSphere
 * Particle (Custom stencil)
 * TODO Maybe only Particle BU neccessary with different stencils (createStencil_X())
 *
 */
#pragma once

#include <CAM/utils.hxx>
#include <algorithm>
#include <array>
#include <chrono>
#include <iostream>
#include <memory>
#include <random>
#include <vector>
namespace CAM
{
/**
 * @brief Calculate maximal feret Diameter
 * Attentation: coords must be in R^n -> particle is not allowed to cross periodic domain boundary
 * Shift particle in the middle of domain
 * @tparam nx
 * @param _fields field indices of particle
 * @return max feret diamater
 */
template <auto nx>
double feretDiameter_max_fields(std::vector<unsigned int> _fields)
{
  unsigned int coord_a, coord_b;
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
            unsigned int bit_a = (c & (1 << i)) >> i;
            unsigned int bit_b = (d & (1 << i)) >> i;

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
/**
 * @brief Calculate maximal feret Diameter
 * Attentation: coords must be in R^n -> particle is not allowed to cross periodic domain boundary
 * Shift particle in the middle of domain
 * @tparam nx
 * @param _fields indices of points of cells
 * @return max feret diamater
 */
template <auto nx>
double feretDiameter_max(std::vector<unsigned int> _points)
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
/**
 * @brief Get the border cells and points of particle
 *TODO Implement convexHull
 * @tparam nx
 * @param _set connected cells
 * @return {borderCells, borderPoints}
 *
 */
template <auto nx>
std::array<std::vector<unsigned int>, 2> getBorder(std::vector<unsigned int> _set)
{
  std::vector<unsigned int> amountNeighbors(_set.size());
  std::fill(amountNeighbors.begin(), amountNeighbors.end(), 0);
  std::vector<unsigned int> borderCells;
  std::vector<unsigned int> borderPoints;
  for (unsigned int a = 0; a < _set.size(); a++)
  {
    unsigned int element = _set[a];
    std::array<unsigned int, n_fieldpoints<nx>()> points = CAM::getPoints<nx>(element);
    for (unsigned int i = 0; i < 2 * nx.size(); ++i)
    {
      unsigned int neigh = aim<nx>(element, direct_neigh<nx>(i));
      if (std::find(_set.begin(), _set.end(), neigh) != _set.end())
        amountNeighbors[a]++;
      else
      {
        unsigned int dim = i / 2;
        unsigned int lr = i % 2;
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
/**
 * @brief Base struct template of building units (BU)
 * \param number index of cells in domain
 * \param jump_parameter How far BU is allowed to jump.
 * \param referencePoints point in domain which is moved by CA and is calculation basis for all
 * cells in BU
 * @tparam nx
 */
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
/**
 * @brief Custom BU with defined field indices in domain
 *
 * @tparam nx
 */
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
/**
 * @brief Hyper sphere (2D: Circle, 3D: Sphere)
 *
 * @tparam nx
 * @param _position field index of center point; no position is given (-1) -> random position
 * @param _radius all cells are completely within the radius
 * @param _jump_parameter How far sphere is allowed to jump.
 */
template <auto nx>
struct HyperSphereBU : public BuildingUnit<nx>
{
  double radius;
  std::vector<unsigned int> stencilBU;

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
/**
 * @brief Limited hyper plane (2D: Rectangle, 3D: cuboid)
 * @param _position field index of center point; no position is given (-1) -> random position
 * @param _extent size of plane in each dimension
 * @param _jump_parameter How far hyper plane is allowed to jump.
 *
 * @tparam nx
 */
template <auto nx>
struct HyperPlaneBU : public BuildingUnit<nx>
{
  std::vector<unsigned int> extent;
  std::vector<unsigned int> stencilBU;
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
          unsigned int newCell = CAM::aim<nx>(stencilBU[i], d * direct_neigh<nx>(2 * dim + 1));
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
    // TODO
    return _index == 42;
  }
};
/**
 * @brief Particle with custom stencil
 * custom stencil defined by min and max feretDiameter and amount of extra cells attached to sphere
 * with radius min_feretDiameter/2
 *
 * @tparam nx
 */
template <auto nx>
struct ParticleBU : public BuildingUnit<nx>
{
  std::vector<unsigned int> stencilBU;
  ParticleBU(unsigned int _number,
             unsigned int _jump_parameter,
             unsigned int _seedPoint,
             std::vector<unsigned int> _stencil)
  {
    this->number = _number;
    this->referencePoints = {_seedPoint};
    this->jump_parameter = _jump_parameter;
    stencilBU = _stencil;
  }
  ~ParticleBU() override {}
  static const std::vector<unsigned int> getStencil(const int _random_seed,
                                                    const double _feretDiameter_max,
                                                    const double _feretDiameter_min,
                                                    const unsigned int _extraCells)
  {
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
    std::vector<unsigned int> stencil = CAM::getPNormedStencil<nx>(_feretDiameter_min / 2.0, 2);
    double feretDiameter;

    // std::array<std::vector<unsigned int>, 2> border_ = CAM::getBorder<nx>(stencil);
    // std::vector<unsigned int> border = border_[0];
    // std::vector<unsigned int> fields(border.size());
    // for (unsigned int i = 0; i < fields.size(); i++)
    //   fields[i] = aim<nx>(centerpoint, border[i]);
    // feretDiameter = CAM::feretDiameter_max_fields<nx>(fields);
    // std::cout << "max_feret1 " << feretDiameter << std::endl;

    // border = border_[1];
    // std::vector<unsigned int> fields1(border.size());
    // for (unsigned int i = 0; i < fields1.size(); i++)
    //   fields1[i] = aim<nx>(centerpoint, border[i]);
    // feretDiameter = CAM::feretDiameter_max<nx>(fields1);
    // std::cout << "max_feret2 " << feretDiameter << std::endl;

    unsigned int axis_min_feret = std::rand() % nx.size();
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

    for (unsigned int c = 0; c < _extraCells; c++)
    {
      std::array<std::vector<unsigned int>, 2> border_ = CAM::getBorder<nx>(stencil);
      std::vector<unsigned int> borderPoints = border_[1];
      std::vector<unsigned int> borderCells = border_[0];
      std::vector<unsigned int> points(borderPoints.size());
      for (unsigned int i = 0; i < borderPoints.size(); i++)
        points[i] = aim<nx>(centerpoint, borderPoints[i]);

      feretDiameter = CAM::feretDiameter_max<nx>(points);
      ;

      std::vector<unsigned int> newStencilCells;
      for (unsigned int a = 0; a < borderCells.size(); a++)
      {
        for (unsigned int i = 0; i < 2 * nx.size(); ++i)
        {
          unsigned int neigh = aim<nx>(borderCells[a], direct_neigh<nx>(i));
          unsigned int coord_min =
            (neigh / (int)direct_neigh<nx>(2 * axis_min_feret + 1) + n_fields<nx>()) %
            nx[axis_min_feret];
          if (std::find(stencil.begin(), stencil.end(), neigh) == stencil.end() &&
              find(min_feret_coord.begin(), min_feret_coord.end(), coord_min) !=
                min_feret_coord.end())
          {
            newStencilCells.push_back(neigh);
          }
        }
      }
      // here possible other metrics as smoothness
      std::shuffle(newStencilCells.begin(), newStencilCells.end(),
                   std::default_random_engine(std::rand()));
      for (unsigned int j = 0; j < newStencilCells.size(); j++)
      {
        std::array<unsigned int, n_fieldpoints<nx>()> newPoints =
          CAM::getPoints<nx>(newStencilCells[j]);
        for (unsigned int p = 0; p < newPoints.size(); p++)
          points.push_back(aim<nx>(centerpoint, (int)newPoints[p]));
        feretDiameter = CAM::feretDiameter_max<nx>(points);

        if (feretDiameter < _feretDiameter_max)
        {
          stencil.push_back(newStencilCells[j]);
          break;
        }
        else
        {
          for (unsigned int p = 0; p < newPoints.size(); p++)
            points.pop_back();
        }
      }
    }
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
  bool isMember(unsigned int _index) override
  {
    // TODO
    return _index == 42;
  }
};
}  // namespace CAM
