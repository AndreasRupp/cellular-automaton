/*!*********************************************************************************************
 * \file building_units.hxx
 * \brief Implementation of different building units
 * CustomBU
 * HyperPlane
 * HyperSphere
 * Particle (Custom stencil)
 * TODO Maybe only Particle BU neccessary with different stencils (createStencil_X())
 * OR non virtual base class
 **********************************************************************************************/
#pragma once

#include <CAM/utils.hxx>
#include <algorithm>
#include <array>
#include <chrono>
#include <iostream>
#include <random>
#include <vector>
namespace CAM
{
/*!*********************************************************************************************
 * \brief Base struct template of building units (BU)
 * \param number index of cells in domain
 * \param jump_parameter How far BU is allowed to jump.
 * \param referenceFields cells in domain which is moved by CA and is calculation basis for all
 * cells in BU
 * \tparam nx
 **********************************************************************************************/
template <auto nx>
struct BuildingUnit
{
 public:
  virtual ~BuildingUnit(){};
  unsigned int number;
  unsigned int jump_parameter;
  std::vector<unsigned int> referenceFields, fields, fieldsBorder, stencilBorder, stencilBU;

  virtual const std::vector<unsigned int>& get_field_indices() = 0;
  virtual const std::vector<unsigned int>& get_border_indices() = 0;
  virtual bool is_member(unsigned int _index) = 0;
};
/*!*********************************************************************************************
 * \brief Custom BU with defined field indices in domain
 * Used for Single Cell BU
 * \TODO reference Point = fieldIndices[0] and stencil
 * \tparam nx
 **********************************************************************************************/
template <auto nx>
struct CustomBU : public BuildingUnit<nx>
{
  CustomBU(unsigned int _number,
           unsigned int _jump_parameter,
           std::vector<unsigned int> _fieldIndices)
  {
    this->number = _number;
    this->referenceFields = _fieldIndices;
    this->jump_parameter = _jump_parameter;
  }
  ~CustomBU() override {}

  const std::vector<unsigned int>& get_field_indices() override { return this->referenceFields; }
  const std::vector<unsigned int>& get_border_indices() override { return this->referenceFields; }
  bool is_member(unsigned int _index) override
  {
    return (std::find(this->referenceFields.begin(), this->referenceFields.end(), _index) !=
            this->referenceFields.end());
  }
};
/*!*********************************************************************************************
 * \brief Hyper sphere (2D: Circle, 3D: Sphere)
 *
 * \tparam nx
 * \param _position field index of center point; no position is given (-1) -> random position
 * \param _radius all cells are completely within the radius
 * \param _jump_parameter How far sphere is allowed to jump.
 **********************************************************************************************/
template <auto nx>
struct HyperSphereBU : public BuildingUnit<nx>
{
  double radius;

  HyperSphereBU(unsigned int _number,
                unsigned int _jump_parameter,
                unsigned int _centerPoint,
                double _radius)
  {
    this->number = _number;
    this->referenceFields = {_centerPoint};
    radius = _radius;
    this->jump_parameter = _jump_parameter;
    this->stencilBU = CAM::get_p_normed_particle<nx, 2>(radius);
    this->stencilBorder = CAM::get_border<nx>(this->stencilBU)[0];
  }
  ~HyperSphereBU() override {}

  const std::vector<unsigned int>& get_field_indices() override
  {
    this->fields.clear();
    for (unsigned int i = 0; i < this->stencilBU.size(); i++)
    {
      this->fields.push_back(CAM::aim<nx>(this->referenceFields[0], this->stencilBU[i]));
    }
    return this->fields;
  }
  const std::vector<unsigned int>& get_border_indices() override
  {
    this->fieldsBorder.clear();
    for (unsigned int i = 0; i < this->stencilBorder.size(); i++)
    {
      this->fieldsBorder.push_back(CAM::aim<nx>(this->referenceFields[0], this->stencilBorder[i]));
    }
    return this->fieldsBorder;
  }
  bool is_member(unsigned int _index) override
  {
    return (CAM::p_norm_distance<nx, 2>(this->referenceFields[0], (int)_index) < radius);
  }
};
/*!*********************************************************************************************
 * \brief Limited hyper plane (2D: Rectangle, 3D: cuboid)
 * \param _position field index of center point; no position is given (-1) -> random position
 * \param _extent size of plane in each dimension
 * \param _jump_parameter How far hyper plane is allowed to jump.
 *
 * \tparam nx
 **********************************************************************************************/
template <auto nx>
struct HyperPlaneBU : public BuildingUnit<nx>
{
  std::vector<unsigned int> extent;
  HyperPlaneBU(unsigned int _number,
               unsigned int _jump_parameter,
               unsigned int _referencePoint,
               std::vector<unsigned int> _extent)
  {
    this->number = _number;
    this->referenceFields = {_referencePoint};
    this->jump_parameter = _jump_parameter;
    extent = _extent;

    this->stencilBU.push_back(0);
    // static_assert(_extent.size() == nx.size(), "Size");
    unsigned int size, newMove;
    for (unsigned int dim = 0; dim < nx.size(); dim++)
    {
      size = this->stencilBU.size();
      for (unsigned int d = 1; d < extent[dim]; d++)
      {
        for (unsigned int i = 0; i < size; i++)
        {
          newMove = CAM::aim<nx>(this->stencilBU[i], d * direct_neigh<nx>(2 * dim + 1));
          this->stencilBU.push_back(newMove);
        }
      }
    }
    this->stencilBorder = CAM::get_border<nx>(this->stencilBU)[0];
  }
  ~HyperPlaneBU() override {}
  const std::vector<unsigned int>& get_field_indices() override
  {
    this->fields.clear();
    for (unsigned int i = 0; i < this->stencilBU.size(); i++)
    {
      this->fields.push_back(CAM::aim<nx>(this->referenceFields[0], this->stencilBU[i]));
    }
    return this->fields;
  }
  const std::vector<unsigned int>& get_border_indices() override
  {
    this->fieldsBorder.clear();
    for (unsigned int i = 0; i < this->stencilBorder.size(); i++)
    {
      this->fieldsBorder.push_back(CAM::aim<nx>(this->referenceFields[0], this->stencilBorder[i]));
    }
    return this->fieldsBorder;
  }
  bool is_member(unsigned int _index) override
  {
    // TODO calculate if its inside by considering extent/refPoints
    for (unsigned int i = 0; i < this->stencilBU.size(); i++)
    {
      if (CAM::aim<nx>(this->referenceFields[0], this->stencilBU[i]) == _index)
        return true;
    }
    return false;
  }
};
/*!*********************************************************************************************
 * \brief Particle with custom stencil
 * custom stencil defined by min and max feretDiameter and amount of extra cells attached to sphere
 * with radius min_feretDiameter/2
 *
 * \tparam nx
 **********************************************************************************************/
template <auto nx>
struct ParticleBU : public BuildingUnit<nx>
{
  ParticleBU(unsigned int _number,
             unsigned int _jump_parameter,
             unsigned int _seedPoint,
             std::vector<unsigned int> _stencil)
  {
    this->number = _number;
    this->referenceFields = {_seedPoint};
    this->jump_parameter = _jump_parameter;
    this->stencilBU = _stencil;
    this->stencilBorder = CAM::get_border<nx>(this->stencilBU)[0];
  }
  ~ParticleBU() override {}
  static const std::vector<unsigned int> get_stencil(const int _random_seed,
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
    std::vector<unsigned int> stencil = CAM::get_p_normed_particle<nx, 2>(_feretDiameter_min / 2.0);
    double feretDiameter;
    unsigned int neigh, coord_min, axis_min_feret = std::rand() % nx.size();
    std::vector<unsigned int> min_feret_coord, borderPoints, borderCells, points;
    std::array<std::vector<unsigned int>, 2> border_;
    std::array<unsigned int, n_fieldpoints<nx>()> newPoints;
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
      border_ = CAM::get_border<nx>(stencil);
      borderPoints = border_[1];
      borderCells = border_[0];
      points.reserve(borderPoints.size());
      for (unsigned int i = 0; i < borderPoints.size(); i++)
        points[i] = aim<nx>(centerpoint, borderPoints[i]);

      feretDiameter = CAM::feretDiameter_max<nx>(points);

      std::vector<unsigned int> newStencilCells;
      for (unsigned int a = 0; a < borderCells.size(); a++)
      {
        for (unsigned int i = 0; i < 2 * nx.size(); ++i)
        {
          neigh = aim<nx>(borderCells[a], direct_neigh<nx>(i));
          coord_min = (neigh / (int)direct_neigh<nx>(2 * axis_min_feret + 1) + n_fields<nx>()) %
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
        newPoints = CAM::get_points<nx>(newStencilCells[j]);
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

  const std::vector<unsigned int>& get_field_indices() override
  {
    this->fields.clear();
    for (unsigned int i = 0; i < this->stencilBU.size(); i++)
    {
      this->fields.push_back(CAM::aim<nx>(this->referenceFields[0], this->stencilBU[i]));
    }
    return this->fields;
  }
  const std::vector<unsigned int>& get_border_indices() override
  {
    this->fieldsBorder.clear();
    for (unsigned int i = 0; i < this->stencilBorder.size(); i++)
    {
      this->fieldsBorder.push_back(CAM::aim<nx>(this->referenceFields[0], this->stencilBorder[i]));
    }
    return this->fieldsBorder;
  }
  bool is_member(unsigned int _index) override
  {
    for (unsigned int i = 0; i < this->stencilBU.size(); i++)
    {
      if (CAM::aim<nx>(this->referenceFields[0], this->stencilBU[i]) == _index)
        return true;
    }
    return false;
  }
};
}  // namespace CAM
