/*!*********************************************************************************************
 * \file building_units.hxx
 * \brief Implementation of different building units
 * Custombu
 * HyperPlane
 * HyperSphere
 **********************************************************************************************/
#pragma once

#include <CAM/utils.hxx>
#include <algorithm>
#include <array>
#include <chrono>
#include <exception>
#include <iostream>
#include <random>
#include <vector>
namespace CAM
{
/*!*********************************************************************************************
 * \brief Class of building units (bu)
 * \param number index of cells in domain
 * \param jump_parameter How far bu is allowed to jump.
 * \param reference_field cells in domain which is moved by CA and is calculation basis for all
 *cells in bu  (for sphere center_point) \param shape shape of bu \param boundary shape boundary of
 *bu \tparam nx
 **********************************************************************************************/
template <auto nx>
class BuildingUnit
{
 private:
  double jump_parameter;
  std::vector<unsigned int> shape, boundary;
  unsigned int reference_field, number;

 public:
  /*!*********************************************************************************************
   * \brief Construct a new Building Unit object
   ************************************************************************************************/
  BuildingUnit(const double _jump_parameter,
               const std::vector<unsigned int>& _shape,
               const unsigned int _reference_field,
               const unsigned int _number)
  : jump_parameter(_jump_parameter),
    shape(_shape),
    reference_field(_reference_field),
    number(_number)
  {
    boundary = CAM::get_boundary_fields<nx>(shape);
  }

  void set_reference_field(const unsigned int _reference_field)
  {
    reference_field = _reference_field;
  }
  unsigned int get_reference_field() const { return reference_field; }
  unsigned int get_number() const { return number; }
  unsigned int get_jump_parameter() const { return jump_parameter; }
  const std::vector<unsigned int>& get_shape() const { return shape; }
  const std::vector<unsigned int>& get_boundary() const { return boundary; }

  bool constexpr is_member(const unsigned int _index)
  {
    for (unsigned int i = 0; i < shape.size(); i++)
    {
      if (CAM::aim<nx>(reference_field, shape[i]) == _index)
        return true;
    }
    return false;
  }
};
/*!*********************************************************************************************
 * \brief Hyper sphere (2D: Circle, 3D: Sphere)
 * \tparam nx
 * \param _radius all cells of bu are completely within the radius
 **********************************************************************************************/
template <auto nx>
static CAM::BuildingUnit<nx> create_hyper_sphere(const double _jump_parameter,
                                                 const double _radius,
                                                 const int _reference_field = -1,
                                                 const int _number = -1)
{
  std::vector<unsigned int> shape = CAM::get_p_normed_shape<nx, 2>(_radius);

  unsigned int reference_field, number;
  if (_reference_field < 0)
    reference_field = CAM::get_random_field_index<nx>();
  else
    reference_field = _reference_field;
  // or what should be done with non valid numbers?
  if (_number <= 0)
    number = CAM::get_random_field_index<nx>();
  else
    number = _number;

  return CAM::BuildingUnit<nx>(_jump_parameter, shape, reference_field, number);
}
/*!*********************************************************************************************
 * \brief Limited hyper plane (2D: Rectangle, 3D: cuboid)
 * \param _extent size of plane in each dimension
 * \tparam nx
 **********************************************************************************************/
template <auto nx>
static CAM::BuildingUnit<nx> create_hyper_plane(const double _jump_parameter,
                                                const std::array<unsigned int, nx.size()>& _extent,
                                                const int _reference_field = -1,
                                                const int _number = -1)
{
  unsigned int size, new_move;
  std::vector<unsigned int> shape = {0};
  for (unsigned int dim = 0; dim < nx.size(); dim++)
  {
    size = shape.size();
    for (unsigned int d = 1; d < _extent[dim]; d++)
    {
      for (unsigned int i = 0; i < size; i++)
      {
        new_move = CAM::aim<nx>(shape[i], d * direct_neigh<nx>(2 * dim + 1));
        shape.push_back(new_move);
      }
    }
  }

  unsigned int reference_field, number;
  if (_reference_field < 0)
    reference_field = CAM::get_random_field_index<nx>();
  else
    reference_field = _reference_field;
  if (_number <= 0)
    number = CAM::get_random_field_index<nx>();
  else
    number = _number;

  return CAM::BuildingUnit<nx>(_jump_parameter, shape, reference_field, number);
}

/*!*********************************************************************************************
 * \brief Bu defined by custom shape
 * \param shape custom shape
 * \tparam nx
 **********************************************************************************************/
template <auto nx>
static CAM::BuildingUnit<nx> create_custom_bu(const double _jump_parameter,
                                              const std::vector<unsigned int>& _shape,
                                              const int _reference_field = -1,
                                              const int _number = -1)
{
  unsigned int reference_field, number;
  if (_reference_field < 0)
    reference_field = CAM::get_random_field_index<nx>();
  else
    reference_field = _reference_field;
  if (_number <= 0)
    number = CAM::get_random_field_index<nx>();
  else
    number = _number;

  return CAM::BuildingUnit<nx>(_jump_parameter, _shape, reference_field, number);
}
/*!*********************************************************************************************
 * \brief bu containing only one cell
 * \tparam nx
 **********************************************************************************************/
template <auto nx>
static constexpr CAM::BuildingUnit<nx> create_single_cell_bu(const double _jump_parameter,
                                                             const int _reference_field = -1,
                                                             const int _number = -1)
{
  unsigned int reference_field, number;
  if (_reference_field < 0)
    reference_field = CAM::get_random_field_index<nx>();
  else
    reference_field = _reference_field;
  if (_number <= 0)
    number = CAM::get_random_field_index<nx>();
  else
    number = _number;

  const std::vector<unsigned int> shape = {0};
  return CAM::BuildingUnit<nx>(_jump_parameter, shape, reference_field, number);
}

/*!*********************************************************************************************
 * \brief
 * custom shape defined by min and max feret_diameter and amount of extra cells attached to
 *sphere with radius min_feret_diameter/2
 *
 * \tparam nx
 **********************************************************************************************/
template <auto nx>
static constexpr std::vector<unsigned int> get_shape_by_feret(const int _random_seed,
                                                              const double _feret_diameter_max,
                                                              const double _feret_diameter_min,
                                                              const unsigned int _extra_cells)
{
  unsigned int rand_seed;
  if (_random_seed == 0)
  {
    rand_seed = std::chrono::system_clock::now().time_since_epoch().count();
  }
  else
    rand_seed = _random_seed;
  std::srand(rand_seed);

  unsigned int center_point = 0;
  for (unsigned int i = 0; i < nx.size(); ++i)
  {
    center_point += nx[i] / 2 * direct_neigh<nx>(2 * i + 1);
  }
  std::vector<unsigned int> shape = CAM::get_p_normed_shape<nx, 2>(_feret_diameter_min / 2.0);
  double feret_diameter;
  unsigned int neigh, coord_min, axis_min_feret = std::rand() % nx.size();
  std::vector<unsigned int> min_feret_coord, boundary_points, boundary_cells, points;
  std::array<std::vector<unsigned int>, 2> boundary_;
  std::array<unsigned int, n_field_corner_points<nx>()> new_points;
  for (unsigned int i = 0; i < shape.size(); i++)
  {
    min_feret_coord.push_back(
      (shape[i] / (int)direct_neigh<nx>(2 * axis_min_feret + 1) + n_fields<nx>()) %
      nx[axis_min_feret]);
  }
  sort(min_feret_coord.begin(), min_feret_coord.end());
  min_feret_coord.erase(unique(min_feret_coord.begin(), min_feret_coord.end()),
                        min_feret_coord.end());

  for (unsigned int c = 0; c < _extra_cells; c++)
  {
    boundary_ = CAM::get_boundary<nx>(shape);
    boundary_points = boundary_[1];
    boundary_cells = boundary_[0];
    points.reserve(boundary_points.size());
    for (unsigned int i = 0; i < boundary_points.size(); i++)
      points[i] = aim<nx>(center_point, boundary_points[i]);

    feret_diameter = CAM::feret_diameter_max<nx>(points);

    std::vector<unsigned int> new_shape_cells;
    for (unsigned int a = 0; a < boundary_cells.size(); a++)
    {
      for (unsigned int i = 0; i < 2 * nx.size(); ++i)
      {
        neigh = aim<nx>(boundary_cells[a], direct_neigh<nx>(i));
        coord_min = (neigh / (int)direct_neigh<nx>(2 * axis_min_feret + 1) + n_fields<nx>()) %
                    nx[axis_min_feret];
        if (std::find(shape.begin(), shape.end(), neigh) == shape.end() &&
            find(min_feret_coord.begin(), min_feret_coord.end(), coord_min) !=
              min_feret_coord.end())
        {
          new_shape_cells.push_back(neigh);
        }
      }
    }
    // here possible other metrics as smoothness
    std::shuffle(new_shape_cells.begin(), new_shape_cells.end(),
                 std::default_random_engine(std::rand()));
    for (unsigned int j = 0; j < new_shape_cells.size(); j++)
    {
      new_points = CAM::get_corner_points<nx>(new_shape_cells[j]);
      for (unsigned int p = 0; p < new_points.size(); p++)
        points.push_back(aim<nx>(center_point, (int)new_points[p]));
      feret_diameter = CAM::feret_diameter_max<nx>(points);

      if (feret_diameter < _feret_diameter_max)
      {
        shape.push_back(new_shape_cells[j]);
        break;
      }
      else
      {
        for (unsigned int p = 0; p < new_points.size(); p++)
          points.pop_back();
      }
    }
  }
  return shape;
}

}  // namespace CAM
