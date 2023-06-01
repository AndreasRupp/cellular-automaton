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

template<auto nx>
static BuildingUnit<nx> create_hyper_sphere(const double jump_parameter, const double _radius, const int number=-1, const int center=-1)
{
  std::vector<unsigned int> shape = ...
  if (center == -1)
    center = // random number in [0,prod(nx)]
  return BuildingUnit<nx>(number, jump_parameter, center, shape);
}

/*!*********************************************************************************************
 * \brief Base struct template of building units (bu)
 * \param number index of cells in domain
 * \param jump_parameter How far bu is allowed to jump.
 * \param reference_field cells in domain which is moved by CA and is calculation basis for all (for
 *sphere center_point) cells in bu \tparam nx
 * \param stencil_bu stencil of all cells of bu
 * \param stencil_bu stencil of all cells at the border of bu
 * \param fields all cells indices in domain covered by bu
 * \param fields_border all cells in domain covered by border of bu
 **********************************************************************************************/
template <auto nx>
class BuildingUnit
{
 private:
  unsigned int number, reference_field;
  double jump_parameter;
  // std::vector<unsigned int> stencil_bu, stencil_border, fields, fields_border;
  std::vector<unsigned int> shape, boundary;
public:
  BuildingUnit(const unsigned int number, const unsigned int center, const double jump_parameter, const std::vector<unsigned int>& shape)

  void set_center(const unsigned int _center)
  {
    center = _center;
  }

  // BuildingUnit();
  BuildingUnit(double _jump_parameter) : jump_parameter(_jump_parameter)
  {
    stencil_bu.push_back(0);
    stencil_border = CAM::get_border<nx>(stencil_bu)[0];
  }

  /*!*********************************************************************************************
   * \brief Hyper sphere (2D: Circle, 3D: Sphere)
   *
   * \tparam nx
   * \param _radius all cells are completely within the radius
   **********************************************************************************************/
  static constexpr BuildingUnit<nx> create_hyper_sphere(double const jump_parameter,
                                                        const double _radius);
  /*!*********************************************************************************************
   * \brief Limited hyper plane (2D: Rectangle, 3D: cuboid)
   * \param _extent size of plane in each dimension
   *
   * \tparam nx
   **********************************************************************************************/
  static constexpr BuildingUnit<nx> create_hyper_plane(
    double const jump_parameter,
    const std::array<unsigned int, nx.size()>& _extent);
  /*!*********************************************************************************************
   * \brief Custom bu defined by custom stencil
   * \param stencil custom stencil
   * \tparam nx
   **********************************************************************************************/
  static constexpr BuildingUnit<nx> create_custom_bu(double const jump_parameter,
                                                     const std::vector<unsigned int>& _stencil);
  /*!*********************************************************************************************
   * \brief bu containing only one cell
   * \tparam nx
   **********************************************************************************************/
  static constexpr BuildingUnit<nx> create_single_cell_bu(double const jump_parameter)
  {
    return BuildingUnit<nx>(jump_parameter);
  }
  /*!*********************************************************************************************
   * \brief
   * custom stencil defined by min and max feret_diameter and amount of extra cells attached to
   *sphere with radius min_feret_diameter/2
   *
   * \tparam nx
   **********************************************************************************************/
  static constexpr std::vector<unsigned int> get_custom_particle_stencil(
    const int _random_seed,
    const double _feret_diameter_max,
    const double _feret_diameter_min,
    const unsigned int _extra_cells);
  const std::vector<unsigned int>& get_field_indices();
  const std::vector<unsigned int>& get_border_indices();
  bool is_member(unsigned int _index);
};

template <auto nx>
const std::vector<unsigned int>& BuildingUnit<nx>::get_field_indices()
{
  fields.clear();
  for (unsigned int i = 0; i < stencil_bu.size(); i++)
  {
    fields.push_back(CAM::aim<nx>(reference_field, stencil_bu[i]));
  }
  return fields;
}
template <auto nx>
const std::vector<unsigned int>& BuildingUnit<nx>::get_border_indices()
{
  fields_border.clear();
  for (unsigned int i = 0; i < stencil_border.size(); i++)
  {
    fields_border.push_back(CAM::aim<nx>(reference_field, stencil_border[i]));
  }
  return fields_border;
}
template <auto nx>
bool BuildingUnit<nx>::is_member(const unsigned int _index)
{
  for (unsigned int i = 0; i < stencil_bu.size(); i++)
  {
    if (CAM::aim<nx>(reference_field, stencil_bu[i]) == _index)
      return true;
  }
  return false;
}
template <auto nx>
constexpr BuildingUnit<nx> BuildingUnit<nx>::create_hyper_sphere(double const jump_parameter,
                                                                 const double _radius)
{
  BuildingUnit<nx> bu = BuildingUnit<nx>(jump_parameter);
  bu.stencil_border.clear();
  bu.stencil_bu.clear();
  bu.stencil_bu = CAM::get_p_normed_particle<nx, 2>(_radius);
  bu.stencil_border = CAM::get_border<nx>(bu.stencil_bu)[0];
  return bu;
}
template <auto nx>
constexpr BuildingUnit<nx> BuildingUnit<nx>::create_hyper_plane(
  double const jump_parameter,
  const std::array<unsigned int, nx.size()>& _extent)
{
  BuildingUnit<nx> bu = BuildingUnit<nx>(jump_parameter);
  bu.stencil_border.clear();
  bu.stencil_bu.clear();
  bu.stencil_bu.push_back(0);
  unsigned int size, new_move;
  for (unsigned int dim = 0; dim < nx.size(); dim++)
  {
    size = bu.stencil_bu.size();
    for (unsigned int d = 1; d < _extent[dim]; d++)
    {
      for (unsigned int i = 0; i < size; i++)
      {
        new_move = CAM::aim<nx>(bu.stencil_bu[i], d * direct_neigh<nx>(2 * dim + 1));
        bu.stencil_bu.push_back(new_move);
      }
    }
  }
  bu.stencil_border = CAM::get_border<nx>(bu.stencil_bu)[0];
  return bu;
}
template <auto nx>
constexpr BuildingUnit<nx> BuildingUnit<nx>::create_custom_bu(
  double const jump_parameter,
  const std::vector<unsigned int>& _stencil)
{
  BuildingUnit<nx> bu = BuildingUnit<nx>(jump_parameter);
  bu.stencil_bu = _stencil;
  bu.stencil_border = CAM::get_border<nx>(bu.stencil_bu)[0];
  return bu;
}
template <auto nx>
constexpr std::vector<unsigned int> BuildingUnit<nx>::get_custom_particle_stencil(
  const int _random_seed,
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
  std::vector<unsigned int> stencil = CAM::get_p_normed_particle<nx, 2>(_feret_diameter_min / 2.0);
  double feret_diameter;
  unsigned int neigh, coord_min, axis_min_feret = std::rand() % nx.size();
  std::vector<unsigned int> min_feret_coord, border_points, border_cells, points;
  std::array<std::vector<unsigned int>, 2> border_;
  std::array<unsigned int, n_field_corner_points<nx>()> new_points;
  for (unsigned int i = 0; i < stencil.size(); i++)
  {
    min_feret_coord.push_back(
      (stencil[i] / (int)direct_neigh<nx>(2 * axis_min_feret + 1) + n_fields<nx>()) %
      nx[axis_min_feret]);
  }
  sort(min_feret_coord.begin(), min_feret_coord.end());
  min_feret_coord.erase(unique(min_feret_coord.begin(), min_feret_coord.end()),
                        min_feret_coord.end());

  for (unsigned int c = 0; c < _extra_cells; c++)
  {
    border_ = CAM::get_border<nx>(stencil);
    border_points = border_[1];
    border_cells = border_[0];
    points.reserve(border_points.size());
    for (unsigned int i = 0; i < border_points.size(); i++)
      points[i] = aim<nx>(center_point, border_points[i]);

    feret_diameter = CAM::feret_diameter_max<nx>(points);

    std::vector<unsigned int> new_stencil_cells;
    for (unsigned int a = 0; a < border_cells.size(); a++)
    {
      for (unsigned int i = 0; i < 2 * nx.size(); ++i)
      {
        neigh = aim<nx>(border_cells[a], direct_neigh<nx>(i));
        coord_min = (neigh / (int)direct_neigh<nx>(2 * axis_min_feret + 1) + n_fields<nx>()) %
                    nx[axis_min_feret];
        if (std::find(stencil.begin(), stencil.end(), neigh) == stencil.end() &&
            find(min_feret_coord.begin(), min_feret_coord.end(), coord_min) !=
              min_feret_coord.end())
        {
          new_stencil_cells.push_back(neigh);
        }
      }
    }
    // here possible other metrics as smoothness
    std::shuffle(new_stencil_cells.begin(), new_stencil_cells.end(),
                 std::default_random_engine(std::rand()));
    for (unsigned int j = 0; j < new_stencil_cells.size(); j++)
    {
      new_points = CAM::get_corner_points<nx>(new_stencil_cells[j]);
      for (unsigned int p = 0; p < new_points.size(); p++)
        points.push_back(aim<nx>(center_point, (int)new_points[p]));
      feret_diameter = CAM::feret_diameter_max<nx>(points);

      if (feret_diameter < _feret_diameter_max)
      {
        stencil.push_back(new_stencil_cells[j]);
        break;
      }
      else
      {
        for (unsigned int p = 0; p < new_points.size(); p++)
          points.pop_back();
      }
    }
  }
  return stencil;
}

}  // namespace CAM
