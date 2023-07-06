
#pragma once
#include <CAM/building_units.hxx>
#include <cmath>
#include <vector>
namespace CAM
{
static double jump_parameter_composite;
template <auto nx>
static constexpr double get_jump_range_composite(const unsigned int _comp_size)
{
  return std::floor(
    jump_parameter_composite / std::sqrt(_comp_size) +
    0.5);  // jump_parameter_composite;  // / std::pow(_comp_size, 1.0 / (double)nx.size());
}
/*!*********************************************************************************************
 * \brief Stores information about composites
 * \param building_units what bus is the composite made of
 * \param jump_parameter How far composite is allowed to jump.
 * \param field_indices what field indices is the composite made of (information also stored in bus)
 * \tparam nx
 **********************************************************************************************/
template <auto nx>
struct Composite
{
  std::vector<CAM::BuildingUnit<nx>*> building_units;
  double jump_parameter;
  std::vector<unsigned int> field_indices;
};

}  // namespace CAM
