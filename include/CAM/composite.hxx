
#pragma once
#include <CAM/building_units.hxx>
#include <vector>
namespace CAM
{
static double jump_parameter_composite;
template <auto nx>
static constexpr double get_jump_range_composite(const unsigned int _comp_size)
{
  return jump_parameter_composite / std::pow(_comp_size, 1.0 / (double)nx.size());
}
/*!*********************************************************************************************
 * \brief Stores information about composites
 *
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
