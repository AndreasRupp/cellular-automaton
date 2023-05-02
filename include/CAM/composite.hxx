
#pragma once
#include <CAM/building_units.hxx>
#include <vector>
namespace CAM
{
/*!*********************************************************************************************
 * \brief Stores information about composites
 *
 * \tparam nx
 **********************************************************************************************/
template <auto nx>
struct Composite
{
  std::vector<CAM::BuildingUnit<nx>*> buildingUnits;
  unsigned int jump_parameter;
  std::vector<unsigned int> fieldIndices;
};

}  // namespace CAM
