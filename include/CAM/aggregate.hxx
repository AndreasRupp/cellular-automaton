#pragma once
#include <CAM/building_units.hxx>
#include <vector>
namespace CAM
{
struct Aggregate
{
  std::vector<CAM::BuildingUnit*> buildingUnits;
  unsigned int jump_parameter;
  // is not updated -> implement fct getFieldIndices iterating over BUs
  std::vector<unsigned int> fieldIndices;
};

}  // namespace CAM