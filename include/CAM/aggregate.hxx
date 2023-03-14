#pragma once
#include <vector>
#include <CAM/building_units.hxx>
namespace CAM
{
    struct Aggregate
    {
        std::vector<CAM::BuildingUnit*> buildingUnits;
        unsigned int jump_parameter;
        //is not updated -> implement fct getFieldIndices iterating over BUs
        std::vector<unsigned int> fieldIndices;     
    };

}