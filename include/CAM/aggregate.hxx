#pragma once
#include <vector>
#include <CAM/building_units.hxx>
namespace CAM
{
    struct Aggregate
    {
        std::vector<CAM::BuildingUnit*> buildingUnits;
        unsigned int jump_parameter;
        std::vector<unsigned int> fieldIndices;
        // std::vector<unsigned int> getFieldIndices()
        // {
        //     std::vector<unsigned int> fieldIndices;
        //     std::for_each(buildingUnits.begin(), buildingUnits.end(), [&](BuildingUnit* unit)
        //     {
        //         fieldIndices.push_back(unit->getFieldIndices());
        //     });
        //     return fieldIndices;
        // }
        
    };

}