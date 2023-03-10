#pragma once 

#include <vector>
#include <algorithm>
#include <memory>
namespace CAM
{
struct BuildingUnit
{  
    ~BuildingUnit(){};
    unsigned int number;
    unsigned int jump_parameter;
    std::vector<unsigned int> referencePoints;

    virtual std::vector<unsigned int> getFieldIndices() = 0;


    virtual bool isMember(unsigned int _index) = 0;
};
struct ParticleBU : public BuildingUnit
{
    ParticleBU(unsigned int _number, unsigned int _jump_parameter, std::vector<unsigned int> _fieldIndices ) : BuildingUnit(){
        number = _number;
        referencePoints = _fieldIndices;
        jump_parameter = _jump_parameter;
        
    }
    ~ParticleBU(){}
   
    std::vector<unsigned int> getFieldIndices() override
    {
        return referencePoints;
    }
    bool isMember(unsigned int _index) override
    {
        return (std::find(referencePoints.begin(), referencePoints.end(), _index) != referencePoints.end());
    }
};
// struct Aggregates : public BuildingUnit
// {
//         Aggregates(unsigned int _number, vector<unsigned int>  unsigned int _jump_parameter, std::vector<unsigned int> _fieldIndices ) : BuildingUnit(){
//         number = _number;
//         referencePoints = _fieldIndices;
//         jump_parameter = _jump_parameter;
        
//     }
// };
}
