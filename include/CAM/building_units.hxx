#pragma once 

#include <vector>
#include <algorithm>
#include <memory>
namespace CAM
{
struct BuildingUnit
{  
    virtual ~BuildingUnit(){};
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
    ~ParticleBU()override{}
   
    std::vector<unsigned int> getFieldIndices() override
    {
        return referencePoints;
    }
    bool isMember(unsigned int _index) override
    {
        return (std::find(referencePoints.begin(), referencePoints.end(), _index) != referencePoints.end());
    }
};
struct SphereBU : public BuildingUnit
{
    double radius;
    SphereBU(unsigned int _number, unsigned int _jump_parameter, unsigned int _centerPoint, double _radius ) : BuildingUnit(){
    number = _number;
    referencePoints = {_centerPoint};
    radius = _radius;
    jump_parameter = _jump_parameter;
        
    }
    ~SphereBU()override{}
   
    std::vector<unsigned int> getFieldIndices() override
    {
        //TODO calculate Sphere
        return {1};
    }
    bool isMember(unsigned int _index) override
    {
        return true;
    }
};

}
