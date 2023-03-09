#pragma once 

#include <vector>
#include <algorithm>
#include <memory>
#include <CAM/domain_new.hxx>
namespace CAM
{
//template<auto nx>
struct BuildingUnit
{
    
    unsigned int number;
    unsigned int jump_parameter;
    virtual std::vector<unsigned int> getFieldIndices() = 0;
    void setFieldIndices(){};
    virtual bool isMember(unsigned int _index) = 0;
    //virtual void move(const int _move) = 0;


};
struct ParticleBU : public BuildingUnit
{
     std::vector<unsigned int> fieldIndices;
    ParticleBU(){}
    ParticleBU(unsigned int _number, std::vector<unsigned int> _fieldIndices ) : BuildingUnit(){
        number = _number;
        fieldIndices = _fieldIndices;
    }
   
    std::vector<unsigned int> getFieldIndices() override 
    {
        return fieldIndices;
    }
    bool isMember(unsigned int _index) override
    {
        return (std::find(fieldIndices.begin(), fieldIndices.end(), _index) != fieldIndices.end());
    }
    void setFieldIndices(std::vector<unsigned int> _fieldIndices) 
    {
        fieldIndices = _fieldIndices;

    }
    // template<auto nx>
    // void move(const int _move)
    // {
    //     std::cout<<"move"<<std::endl;
    //     std::for_each(fieldIndices.begin(), fieldIndices.end(), [&](unsigned int& field)
    //     {
    //         field = aim<nx>(field, _move);     
    //     });
    // }
};

}
