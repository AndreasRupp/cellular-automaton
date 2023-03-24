#pragma once

#include <algorithm>
#include <memory>
#include <vector>
#include <array>
#include <CAM/utils.hxx>
namespace CAM
{
struct BuildingUnit
{
  virtual ~BuildingUnit(){};
  CAM::fieldNumbers number;
  unsigned int jump_parameter;
  std::vector<unsigned int> referencePoints;

  virtual std::vector<unsigned int> getFieldIndices() = 0;

  virtual bool isMember(unsigned int _index) = 0;
};
struct ParticleBU : public BuildingUnit
{
  ParticleBU(unsigned int _number,
             unsigned int _jump_parameter,
             std::vector<unsigned int> _fieldIndices)
  : BuildingUnit()
  {
    number = _number;
    referencePoints = _fieldIndices;
    jump_parameter = _jump_parameter;
  }
  ~ParticleBU() override {}

  std::vector<unsigned int> getFieldIndices() override { return referencePoints; }
  bool isMember(unsigned int _index) override
  {
    return (std::find(referencePoints.begin(), referencePoints.end(), _index) !=
            referencePoints.end());
  }
};
struct HyperSphereBU : public BuildingUnit
{
  double radius;
  std::vector<int> stencilBU;
  static constexpr std::array<unsigned int, 2> nx = {10, 10};

  HyperSphereBU(unsigned int _number,
           unsigned int _jump_parameter,
           unsigned int _centerPoint,
           double _radius)
  : BuildingUnit()
  {
    number = _number;
    referencePoints = {_centerPoint};
    radius = _radius;
    jump_parameter = _jump_parameter;
    stencilBU = CAM::getPNormedStencil<nx>(radius, 2);//CAM::getPNormedStencil<nx>(radius, 2);
  }
  ~HyperSphereBU() override {}

  std::vector<unsigned int> getFieldIndices() override
  {
    std::vector<unsigned int> fields;
   // std::cout<<"ref "<<referencePoints[0]<<std::endl;
    for(unsigned int i = 0; i < stencilBU.size(); i++)
    {
      
      fields.push_back(CAM::aim<nx>(referencePoints[0], stencilBU[i]));
      //std::cout<<fields[i]<< " ";
    }
    //std::cout<<std::endl;
    // std::cout<< "fie_siz "<<fields.size()<<std::endl;
    // TODO calculate Sphere
    return fields;
  }
  bool isMember(unsigned int _index) override { return _index == 42; }
};

}  // namespace CAM
