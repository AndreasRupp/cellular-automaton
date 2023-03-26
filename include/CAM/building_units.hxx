#pragma once

#include <CAM/utils.hxx>
#include <algorithm>
#include <array>
#include <memory>
#include <vector>
namespace CAM
{
template <auto nx>
struct BuildingUnit
{
 public:
  virtual ~BuildingUnit(){};
  CAM::fieldNumbers_t number;
  unsigned int jump_parameter;
  std::vector<unsigned int> referencePoints;

  virtual std::vector<unsigned int> getFieldIndices() = 0;

  virtual bool isMember(unsigned int _index) = 0;
};
template <auto nx>
struct ParticleBU : public BuildingUnit<nx>
{
  ParticleBU(unsigned int _number,
             unsigned int _jump_parameter,
             std::vector<unsigned int> _fieldIndices)
  {
    this->number = _number;
    this->referencePoints = _fieldIndices;
    this->jump_parameter = _jump_parameter;
  }
  ~ParticleBU() override {}

  std::vector<unsigned int> getFieldIndices() override { return this->referencePoints; }
  bool isMember(unsigned int _index) override
  {
    return (std::find(this->referencePoints.begin(), this->referencePoints.end(), _index) !=
            this->referencePoints.end());
  }
};
template <auto nx>
struct HyperSphereBU : public BuildingUnit<nx>
{
  double radius;
  std::vector<int> stencilBU;

  HyperSphereBU(unsigned int _number,
                unsigned int _jump_parameter,
                unsigned int _centerPoint,
                double _radius)
  {
    this->number = _number;
    this->referencePoints = {_centerPoint};
    radius = _radius;
    this->jump_parameter = _jump_parameter;
    stencilBU = CAM::getPNormedStencil<nx>(radius, 2);
  }
  ~HyperSphereBU() override {}

  std::vector<unsigned int> getFieldIndices() override
  {
    std::vector<unsigned int> fields;
    for (unsigned int i = 0; i < stencilBU.size(); i++)
    {
      fields.push_back(CAM::aim<nx>(this->referencePoints[0], stencilBU[i]));
    }
    return fields;
  }
  bool isMember(unsigned int _index) override { return _index == 42; }
};

}  // namespace CAM
