#pragma once

#include <particle.hxx>

#include <array>
#include <vector>

namespace CAM
{

template <unsigned int nx, unsigned int ny>
class domain
{
 private:
  std::array<nx * ny, unsigned int> fields_;
  std::vector<CAM::particle<nx, ny>> particles_;
  std::vector<unsigned int> particle_list_;
  domain(const unsigned int nx, const unsigned int ny) { fields_.fill(0); }

 public:
  std::vector<unsigned int>& fields() { return fields; }
};

}  // namespace CAM
