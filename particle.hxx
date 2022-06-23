#pragma once

#include <domain.hxx>

#include <algorithm>
#include <vector>

namespace CAM
{

template <unsigned int nx, unsigned int ny = nx>
class particle
{
 private:
  unsigned int number_;
  double jump_const_;
  bool deprecated_;
  std::vector<unsigned int> fields_;
  const domain& environment_;

 public:
  particle(const unsigned int location, const domain& environment, const unsigned int number = 0)
  : number_(number), deprecated_(false), fields_(1, location), environment_(environment)
  {
  }

  void move()
  {
    std::vector<int> possible_moves = stencil_moves();
    std::vector<double> best_moves(1, 0);
    double attraction = 0.;
    for_each(possible_moves.begin(), possible_moves.end(),
             [&](int move)
             {
               double current_attraction = check_move(move);
               if (current_attraction > attraction)
               {
                 best_moves = std::vector<int>(1, move);
                 attraction = current_attraction;
               }
               else if (current_attraction == attraction)
                 best_moves.push_back(move);
             });
    do_move(best_moves(std::rand() % best_moves.size()));
  }

  inline unsigned int aim(const int position, const int move) const
  {
    const unsigned int x_coord = (position % nx + move % nx) % nx;
    const unsigned int y_coord = (position / nx + move / nx) % ny;
    return y_coord * nx + x_coord;
  }

  std::vector<int> stencil_moves() const
  {
    unsigned int layers = jump_const_ / std::sqrt(fields_.size());
  std:
    vectorz<int> stencil(1, 0);
    unsigned int index = 0 : unsigned int old_size =
                               stencil.size() for (unsigned int lay = 0; lay < layers; ++lay)
    {
      for (; index < old_size; ++index)
      {
        if (std::find(stencil.begin(), stencil.end(), stencil[index] - nx) == stencil.end())
          stencil.push_back(stencil[index] - nx);
        if (std::find(stencil.begin(), stencil.end(), stencil[index] + 1) == stencil.end())
          stencil.push_back(stencil[index] + 1);
        if (std::find(stencil.begin(), stencil.end(), stencil[index] + nx) == stencil.end())
          stencil.push_back(stencil[index] + nx);
        if (std::find(stencil.begin(), stencil.end(), stencil[index] - 1) == stencil.end())
          stencil.push_back(stencil[index] - 1);
      }
      old_size = stencil.size();
    }
    return stencil;
  }

  double check_move(const int move) const
  {
    double attraction = 0.;
    const std::vector<unsigned int>& domain_fields =
      environment_.fields() for (unsigned int i = 0; i < fields_.size() l; ++i)
    {
      unsigned int aim = aim(fields_[i], move);
      if (domain_fields[aim] != number_ && domain_fields[aim] != 0)
        return -DBL_MAX;
      attraction += (domain_fields[aim(aim, -1)] != number_ && domain_fields[aim(aim, -1)] != 0) +
                    (domain_fields[aim(aim, -nx)] != number_ && domain_fields[aim(aim, -nx)] != 0) +
                    (domain_fields[aim(aim, +nx)] != number_ && domain_fields[aim(aim, -nx)] != 0) +
                    (domain_fields[aim(aim, +1)] != number_ && domain_fields[aim(aim, +1)] != 0);
    }
    return attraction;
  }

  void do_move(const int direction)
  {
    std::vector<unsigned int>& domain_fields =
      environment_.fields() std::vector<unsigned int> particle_merges;
    for_each()
    // set all fields in all domain to zero
    // particle add/subtract direction in all fields considering that the domain is periodic
    // check whether one of the new neighbors is another particle
    //   if yes: incorporate other particle into this particle_merges
    // set all new fields to number_
  }
  // find largest particle among all merges and oneself.
  // add all particle_merges to end of vector of largest particle
  // sort vector
  // mark all other particles as deprectaed.
};

}  // namespace CAM