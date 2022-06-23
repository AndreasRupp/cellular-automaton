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
  static constexpr std::array<int, 4> direct_neigh_ = {-nx, 1, nx, -1};
  const double jump_const_;
  unsigned int number_;
  bool deprecated_;
  std::vector<unsigned int> fields_;
  domain<nx, ny>& environment_;

 public:
  particle(const unsigned int location,
           const domain& environment,
           const double jump_const = 5,
           const unsigned int number = 0)
  : jump_const_(jump_const),
    number_(number),
    deprecated_(false),
    fields_(1, location),
    environment_(environment)
  {
  }

  bool is_depreacted() { return deprecated_; }

  inline unsigned int aim(const int position, const int move) const
  {
    const unsigned int x_coord = (nx + position % nx + move % nx) % nx;
    const unsigned int y_coord = (ny + position / nx + move / nx) % ny;
    return y_coord * nx + x_coord;
  }

  void move_all()
  {
    std::vector<int> possible_moves = stencil_moves();
    std::vector<double> best_moves(1, 0);
    double attraction = 0.;
    std::for_each(possible_moves.begin(), possible_moves.end(),
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

  std::vector<int> stencil_moves_all() const
  {
    unsigned int layers = std::max(1, jump_const_ / std::sqrt(fields_.size()));
    std::vector<int> stencil(1, 0);
    unsigned int index = 0;
    unsigned int old_size = stencil.size();

    for (unsigned int lay = 0; lay < layers; ++lay)
    {
      for (; index < old_size; ++index)
        for (unsigned int i = 0; i < 4; ++i)
          if (std::find(stencil.begin(), stencil.end(), stencil[index] + direct_neigh_[i]) ==
              stencil.end())
            stencil.push_back(stencil[index] + direct_neigh_[i]);
      old_size = stencil.size();
    }
    return stencil;
  }

  double check_move_all(const int move) const
  {
    double attraction = 0.;
    const std::vector<unsigned int>& domain_fields =
      environment_.fields() for (unsigned int i = 0; i < fields_.size() l; ++i)
    {
      unsigned int aim = aim(fields_[i], move);
      if (domain_fields[aim] != number_ && domain_fields[aim] != 0)
        return -DBL_MAX;
      for (unsigned int i = 0; i < 4; ++i)
        attraction += domain_fields[aim(aim, direct_neigh_[i])] != number_ &&
                      domain_fields[aim(aim, direct_neigh_[i])] != 0;
    }
    return attraction;
  }

  void do_move_all(const int move)
  {
    std::vector<unsigned int>& domain_fields = environment_.fields();
    std::vector<unsigned int> particle_merges;
    std::for_each(
      fields_.begin(), fields_.end(),
      [&](unsigned int index)
      {
        domain_fields[fields_[index]] -= number_;
        unsigned int aim = aim(fields_[index], move);
        domain_fields[aim] += number_;

        for (unsigned int i = 0; i < 4; ++i)
          if (domain_fields[aim(aim, direct_neigh_[i])] != number_ &&
              domain_fields[aim(aim, direct_neigh_[i])] != 0 &&
              std::find(particle_merges.begin(), particle_merges.end(),
                        domain_fields[aim(aim, direct_neigh_[i])]) == particle_merges.end())
            particle_merges.push_back(domain_fields[aim(aim, direct_neigh_[i])]);
      });
    if (particle_merges.empty())
      return;
    // unsigned int max_ind = std::distance(
    //   particle_merges.begin(), std::max_element(particle_merges.begin(), particle_merges.end()));
    std::vector<particle<nx, ny> > partciles = environment_.particles();
    std::for_each(
      particle_merges.begin(), particle_merges.end()[&](unsigned int index) {
        auto other_particle = std::find(particles.begin(), particles.end(), index);
        std::vector<unsigned int>& other_fields = other_particle->fields_;
        for_each(other_fields.begin(), other_fields.end(),
                 [&](unsigned int ind) { domain_fields[ind] = number_; });
        particle->deprecated_ = true;

        fields_.insert(fields_.end(), other_fields.begin(), other_fields.end());
      });
  }

  bool operator==(const particle<nx, ny>& other) const { return number_ == other.number_; }

  bool operator<(const particle<nx, ny>& other) const
  {
    return fields_.size() < other.fields_.size();
  }
};

}  // namespace CAM