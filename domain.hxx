#pragma once

#include <algorithm>
#include <array>
#include <cfloat>
#include <cmath>
#include <random>
#include <vector>

template <unsigned int nx, unsigned int ny = nx>
class cellular_automaton
{
 private:
  // -----------------------------------------------------------------------------------------------

  class particle
  {
   private:
    static constexpr std::array<int, 4> direct_neigh_ = {-1 * (int)nx, 1, nx, -1};
    unsigned int number_;
    bool deprecated_;
    std::vector<unsigned int> fields_;
    cellular_automaton<nx, ny>& domain_;

   public:
    particle(const unsigned int location,
             cellular_automaton<nx, ny>& domain,
             const unsigned int number)
    : number_(number), deprecated_(false), fields_(1, location), domain_(domain)
    {
    }

    particle(particle&& other) noexcept
    : number_(other.number_),
      deprecated_(other.deprecated_),
      fields_(std::move(other.fields_)),
      domain_(other.domain_)
    {
    }

    particle& operator=(particle&& other) noexcept
    {
      number_ = other.number_;
      deprecated_ = other.deprecated_;
      std::swap(fields_, other.fields_);
      // domain_ = other.domain_;
      return *this;
    }

    bool is_deprecated() const { return deprecated_; }

    inline unsigned int aim(const int position, const int move) const
    {
      const unsigned int x_coord = (position % nx + (nx + move) % nx) % nx;
      const unsigned int y_coord = (position / nx + (nx + move) / nx - 1 + ny) % ny;
      return y_coord * nx + x_coord;
    }

    void move_all()
    {
      std::vector<int> possible_moves = stencil_moves_all();
      std::vector<int> best_moves(1, 0);
      double attraction = 0.;
      std::for_each(possible_moves.begin(), possible_moves.end(),
                    [&](int move)
                    {
                      double current_attraction = check_move_all(move);
                      if (current_attraction > attraction)
                      {
                        best_moves = std::vector<int>(1, move);
                        attraction = current_attraction;
                      }
                      else if (current_attraction == attraction)
                        best_moves.push_back(move);
                    });
      do_move_all(best_moves[std::rand() % best_moves.size()]);
    }

    std::vector<int> stencil_moves_all() const
    {
      unsigned int layers = std::max(1., domain_.jump_parameter_ / std::sqrt(fields_.size()));
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
      const std::array<unsigned int, nx* ny>& domain_fields = domain_.fields();
      for (unsigned int i = 0; i < fields_.size(); ++i)
      {
        unsigned int aiming = aim(fields_[i], move);
        if (domain_fields[aiming] != number_ && domain_fields[aiming] != 0)
          return -DBL_MAX;
        for (unsigned int i = 0; i < 4; ++i)
          attraction += domain_fields[aim(aiming, direct_neigh_[i])] != number_ &&
                        domain_fields[aim(aiming, direct_neigh_[i])] != 0;
      }
      return attraction;
    }

    void do_move_all(const int move)
    {
      std::array<unsigned int, nx* ny>& domain_fields = domain_.fields();
      std::vector<unsigned int> particle_merges;
      std::for_each(
        fields_.begin(), fields_.end(),
        [&](unsigned int& field)
        {
          domain_fields[field] -= number_;
          field = aim(field, move);
          domain_fields[field] += number_;

          for (unsigned int i = 0; i < 4; ++i)
          {
            if (domain_fields[aim(field, direct_neigh_[i])] != number_ &&
                domain_fields[aim(field, direct_neigh_[i])] != 0 &&
                std::find(particle_merges.begin(), particle_merges.end(),
                          domain_fields[aim(field, direct_neigh_[i])]) == particle_merges.end())
              particle_merges.push_back(domain_fields[aim(field, direct_neigh_[i])]);
          }
        });
      if (particle_merges.empty())
        return;
      std::vector<particle>& particles = domain_.particles();
      std::for_each(particle_merges.begin(), particle_merges.end(),
                    [&](unsigned int index)
                    {
                      auto other_particle = std::find(particles.begin(), particles.end(), index);
                      std::vector<unsigned int>& other_fields = other_particle->fields_;
                      for_each(other_fields.begin(), other_fields.end(),
                               [&](unsigned int ind) { domain_fields[ind] = number_; });
                      other_particle->deprecated_ = true;

                      fields_.insert(fields_.end(), other_fields.begin(), other_fields.end());
                    });
    }

    bool operator==(const unsigned int number) const { return number_ == number; }

    bool operator==(const particle& other) const { return number_ == other.number_; }

    bool operator<(const particle& other) const { return fields_.size() < other.fields_.size(); }
  };

  // -----------------------------------------------------------------------------------------------

  const double jump_parameter_;
  unsigned int n_particles;
  std::array<unsigned int, nx * ny> fields_;
  std::vector<particle> particles_;
  std::mt19937 random_seed;

 public:
  std::array<unsigned int, nx * ny>& fields() { return fields_; }
  std::vector<particle>& particles() { return particles_; }

  cellular_automaton(const double porosity, const double jump_parameter)
  : jump_parameter_(jump_parameter)
  {
    n_particles = (1. - porosity) * nx * ny;
    fields_.fill(0);
    unsigned int position = std::rand() % (nx * ny);
    for (unsigned int i = 0; i < n_particles; ++i)
    {
      while (fields_[position] != 0)
        position = std::rand() % (nx * ny);

      particles_.push_back(particle(position, *this, i + 1));
      fields_[position] = i + 1;
    }
    std::random_device rd;
    random_seed = std::mt19937(rd());
  }

  const std::array<unsigned int, nx * ny>& move_particles()
  {
    std::shuffle(particles_.begin(), particles_.end(), random_seed);
    std::for_each(particles_.begin(), particles_.end(),
                  [&](particle& part)
                  {
                    if (!part.is_deprecated())
                      part.move_all();
                  });
    particles_.erase(
      std::remove_if(particles_.begin(), particles_.end(),
                     [&](const particle& particle) -> bool { return particle.is_deprecated(); }),
      particles_.end());
    return fields_;
  }
};
