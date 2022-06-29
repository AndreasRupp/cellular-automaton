#pragma once

#include <algorithm>
#include <array>
#include <cfloat>
#include <chrono>
#include <cmath>
#include <random>
#include <vector>

#include <iostream>

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
      domain_.fields_[fields_[0]] = number_;
    }

    particle(const std::vector<unsigned int>& fields,
             cellular_automaton<nx, ny>& domain,
             const unsigned int number)
    : number_(number), deprecated_(false), fields_(fields), domain_(domain)
    {
      std::for_each(fields_.begin(), fields_.end(),
                    [&](const unsigned int field) { domain_.fields_[field] = number_; });
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

    static inline unsigned int aim(const unsigned int position, const int move)
    {
      const unsigned int x_coord = (position % nx + (nx * ny + move) % nx) % nx;
      const unsigned int y_coord = (position / nx + (move / (int)nx) + ny) % ny;
      return y_coord * nx + x_coord;
    }

    bool is_deprecated() const { return deprecated_; }

    void move_all()
    {
      if (fields_.size() == 1)
        return;
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

    void move_singles()
    {
      std::vector<int> possible_moves = stencil_moves_single();
      std::shuffle(fields_.begin(), fields_.end(), std::default_random_engine(std::rand()));
      std::vector<int> best_moves(1, 0);
      double attraction;

      std::for_each(fields_.begin(), fields_.end(),
                    [&](unsigned int& field)
                    {
                      attraction = 0.;
                      best_moves = std::vector<int>(1, 0);
                      std::for_each(possible_moves.begin(), possible_moves.end(),
                                    [&](int move)
                                    {
                                      double current_attraction = check_move_single(field, move);
                                      if (current_attraction > attraction)
                                      {
                                        best_moves = std::vector<int>(1, move);
                                        attraction = current_attraction;
                                      }
                                      else if (current_attraction == attraction)
                                        best_moves.push_back(move);
                                    });
                      do_move_single(field, best_moves[std::rand() % best_moves.size()]);
                    });
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

    std::vector<int> stencil_moves_single() const
    {
      unsigned int layers = std::max(1., domain_.jump_parameter_);
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
      const std::array<unsigned int, nx* ny>& domain_fields = domain_.fields_;
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

    inline double check_move_single(const unsigned int field, const int move)
    {
      double attraction = 0.;
      std::array<unsigned int, nx* ny>& domain_fields = domain_.fields_;
      unsigned int aiming = aim(field, move);
      if (domain_fields[aiming] != 0 && move != 0)
        return -DBL_MAX;
      domain_fields[field] = 0;
      for (unsigned int i = 0; i < 4; ++i)
        attraction += domain_fields[aim(aiming, direct_neigh_[i])] != 0;
      domain_fields[field] = number_;
      return attraction;
    }

    void do_move_all(const int move)
    {
      std::array<unsigned int, nx* ny>& domain_fields = domain_.fields_;

      std::for_each(fields_.begin(), fields_.end(),
                    [&](unsigned int& field)
                    {
                      domain_fields[field] -= number_;
                      field = aim(field, move);
                      domain_fields[field] += number_;
                    });
    }

    inline void do_move_single(unsigned int& field, const int move)
    {
      std::array<unsigned int, nx* ny>& domain_fields = domain_.fields_;

      domain_fields[field] -= number_;
      field = aim(field, move);
      domain_fields[field] += number_;
    }

    void update_particle()
    {
      std::array<unsigned int, nx* ny>& domain_fields = domain_.fields_;
      auto first_with_number = std::find_if(fields_.begin(), fields_.end(),
                                            [&](const unsigned int field) -> bool
                                            { return domain_fields[field] == number_; });

      if (first_with_number == fields_.end())
      {
        deprecated_ = true;
        return;
      }

      std::vector<unsigned int> new_fields(1, *first_with_number);
      unsigned int neigh_field, neigh_num, field;
      unsigned int fields_size = new_fields.size();

      for (unsigned int k = 0; k < fields_size; ++k, fields_size = new_fields.size())
      {
        field = new_fields[k];
        for (unsigned int i = 0; i < 4; ++i)
        {
          neigh_field = aim(field, direct_neigh_[i]);
          if (domain_fields[neigh_field] == number_ &&
              std::find(new_fields.begin(), new_fields.end(), neigh_field) == new_fields.end())
            new_fields.push_back(neigh_field);
        }
      }

      if (new_fields.size() != fields_.size())
      {
        std::vector<unsigned int> lost_fields(fields_.size() - new_fields.size());
        std::sort(fields_.begin(), fields_.end());
        std::sort(new_fields.begin(), new_fields.end());
        std::set_difference(fields_.begin(), fields_.end(), new_fields.begin(), new_fields.end(),
                            lost_fields.begin());
        lost_fields.erase(std::remove_if(lost_fields.begin(), lost_fields.end(),
                                         [&](const unsigned int lost_field) -> bool
                                         { return domain_fields[lost_field] != number_; }),
                          lost_fields.end());

        domain_.add_particle(lost_fields);
        std::swap(fields_, new_fields);
      }

      fields_size = fields_.size();

      for (unsigned int k = 0; k < fields_size; ++k, fields_size = fields_.size())
      {
        field = fields_[k];
        for (unsigned int i = 0; i < 4; ++i)
        {
          neigh_field = aim(field, direct_neigh_[i]);
          neigh_num = domain_fields[neigh_field];
          if (neigh_num == number_ || neigh_num == 0)
            continue;

          fields_.push_back(neigh_field);
          domain_fields[neigh_field] = number_;
        }
      }
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
  unsigned int random_seed;

 public:
  const std::array<unsigned int, nx * ny>& fields() const { return fields_; }

  cellular_automaton(const double porosity,
                     const double jump_parameter,
                     const unsigned int rand_seed = 0)
  : jump_parameter_(jump_parameter)
  {
    if (rand_seed == 0)
    {
      random_seed = std::chrono::system_clock::now().time_since_epoch().count();
    }
    else
      random_seed = rand_seed;
    std::srand(random_seed);
    std::cout << "The random seed is set to be: " << random_seed << std::endl;

    n_particles = (1. - porosity) * nx * ny;
    fields_.fill(0);
    unsigned int position = std::rand() % (nx * ny);
    for (unsigned int i = 0; i < n_particles; ++i)
    {
      while (fields_[position] != 0)
        position = std::rand() % (nx * ny);

      particles_.push_back(particle(position, *this, i + 1));
    }
  }

  const std::array<unsigned int, nx * ny>& move_particles()
  {
    std::shuffle(particles_.begin(), particles_.end(), std::default_random_engine(std::rand()));
    std::for_each(particles_.begin(), particles_.end(), [&](particle& part) { part.move_all(); });
    std::shuffle(particles_.begin(), particles_.end(), std::default_random_engine(std::rand()));
    std::for_each(particles_.begin(), particles_.end(),
                  [&](particle& part) { part.move_singles(); });
    unsigned int particles_size = particles_.size();
    for (unsigned int k = 0; k < particles_size; ++k, particles_size = particles_.size())
      particles_[k].update_particle();
    particles_.erase(
      std::remove_if(particles_.begin(), particles_.end(),
                     [&](const particle& particle) -> bool { return particle.is_deprecated(); }),
      particles_.end());
    return fields_;
  }

  void add_particle(const std::vector<unsigned int>& fields)
  {
    particles_.push_back(particle(fields, *this, ++n_particles));
  }
};
