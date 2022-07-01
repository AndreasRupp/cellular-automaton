#pragma once

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <limits>
#include <random>
#include <vector>

/*!*************************************************************************************************
 * \brief   This class implements a cellular automaton.
 *
 * TODO: Detailed description goes here.
 *
 * \tparam  nx       Number of rows of the matrix.
 * \tparam  ny       Number of columns of the matrix. Defaults to create square matrix.
 *
 * \authors   Andreas Rupp, Lappeenranta-Lahti University of Technology LUT, 2022.
 * \authors   Joona Lappalainen, Lappeenranta-Lahti University of Technology LUT, 2022.
 * \authors   Simon Zech, University of Erlangen–Nuremberg, 2022.
 **************************************************************************************************/
template <unsigned int nx, unsigned int ny = nx>
class cellular_automaton
{
 private:
  /*!***********************************************************************************************
   * \brief   Array containing tentative index shifts of direct neighbors.
   ************************************************************************************************/
  static constexpr std::array<int, 4> direct_neigh_ = {-1 * (int)nx, 1, nx, -1};
  /*!***********************************************************************************************
   * \brief   Maximum unsigned integer.
   ************************************************************************************************/
  static constexpr unsigned int uint_max = std::numeric_limits<unsigned int>::max();
  /*!***********************************************************************************************
   * \brief   Smallest (negative) double.
   ************************************************************************************************/
  static constexpr double double_min = std::numeric_limits<double>::lowest();

  /*!***********************************************************************************************
   * \brief   Find field if one moves from position to move.
   *
   * \param   position  Current position of field that may move.
   * \param   move      Index shift induced by possible move.
   * \retval  index     Index of move target.
   ************************************************************************************************/
  static inline unsigned int aim(const unsigned int position, const int move)
  {
    const unsigned int x_coord = (position % nx + (nx * ny + move) % nx) % nx;
    const unsigned int y_coord = (position / nx + (move / (int)nx) + ny) % ny;
    return y_coord * nx + x_coord;
  }

  // -----------------------------------------------------------------------------------------------

  /*!***********************************************************************************************
   * \brief   Describes one connected set of bulk cells/fields.
   *
   * More info
   ************************************************************************************************/
  class particle
  {
   private:
    /*!*********************************************************************************************
     * \brief   Describes one connected set of bulk cells/fields.
     *
     * More info
     **********************************************************************************************/
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
      return *this;
    }
    /*!*********************************************************************************************
     * \brief   Checks if particle is deprecated.
     *
     * \retval  deprecated      True of false.
     **********************************************************************************************/
    bool is_deprecated() const { return deprecated_; }
    /*!*********************************************************************************************
     * \brief   Return size of the particle.
     *
     * \retval  size            Size of the particle.
     **********************************************************************************************/
    unsigned int size() const { return fields_.size(); }
    /*!*********************************************************************************************
     * \brief   Counts surfaces of all particles.
     *
     * \retval  n_surfaces      Surfaces of all particles.
     **********************************************************************************************/
    unsigned int n_surfaces() const
    {
      unsigned int n_surfaces = 0;
      const std::array<unsigned int, nx* ny>& domain_fields = domain_.fields_;

      std::for_each(fields_.begin(), fields_.end(),
                    [&](const unsigned int field)
                    {
                      for (unsigned int i = 0; i < 4; ++i)
                        if (domain_fields[aim(field, direct_neigh_[i])] == 0)
                          ++n_surfaces;
                    });
      return n_surfaces;
    }
    /*!*********************************************************************************************
     * \brief   Moves merged particles.
     **********************************************************************************************/
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
    /*!*********************************************************************************************
     * \brief   Moves single particles.
     **********************************************************************************************/
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
    /*!*********************************************************************************************
     * \brief   Checks possible moves.
     *
     * Order of possible moves: down, right, up, left.
     *
     * \retval stencil
     **********************************************************************************************/
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
    /*!*********************************************************************************************
     * \brief   Checks possible moves for a single particle.
     *
     * Order of possible moves: down, right, up, left.
     *
     * \retval stencil
     **********************************************************************************************/
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
    /*!*********************************************************************************************
     * \brief   Check attraction of the possible moves for merged particles.
     *
     * \param   move            Index shift induced by possible move.
     * \retval  attraction      Amount of neighbours.
     **********************************************************************************************/
    double check_move_all(const int move) const
    {
      double attraction = 0.;
      const std::array<unsigned int, nx* ny>& domain_fields = domain_.fields_;
      for (unsigned int i = 0; i < fields_.size(); ++i)
      {
        unsigned int aiming = aim(fields_[i], move);
        if (domain_fields[aiming] != number_ && domain_fields[aiming] != 0)
          return double_min;
        for (unsigned int i = 0; i < 4; ++i)
          attraction += domain_fields[aim(aiming, direct_neigh_[i])] != number_ &&
                        domain_fields[aim(aiming, direct_neigh_[i])] != 0;
      }
      return attraction;
    }
    /*!*********************************************************************************************
     * \brief   Check attraction of the possible moves for a single particle.
     *
     * \param   field           Particle
     * \param   move            Index shift induced by possible move.
     * \retval  attraction      Amount of neighbours.
     **********************************************************************************************/
    inline double check_move_single(const unsigned int field, const int move)
    {
      double attraction = 0.;
      std::array<unsigned int, nx* ny>& domain_fields = domain_.fields_;
      unsigned int aiming = aim(field, move);
      if (domain_fields[aiming] != 0 && move != 0)
        return double_min;
      domain_fields[field] = 0;
      for (unsigned int i = 0; i < 4; ++i)
        attraction += domain_fields[aim(aiming, direct_neigh_[i])] != 0;
      domain_fields[field] = number_;
      return attraction;
    }
    /*!*********************************************************************************************
     * \brief   Moves merged particles.
     *
     * \param   move      Index shift induced by possible move.
     **********************************************************************************************/
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
    /*!*********************************************************************************************
     * \brief   Moves single particles.
     *
     * \param   field     Particle
     * \param   move      Index shift induced by possible move.
     **********************************************************************************************/
    inline void do_move_single(unsigned int& field, const int move)
    {
      std::array<unsigned int, nx* ny>& domain_fields = domain_.fields_;

      domain_fields[field] -= number_;
      field = aim(field, move);
      domain_fields[field] += number_;
    }
    /*!*********************************************************************************************
     * \brief   Merges or deprecates particles.
     *
     * More info
     **********************************************************************************************/
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
    /*!*********************************************************************************************
     * \brief   Find out whether two numbers are the same.
     *
     * \param   number     Number to compare.
     * \retval  isEqual    True of false.
     **********************************************************************************************/
    bool operator==(const unsigned int number) const { return number_ == number; }
    /*!*********************************************************************************************
     * \brief   Find out whether number_ is the same the same than other-number_.
     *
     * \param   other      Number to compare.
     * \retval  isEqual    True of false.
     **********************************************************************************************/
    bool operator==(const particle& other) const { return number_ == other.number_; }
    /*!*********************************************************************************************
     * \brief   Find out whether fields_.size() is smaller than other.fields_.size().
     *
     * \param   other      Size to compare.
     * \retval  smaller    True of false.
     **********************************************************************************************/
    bool operator<(const particle& other) const { return fields_.size() < other.fields_.size(); }
  };

  // -----------------------------------------------------------------------------------------------

  /*!***********************************************************************************************
   * \brief   Jump parameter.
   *
   * More info
   ************************************************************************************************/
  const double jump_parameter_;
  /*!***********************************************************************************************
   * \brief   Particle class.
   *
   * More info
   ************************************************************************************************/
  unsigned int n_particles;
  /*!***********************************************************************************************
   * \brief   Particle class.
   *
   * More info
   ************************************************************************************************/
  std::array<unsigned int, nx * ny> fields_;
  /*!***********************************************************************************************
   * \brief   Particle class.
   *
   * More info
   ************************************************************************************************/
  std::vector<particle> particles_;
  /*!***********************************************************************************************
   * \brief   Particle class.
   *
   * More info
   ************************************************************************************************/
  unsigned int rand_seed;

 public:
  /*!***********************************************************************************************
   * \brief   Returns fields
   *
   * \retval  fields_
   ************************************************************************************************/
  const std::array<unsigned int, nx * ny>& fields() const { return fields_; }
  /*!***********************************************************************************************
   * \brief   Cellular automaton.
   *
   * \param   porosity          The percentage of void space, not occupied by solid.
   * \param   jump_parameter    How far individual particles are allowed to jump.
   * \param   random_seed       If given, sets random seed to given seed.
   ************************************************************************************************/
  cellular_automaton(const double porosity,
                     const double jump_parameter,
                     const unsigned int random_seed = 0)
  : jump_parameter_(jump_parameter)
  {
    if (random_seed == 0)
    {
      rand_seed = std::chrono::system_clock::now().time_since_epoch().count();
    }
    else
      rand_seed = random_seed;
    std::srand(rand_seed);

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
  /*!***********************************************************************************************
   * \brief   Moves all particles
   *
   * \retval  fields_     Domain
   ************************************************************************************************/
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
  /*!***********************************************************************************************
   * \brief   Creates particles.
   *
   * \param   fields
   ************************************************************************************************/
  void add_particle(const std::vector<unsigned int>& fields)
  {
    particles_.push_back(particle(fields, *this, ++n_particles));
  }
  /*!***********************************************************************************************
   * \brief   Sets a random seed.
   *
   * \retval  rand_seed     Random seed.
   ************************************************************************************************/
  unsigned int random_seed() const { return rand_seed; }
  /*!***********************************************************************************************
   * \brief   Evaluates measure parameters.
   *
   * Measure parameters:
   * n_single_cells        (number of single solid pixels without solid neighbours)
   * n_particles           (number of solid particles, including single solid pixels
   *                        and agglomorates of solid pixels)
   * n_solids              (total number of solid pixels)
   * n_surfaces            (total solid surface)
   * mean_particle_size    (mean particle size)
   * n_connected_fluids    (number of connected fluid)
   *
   * \retval  array     Array of measure parameters.
   ************************************************************************************************/
  std::array<double, 6> eval_measures()
  {
    unsigned int n_single_cells =
      std::count_if(particles_.begin(), particles_.end(),
                    [](const particle& part) -> bool { return part.size() == 1; });

    unsigned int n_particles = particles_.size();

    unsigned int n_solids = 0;
    unsigned int n_surfaces = 0;

    std::for_each(particles_.begin(), particles_.end(),
                  [&](const particle& part)
                  {
                    n_solids += part.size();
                    n_surfaces += part.n_surfaces();
                  });

    double mean_particle_size = (double)n_solids / (double)n_particles;

    unsigned int n_connected_fluids = n_fluid_comp();

    return {(double)n_single_cells, (double)n_particles, (double)n_solids,
            (double)n_surfaces,     mean_particle_size,  (double)n_connected_fluids};
  }
  /*!***********************************************************************************************
   * \brief   Computes connected fluid areas.
   *
   * \retval  n_connected_fluids     Number of connected fluid.
   ************************************************************************************************/
  unsigned int n_fluid_comp()
  {
    unsigned int n_connected_fluids = 0;
    unsigned int fluids_size, field, neigh_field;
    std::vector<unsigned int> found_fluids;

    for (auto first_fluid = std::find(fields_.begin(), fields_.end(), 0);
         first_fluid != fields_.end(); first_fluid = std::find(first_fluid, fields_.end(), 0))
    {
      found_fluids = std::vector<unsigned int>(1, std::distance(fields_.begin(), first_fluid));
      fields_[found_fluids[0]] = uint_max;
      fluids_size = 1;
      for (unsigned int k = 0; k < fluids_size; ++k, fluids_size = found_fluids.size())
      {
        field = found_fluids[k];
        for (unsigned int i = 0; i < 4; ++i)
        {
          neigh_field = aim(field, direct_neigh_[i]);
          if (fields_[neigh_field] == 0)
          {
            fields_[neigh_field] = uint_max;
            found_fluids.push_back(neigh_field);
          }
        }
      }
      ++n_connected_fluids;
    }

    std::for_each(fields_.begin(), fields_.end(),
                  [](unsigned int& field)
                  {
                    if (field == uint_max)
                      field = 0;
                  });

    return n_connected_fluids;
  }
};
