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
 * CAM illustrates the behaviour of particles in domain.
 *
 * \tparam  nx       The size of a row for each dimension of the matrix.
 *
 * \authors   Andreas Rupp, Lappeenranta-Lahti University of Technology LUT, 2022.
 * \authors   Joona Lappalainen, Lappeenranta-Lahti University of Technology LUT, 2022.
 * \authors   Simon Zech, University of Erlangen–Nuremberg, 2022.
 **************************************************************************************************/
template <auto nx>
class cellular_automaton
{
 private:
  static constexpr unsigned int dim = nx.size();
  static constexpr unsigned int n_fields()
  {
    unsigned int n_field = 1;
    for (unsigned int i = 0; i < dim; ++i)
    {
      n_field *= nx[i];
    }
    return n_field;
  }

 public:
  static constexpr unsigned int n_fields_ = n_fields();

 private:
  static constexpr std::array<int, 2 * dim> find_neigh()
  {
    static_assert(dim != 0, "Dimension of zero does not make sense.");
    std::array<int, 2 * dim> direct_neigh;
    direct_neigh[0] = -1;
    direct_neigh[1] = 1;
    for (unsigned int i = 0; i < dim - 1; ++i)
    {
      direct_neigh[2 * i + 2] = direct_neigh[2 * i] * nx[i];
      direct_neigh[2 * i + 3] = direct_neigh[2 * i + 1] * nx[i];
    }
    return direct_neigh;
  }
  /*!***********************************************************************************************
   * \brief   Array containing tentative index shifts of direct neighbors.
   ************************************************************************************************/
  static constexpr std::array<int, 2 * dim> direct_neigh_ = find_neigh();
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
    unsigned int coord;
    unsigned int new_pos = 0;
    for (unsigned int i = 0; i < dim; ++i)
    {
      coord =
        (position / direct_neigh_[2 * i + 1] + move / (int)direct_neigh_[2 * i + 1] + n_fields_) %
        nx[i];
      new_pos += coord * direct_neigh_[2 * i + 1];
    }
    return new_pos;
  }

  // -----------------------------------------------------------------------------------------------

  /*!***********************************************************************************************
   * \brief   Describes one connected set of bulk cells/fields.
   *
   * Particles can move, merge, split and deprecate.
   ************************************************************************************************/
  class particle
  {
   private:
    /*!*********************************************************************************************
     * \brief   Index of the particle.
     **********************************************************************************************/
    unsigned int number_;
    /*!*********************************************************************************************
     * \brief   If particle is deprecated.
     **********************************************************************************************/
    bool deprecated_;
    /*!*********************************************************************************************
     * \brief   Location of the particle.
     **********************************************************************************************/
    std::vector<unsigned int> fields_;
    /*!*********************************************************************************************
     * \brief   Domain of the particle.
     **********************************************************************************************/
    cellular_automaton<nx>& domain_;

    unsigned int max_min_distance(unsigned int field)
    {
      unsigned int max_distance = 0;
      std::array<unsigned int, n_fields_>& domain_fields = domain_.fields_;
      std::vector<unsigned int> visits(1, field);
      domain_fields[field] = uint_max;

      unsigned int neigh_field, visits_size = visits.size();
      for (unsigned int i = 0; i < visits_size; ++i, visits_size = visits.size())
      {
        field = visits[i];
        for (unsigned int j = 0; j < 2 * dim; ++j)
        {
          neigh_field = aim(field, direct_neigh_[j]);
          if (domain_fields[neigh_field] == number_)
          {
            visits.push_back(neigh_field);
            domain_fields[neigh_field] = domain_fields[field] - 1;
          }
        }
      }
      max_distance = uint_max - domain_fields[field];
      std::for_each(fields_.begin(), fields_.end(),
                    [&](const unsigned int field_loc) { domain_fields[field_loc] = number_; });
      return max_distance;
    }

    unsigned int directed_max_min_distance(unsigned int dir_dim)
    {
      unsigned int max_distance = 0, field = fields_[0];
      std::array<unsigned int, n_fields_>& domain_fields = domain_.fields_;
      std::vector<unsigned int> visits(1, field);
      domain_fields[field] = uint_max - n_fields_;
      unsigned int min_val = domain_fields[field], max_val = domain_fields[field];

      unsigned int neigh_field, visits_size = visits.size();
      for (unsigned int i = 0; i < visits_size; ++i, visits_size = visits.size())
      {
        field = visits[i];
        for (unsigned int j = 0; j < 2 * dim; ++j)
        {
          neigh_field = aim(field, direct_neigh_[j]);
          if (domain_fields[neigh_field] == number_)
          {
            visits.push_back(neigh_field);
            if (j / 2 == dir_dim && j % 2 == 0)
            {
              domain_fields[neigh_field] = domain_fields[field] - 1;
              if (domain_fields[neigh_field] < min_val)
                min_val = domain_fields[neigh_field];
            }
            else if (j / 2 == dir_dim)  // && j % 2 == 1
            {
              domain_fields[neigh_field] = domain_fields[field] + 1;
              if (domain_fields[neigh_field] > max_val)
                max_val = domain_fields[neigh_field];
            }
            else  // j / 2 != dir_dim
              domain_fields[neigh_field] = domain_fields[field];
          }
        }
      }
      max_distance = max_val - min_val;
      std::for_each(fields_.begin(), fields_.end(),
                    [&](const unsigned int field_loc) { domain_fields[field_loc] = number_; });
      return max_distance;
    }

   public:
    particle(const unsigned int location, cellular_automaton<nx>& domain, const unsigned int number)
    : number_(number), deprecated_(false), fields_(1, location), domain_(domain)
    {
      domain_.fields_[fields_[0]] = number_;
    }

    particle(const std::vector<unsigned int>& fields,
             cellular_automaton<nx>& domain,
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
      const std::array<unsigned int, n_fields_>& domain_fields = domain_.fields_;

      std::for_each(fields_.begin(), fields_.end(),
                    [&](const unsigned int field)
                    {
                      for (unsigned int i = 0; i < 2 * dim; ++i)
                        if (domain_fields[aim(field, direct_neigh_[i])] == 0)
                          ++n_surfaces;
                    });
      return n_surfaces;
    }

    unsigned int max_min_distance()
    {
      std::array<unsigned int, n_fields_>& domain_fields = domain_.fields_;
      std::vector<unsigned int> starts(1, fields_[0]);
      std::vector<bool> bigger(1, true);
      unsigned int start_index = 0, end_index = 1, neigh_field;
      unsigned int max_distance = max_min_distance(starts[0]);
      bool found_larger = true;

      while (found_larger)
      {
        found_larger = false;

        for (unsigned int i = start_index; i < end_index; ++i)
        {
          if (!bigger[i - start_index])
            continue;
          for (unsigned int j = 0; j < 2 * dim; ++j)
          {
            neigh_field = aim(starts[i], direct_neigh_[j]);
            if (domain_fields[neigh_field] == number_ &&
                std::find(starts.begin(), starts.end(), neigh_field) == starts.end())
              starts.push_back(neigh_field);
          }
        }
        start_index = end_index;
        end_index = starts.size();
        bigger = std::vector<bool>(end_index - start_index, true);

        for (unsigned int i = start_index; i < end_index; ++i)
          if (max_min_distance(starts[i]) > max_distance)
            found_larger = true;
          else
            bigger[i - start_index] = false;

        if (found_larger)
          ++max_distance;
      }

      return max_distance;
    }

    double max_dimension_ratio()
    {
      std::array<unsigned int, dim> width_dim;
      for (unsigned int i = 0; i < dim; ++i)
        width_dim[i] = directed_max_min_distance(i);
      unsigned int max_diameter = *std::max_element(width_dim.begin(), width_dim.end());
      unsigned int min_diameter = *std::min_element(width_dim.begin(), width_dim.end());
      return (double)max_diameter / (double)min_diameter;
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
     * \brief   Checks possible moves for merged particles.
     *
     * Order of possible moves: down, right, up, left.
     *
     * \retval stencil
     **********************************************************************************************/
    std::vector<int> stencil_moves_all() const
    {
      unsigned int layers =
        std::max(1., domain_.jump_parameter_ / std::pow(fields_.size(), 1.0 / (double)dim));
      std::vector<int> stencil(1, 0);
      unsigned int index = 0;
      unsigned int old_size = stencil.size();

      for (unsigned int lay = 0; lay < layers; ++lay)
      {
        for (; index < old_size; ++index)
          for (unsigned int i = 0; i < 2 * dim; ++i)
            if (std::find(stencil.begin(), stencil.end(), stencil[index] + direct_neigh_[i]) ==
                stencil.end())
              stencil.push_back(stencil[index] + direct_neigh_[i]);
        old_size = stencil.size();
      }
      return stencil;
    }
    /*!*********************************************************************************************
     * \brief   Checks possible moves for single particles.
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
          for (unsigned int i = 0; i < 2 * dim; ++i)
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
      const std::array<unsigned int, n_fields_>& domain_fields = domain_.fields_;
      for (unsigned int i = 0; i < fields_.size(); ++i)
      {
        unsigned int aiming = aim(fields_[i], move);
        if (domain_fields[aiming] != number_ && domain_fields[aiming] != 0)
          return double_min;
        for (unsigned int i = 0; i < 2 * dim; ++i)
          attraction += domain_fields[aim(aiming, direct_neigh_[i])] != number_ &&
                        domain_fields[aim(aiming, direct_neigh_[i])] != 0;
      }
      return attraction;
    }
    /*!*********************************************************************************************
     * \brief   Check attraction of the possible moves for single particles.
     *
     * \param   field           Location of the single particle
     * \param   move            Index shift induced by possible move.
     * \retval  attraction      Amount of neighbours.
     **********************************************************************************************/
    inline double check_move_single(const unsigned int field, const int move)
    {
      double attraction = 0.;
      std::array<unsigned int, n_fields_>& domain_fields = domain_.fields_;
      unsigned int aiming = aim(field, move);
      if (domain_fields[aiming] != 0 && move != 0)
        return double_min;
      domain_fields[field] = 0;
      for (unsigned int i = 0; i < 2 * dim; ++i)
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
      std::array<unsigned int, n_fields_>& domain_fields = domain_.fields_;

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
     * \param   field     Location of the particle.
     * \param   move      Index shift induced by possible move.
     **********************************************************************************************/
    inline void do_move_single(unsigned int& field, const int move)
    {
      std::array<unsigned int, n_fields_>& domain_fields = domain_.fields_;

      domain_fields[field] -= number_;
      field = aim(field, move);
      domain_fields[field] += number_;
    }
    /*!*********************************************************************************************
     * \brief   Merges or splits particles.
     **********************************************************************************************/
    void update_particle()
    {
      std::array<unsigned int, n_fields_>& domain_fields = domain_.fields_;
      auto first_with_number = std::find_if(fields_.begin(), fields_.end(),
                                            [&](const unsigned int field) -> bool
                                            { return domain_fields[field] == number_; });
      // Checks if there is a particle inside the domain.
      if (first_with_number == fields_.end())  // If the particle is not inside the domain.
      {
        deprecated_ = true;  // It is deprecated.
        return;
      }

      std::vector<unsigned int> new_fields(1, *first_with_number);
      unsigned int neigh_field, neigh_num, field;
      unsigned int fields_size = new_fields.size();

      for (unsigned int k = 0; k < fields_size; ++k, fields_size = new_fields.size())
      {
        field = new_fields[k];
        for (unsigned int i = 0; i < 2 * dim; ++i)
        {
          neigh_field = aim(field, direct_neigh_[i]);
          // If neigh_field is part of the same particle and it's already in new_fields.
          if (domain_fields[neigh_field] == number_ &&
              std::find(new_fields.begin(), new_fields.end(), neigh_field) == new_fields.end())
            new_fields.push_back(neigh_field);  // Puts neigh_field in the end of new_fields.
        }
      }

      if (new_fields.size() != fields_.size())
      {
        std::vector<unsigned int> lost_fields(fields_.size() - new_fields.size());
        std::sort(fields_.begin(), fields_.end());
        std::sort(new_fields.begin(), new_fields.end());
        // Checks difference between fields_ and new_fields and puts it in lost_fields.
        std::set_difference(fields_.begin(), fields_.end(), new_fields.begin(), new_fields.end(),
                            lost_fields.begin());
        // Erases particles from lost_fields if they have already merged to other existing particle.
        lost_fields.erase(std::remove_if(lost_fields.begin(), lost_fields.end(),
                                         [&](const unsigned int lost_field) -> bool
                                         { return domain_fields[lost_field] != number_; }),
                          lost_fields.end());

        domain_.add_particle(lost_fields);  // Particle splits into a new particle.
        std::swap(fields_, new_fields);
      }

      fields_size = fields_.size();

      for (unsigned int k = 0; k < fields_size; ++k, fields_size = fields_.size())
      {
        field = fields_[k];
        for (unsigned int i = 0; i < 2 * dim; ++i)
        {
          neigh_field = aim(field, direct_neigh_[i]);
          neigh_num = domain_fields[neigh_field];
          // If particle has a neighbour which is a part of a different particle.
          if (neigh_num == number_ || neigh_num == 0)
            continue;

          fields_.push_back(neigh_field);        // Puts neigh_field in the end of fields_.
          domain_fields[neigh_field] = number_;  // Particles merge.
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
     * \brief   Find out whether number_ is the same the same than other.number_.
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
   ************************************************************************************************/
  const double jump_parameter_;
  /*!***********************************************************************************************
   * \brief   Number of particles.
   ************************************************************************************************/
  unsigned int n_particles;
  /*!***********************************************************************************************
   * \brief   Array of particle locations.
   ************************************************************************************************/
  std::array<unsigned int, n_fields_> fields_;
  /*!***********************************************************************************************
   * \brief   Vector of particles.
   ************************************************************************************************/
  std::vector<particle> particles_;
  /*!***********************************************************************************************
   * \brief   Random seed.
   ************************************************************************************************/
  unsigned int rand_seed;

 public:
  /*!***********************************************************************************************
   * \brief   Returns fields
   *
   * \retval  fields_         Location of the particle.
   ************************************************************************************************/
  const std::array<unsigned int, n_fields_>& fields() const { return fields_; }
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

    n_particles = (1. - porosity) * n_fields_;
    fields_.fill(0);
    unsigned int position = std::rand() % (n_fields_);
    for (unsigned int i = 0; i < n_particles; ++i)
    {
      while (fields_[position] != 0)
        position = std::rand() % (n_fields_);

      particles_.push_back(particle(position, *this, i + 1));
    }
  }
  /*!***********************************************************************************************
   * \brief   Moves all particles
   *
   * \retval  fields_     Domain
   ************************************************************************************************/
  const std::array<unsigned int, n_fields_>& move_particles()
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
   * n_single_cells        number of single solid pixels without solid neighbours
   * n_particles           number of solid particles, including single solid pixels
   *                       and agglomorates of solid pixels
   * n_solids              total number of solid pixels
   * n_surfaces            total solid surface
   * mean_particle_size    mean particle size
   * n_connected_fluids    number of connected fluid
   *
   * \retval  array     Array of measure parameters.
   ************************************************************************************************/
  std::array<double, 12> eval_measures()
  {
    unsigned int n_single_cells =
      std::count_if(particles_.begin(), particles_.end(),
                    [](const particle& part) -> bool { return part.size() == 1; });

    unsigned int n_particles = particles_.size();

    unsigned int n_solids = 0;
    unsigned int n_surfaces = 0;
    unsigned int n_solids_part = 0;
    unsigned int n_surfaces_part = 0;
    unsigned int max_max_min_distance = 0;
    unsigned int local_max_min_distance = 0;

    double mean_sphericity = 0;
    double local_sphericity = 0;
    double max_diameters_ratio = 0;
    double local_diameters_ratio = 0;

    std::for_each(particles_.begin(), particles_.end(),
                  [&](particle& part)
                  {
                    n_solids_part = part.size();
                    n_surfaces_part = part.n_surfaces();

                    n_solids += n_solids_part;
                    n_surfaces += n_surfaces_part;

                    local_sphericity =
                      std::pow((double)n_solids_part, (double)(dim - 1) / (double)dim) /
                      (double)n_surfaces_part;
                    mean_sphericity += local_sphericity;

                    local_max_min_distance = part.max_min_distance();
                    max_max_min_distance = std::max(max_max_min_distance, local_max_min_distance);

                    local_diameters_ratio = part.max_dimension_ratio();
                    max_diameters_ratio = std::max(max_diameters_ratio, local_diameters_ratio);
                  });
    mean_sphericity /= particles_.size();

    double mean_particle_size = (double)n_solids / (double)n_particles;

    double compactness =
      std::pow(((double)n_surfaces / (double)fields_.size()), ((double)dim / (double)(dim - 1)));

    double variance_particle_sizes = 0;
    for (unsigned int i = 0; i < particles_.size(); ++i)
    {
      variance_particle_sizes += std::pow(particles_[i].size() - mean_particle_size, 2);
    }
    variance_particle_sizes /= (double)particles_.size();

    std::array<unsigned int, 2> n_fluid_components = n_fluid_comp();
    unsigned int n_connected_fluids = n_fluid_components[0];
    unsigned int n_periodic_fluid_components = n_fluid_components[1];

    return {(double)n_single_cells,
            (double)n_particles,
            (double)n_solids,
            (double)n_surfaces,
            (double)n_connected_fluids,
            (double)n_periodic_fluid_components,
            mean_particle_size,
            variance_particle_sizes,
            compactness,
            (double)max_max_min_distance,
            mean_sphericity,
            max_diameters_ratio};
  }
  /*!***********************************************************************************************
   * \brief   Computes connected fluid areas.
   *
   * \retval  n_connected_fluids     Number of connected fluid areas.
   ************************************************************************************************/
  std::array<unsigned int, 2> n_fluid_comp()
  {
    unsigned int n_connected_fluids = 0;
    unsigned int n_periodic_fluids = 0;
    unsigned int fluids_size, field, neigh_field;
    std::vector<unsigned int> found_fluids;

    bool periodic;
    std::vector<std::array<int, dim> > ref_dist(n_fields_);
    std::for_each(ref_dist.begin(), ref_dist.end(),
                  [](std::array<int, dim> dist) { dist.fill(0); });

    for (auto first_fluid = std::find(fields_.begin(), fields_.end(), 0);
         first_fluid != fields_.end(); first_fluid = std::find(first_fluid, fields_.end(), 0))
    {
      periodic = false;
      found_fluids = std::vector<unsigned int>(1, std::distance(fields_.begin(), first_fluid));
      fields_[found_fluids[0]] = uint_max;
      fluids_size = 1;
      for (unsigned int k = 0; k < fluids_size; ++k, fluids_size = found_fluids.size())
      {
        field = found_fluids[k];
        for (unsigned int i = 0; i < 2 * dim; ++i)
        {
          neigh_field = aim(field, direct_neigh_[i]);
          if (fields_[neigh_field] == 0)
          {
            fields_[neigh_field] = uint_max;
            found_fluids.push_back(neigh_field);
            ref_dist[neigh_field] = ref_dist[field];
            ref_dist[neigh_field][i / 2] += 2 * (i % 2) - 1;
          }
          else if (fields_[neigh_field] == uint_max)
            for (unsigned int j = 0; j < dim; ++j)
              if ((unsigned int)std::abs(ref_dist[field][j] - ref_dist[neigh_field][j]) > nx[j] - 2)
                periodic = true;
        }
      }
      ++n_connected_fluids;
      n_periodic_fluids += periodic;
    }

    std::for_each(fields_.begin(), fields_.end(),
                  [](unsigned int& field)
                  {
                    if (field == uint_max)
                      field = 0;
                  });

    return {n_connected_fluids, n_periodic_fluids};
  }
};
