/*!*********************************************************************************************
 * \file evaluation.hxx
 * \brief Implementation of different norms and metrics for evaluation the domain
 **********************************************************************************************/
#pragma once
#include <CAM/domain.hxx>
#include <CAM/utils.hxx>

namespace CAM
{
template <auto nx, typename fields_array_t>
class Evaluation
{
 private:
  CAM::Domain<nx, fields_array_t>* domain;

 public:
  Evaluation() {}
  Evaluation(Domain<nx, fields_array_t>* _domain) { domain = _domain; }

  /*!*********************************************************************************************
   * \brief Compares bulks of domain A and B
   *
   * \param domain_a
   * \param domain_b
   * \return amount of cells which alternated from or to 0
   ************************************************************************************************/
  static constexpr unsigned int bulk_distance(const fields_array_t& domain_a,
                                              const fields_array_t& domain_b)
  {
    unsigned int distance = 0;
    for (unsigned int field = 0; field < domain_a.size(); ++field)
      distance += (domain_a[field] == 0) != (domain_b[field] == 0);
    return distance;
  }

  static constexpr unsigned int skeleton_distance(const fields_array_t& domain_a,
                                                  const fields_array_t& domain_b)
  {
    unsigned int neigh_field, distance = 0;
    for (unsigned int field = 0; field < domain_a.size(); ++field)
      for (unsigned int j = 0; j < 2 * nx.size(); j += 2)
      {
        neigh_field = aim<nx>(field, direct_neigh<nx>(j));
        distance += std::abs((domain_a[field] == 0) - (domain_a[neigh_field] == 0) +
                             (domain_b[neigh_field] == 0) - (domain_b[field] == 0));
      }
    return distance;
  }
  /*!*********************************************************************************************
   * \brief Gives size of each particle in a sorted way
   *
   * \return vector of sizes
   ************************************************************************************************/
  std::vector<unsigned int> particle_size_distribution()
  {
    // \deprecated

    // unsigned int fluids_size, field, neigh_field;
    // std::vector<unsigned int> found_solids, distribution1;
    // fields_array_t domain_fields =domain->domain_fields;
    // std::for_each(domain_fields.begin(), domain_fields.end(),
    //               [](unsigned int& field) { field = (field == 0); });

    // for (auto first_fluid = std::find(domain_fields.begin(), domain_fields.end(), 0);
    //      first_fluid != domain_fields.end();
    //      first_fluid = std::find(first_fluid, domain_fields.end(), 0))
    // {
    //   found_solids = std::vector<unsigned int>(1, std::distance(domain_fields.begin(),
    //   first_fluid)); domain_fields[found_solids[0]] = CAM::uint_max; fluids_size = 1; for
    //   (unsigned int k = 0; k < fluids_size; ++k, fluids_size = found_solids.size())
    //   {
    //     field = found_solids[k];
    //     for (unsigned int i = 0; i < 2 * nx.size(); ++i)
    //     {
    //       neigh_field = aim<nx>(field, direct_neigh<nx>(i));
    //       if (domain_fields[neigh_field] == 0)
    //       {
    //         domain_fields[neigh_field] = uint_max;
    //         found_solids.push_back(neigh_field);
    //       }
    //     }
    //   }
    //   distribution1.push_back(found_solids.size());
    // }
    // std::sort(distribution.begin(), distribution.end());

    domain->find_composites_via_bu_border();
    std::vector<unsigned int> distribution;
    distribution.resize(domain->particles.size());
    for (unsigned int i = 0; i < domain->particles.size(); i++)
    {
      distribution.at(i) = domain->particles.at(i).field_indices.size();
    }
    std::sort(distribution.begin(), distribution.end());

    return distribution;
  }
  /*!*********************************************************************************************
   * \brief Gives the number of particles (connected component) inside the domain
   *
   * \return amount of particles
   ***************************************************************************************************/
  unsigned int n_solid_comp() { return particle_size_distribution().size(); }
  /*!*********************************************************************************************
   * \brief   Counts surfaces of particle.
   * \param _fields fields of particle
   * \retval  n_surfaces      Surfaces of particle.
   **********************************************************************************************/
  unsigned int getNSurfaces(const std::vector<unsigned int>& _fields)
  {
    unsigned int n_surfaces = 0;
    std::for_each(_fields.begin(), _fields.end(),
                  [&](const unsigned int field)
                  {
                    for (unsigned int i = 0; i < 2 * nx.size(); ++i)
                      if (domain->domain_fields[aim<nx>(field, direct_neigh<nx>(i))] == 0)
                        ++n_surfaces;
                  });
    return n_surfaces;
  }
  /*!*********************************************************************************************
   * \brief   Find the longest shortest path of particle inside domain->
   *
   * \retval  max_distance     Length of the path.
   **********************************************************************************************/
  unsigned int max_min_distance(const Particle& _particle)
  {
    std::vector<unsigned int> starts(1, _particle.field_indices[0]);
    std::vector<bool> bigger(1, true);
    unsigned int start_index = 0, end_index = 1, neigh_field;
    unsigned int max_distance = max_min_distance(starts[0], _particle);
    bool found_larger = true;

    while (found_larger)
    {
      found_larger = false;

      for (unsigned int i = start_index; i < end_index; ++i)
      {
        if (!bigger[i - start_index])
          continue;
        for (unsigned int j = 0; j < 2 * nx.size(); ++j)
        {
          neigh_field = aim<nx>(starts[i], direct_neigh<nx>(j));
          if (std::find(_particle.numbers.begin(), _particle.numbers.end(),
                        domain->domain_fields[neigh_field]) != _particle.numbers.end() &&
              std::find(starts.begin(), starts.end(), neigh_field) == starts.end())
            starts.push_back(neigh_field);
        }
      }
      start_index = end_index;
      end_index = starts.size();
      bigger = std::vector<bool>(end_index - start_index, true);

      for (unsigned int i = start_index; i < end_index; ++i)
        if (max_min_distance(starts[i], _particle) > max_distance)
          found_larger = true;
        else
          bigger[i - start_index] = false;

      if (found_larger)
        ++max_distance;
    }

    return max_distance;
  }
  /*!*********************************************************************************************
   * \brief   Calculate the ratio of diameters for every dimension.
   *
   * \retval  max_diameter / min_diameter      Max dimension ratio.
   **********************************************************************************************/
  double max_dimension_ratio(const Particle& _particle)
  {
    std::array<unsigned int, nx.size()> width_dim;
    for (unsigned int i = 0; i < nx.size(); ++i)
      width_dim[i] = directed_max_min_distance(i, _particle);
    unsigned int max_diameter = *std::max_element(width_dim.begin(), width_dim.end());
    unsigned int min_diameter = *std::min_element(width_dim.begin(), width_dim.end());
    return (double)max_diameter / (double)min_diameter;
  }

  /*!*********************************************************************************************
   * \brief   Find the longest shortest path inside particle.
   *
   * \param   field            Location of the single particle
   * \retval  max_distance     Length of the path
   **********************************************************************************************/
  unsigned int max_min_distance(unsigned int field, const Particle& _particle)
  {
    unsigned int max_distance = 0;
    fields_array_t domain_fields = domain->domain_fields;
    std::vector<unsigned int> visits(1, field);
    domain_fields[field] = uint_max;

    unsigned int neigh_field, visits_size = visits.size();
    for (unsigned int i = 0; i < visits_size; ++i, visits_size = visits.size())
    {
      field = visits[i];
      for (unsigned int j = 0; j < 2 * nx.size(); ++j)
      {
        neigh_field = aim<nx>(field, direct_neigh<nx>(j));
        if (std::find(_particle.numbers.begin(), _particle.numbers.end(),
                      domain_fields[neigh_field]) != _particle.numbers.end())
        {
          visits.push_back(neigh_field);
          domain_fields[neigh_field] = domain_fields[field] - 1;
        }
      }
    }
    max_distance = uint_max - domain_fields[field];
    return max_distance;
  }
  /*!*********************************************************************************************
   * \brief   Find the longest amount of moves in certain direction inside particle.
   *
   * \param   dir_dim          Certain dimension
   * \retval  max_distance     Amount of moves
   **********************************************************************************************/
  unsigned int directed_max_min_distance(unsigned int dir_dim, const Particle& _particle)
  {
    unsigned int max_distance = 0, field = _particle.field_indices[0];
    fields_array_t domain_fields = domain->domain_fields;
    std::vector<unsigned int> visits(1, field);
    domain_fields[field] = uint_max - n_fields<nx>();
    unsigned int min_val = domain_fields[field], max_val = domain_fields[field];

    unsigned int neigh_field, visits_size = visits.size();
    for (unsigned int i = 0; i < visits_size; ++i, visits_size = visits.size())
    {
      field = visits[i];
      for (unsigned int j = 0; j < 2 * nx.size(); ++j)
      {
        neigh_field = aim<nx>(field, direct_neigh<nx>(j));
        if (std::find(_particle.numbers.begin(), _particle.numbers.end(),
                      domain_fields[neigh_field]) != _particle.numbers.end())
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
    max_distance = max_val - min_val + 1;
    return max_distance;
  }
  /*!*********************************************************************************************
   * \brief Gives the avarage particle size
   *
   * \return avarage particle size
   ************************************************************************************************/
  constexpr double average_particle_size()
  {
    unsigned int n_solids = 0;
    for (unsigned int i = 0; i < domain->domain_fields.size(); ++i)
      n_solids += (domain->domain_fields[i] != 0);
    return (double)n_solids / (double)n_solid_comp();
  }
  /*!***********************************************************************************************
   * \brief   Evaluates measure parameters.
   *
   * Measure parameters:
   * n_single_cells                 number of cells that are not part of any larger agglomerate
   * n_particles                    number of connected solid cells, including single cells and
   *                                larger agglomerates
   * n_solids                       total number of solid cells
   * n_surfaces                     total number of faces between solid and fluid
   * n_connected_fluids             number of connected fluids
   * n_periodic_fluid_components    number of connected fluids which are periodic
   * mean_particle_size             average particle size (n_solids / n_particles)
   * variance_particle_sizes        variance of particle sizes
   * compactness                    evaluates the degree of which a particle is compact
   * max_max_min_distance           maximum over all particles with respect to the maximum of
   *                                shortest distances between two solid cells within a particle
   * mean_sphericity                sphericity rates how close a shape is to the perfect sphere
   * max_diameters_ratio            maximum ratio of diameters with respect to one dimension
   *
   * \retval  array                 Array of measure parameters.
   ************************************************************************************************/
  std::array<double, 12> eval_measures()
  {
    domain->find_composites_via_bu_border();
    unsigned int n_single_cells =
      std::count_if(domain->particles.begin(), domain->particles.end(),
                    [](Particle particle) -> bool { return particle.field_indices.size() == 1; });

    unsigned int n_particles = domain->particles.size();

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

    std::for_each(domain->particles.begin(), domain->particles.end(),
                  [&](Particle particle)
                  {
                    n_solids_part = particle.field_indices.size();
                    n_surfaces_part = getNSurfaces(particle.field_indices);

                    n_solids += n_solids_part;
                    n_surfaces += n_surfaces_part;

                    local_sphericity =
                      std::pow((double)n_solids_part, (double)(nx.size() - 1) / (double)nx.size()) /
                      (double)n_surfaces_part;
                    mean_sphericity += local_sphericity;

                    local_max_min_distance = max_min_distance(particle);
                    max_max_min_distance = std::max(max_max_min_distance, local_max_min_distance);

                    local_diameters_ratio = max_dimension_ratio(particle);
                    max_diameters_ratio = std::max(max_diameters_ratio, local_diameters_ratio);
                  });
    mean_sphericity /= n_particles;

    double mean_particle_size = (double)n_solids / (double)n_particles;

    double compactness =
      std::pow((double)n_surfaces, ((double)nx.size() / (double)(nx.size() - 1))) /
      (double)n_solids;

    double variance_particle_sizes = 0;
    for (unsigned int i = 0; i < n_particles; ++i)
    {
      variance_particle_sizes +=
        std::pow(domain->particles[i].field_indices.size() - mean_particle_size, 2);
    }
    variance_particle_sizes /= (double)n_particles;

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
   * \retval  n_fluid_comp      1st entry of array is connected fluid areas.
   *                            2nd entry of array is periodic connected fluid areas.
   ************************************************************************************************/
  std::array<unsigned int, 2> n_fluid_comp()
  {
    unsigned int n_connected_fluids = 0;
    unsigned int n_periodic_fluids = 0;
    unsigned int fluids_size, field, neigh_field;
    std::vector<unsigned int> found_fluids;
    fields_array_t domain_fields = domain->domain_fields;
    bool periodic;
    std::vector<std::array<int, nx.size()>> ref_dist(n_fields<nx>());
    std::for_each(ref_dist.begin(), ref_dist.end(),
                  [](std::array<int, nx.size()> dist) { dist.fill(0); });

    for (auto first_fluid = std::find(domain_fields.begin(), domain_fields.end(), 0);
         first_fluid != domain_fields.end();
         first_fluid = std::find(first_fluid, domain_fields.end(), 0))
    {
      periodic = false;
      found_fluids =
        std::vector<unsigned int>(1, std::distance(domain_fields.begin(), first_fluid));
      domain_fields[found_fluids[0]] = uint_max;
      fluids_size = 1;
      for (unsigned int k = 0; k < fluids_size; ++k, fluids_size = found_fluids.size())
      {
        field = found_fluids[k];
        for (unsigned int i = 0; i < 2 * nx.size(); ++i)
        {
          neigh_field = aim<nx>(field, direct_neigh<nx>(i));
          if (domain_fields[neigh_field] == 0)
          {
            domain_fields[neigh_field] = uint_max;
            found_fluids.push_back(neigh_field);
            ref_dist[neigh_field] = ref_dist[field];
            ref_dist[neigh_field][i / 2] += 2 * (i % 2) - 1;
          }
          else if (domain_fields[neigh_field] == uint_max)
            for (unsigned int j = 0; j < nx.size(); ++j)
              if ((unsigned int)std::abs(ref_dist[field][j] - ref_dist[neigh_field][j]) > nx[j] - 2)
                periodic = true;
        }
      }
      ++n_connected_fluids;
      n_periodic_fluids += periodic;
    }

    return {n_connected_fluids, n_periodic_fluids};
  }
};
}  // namespace CAM