/**
 * @file domain.hxx
 *
 * @brief This class implements a domain and holds all its properties
 * TODO Mayber outsource all evaluation functions to own class
 *
 * \tparam  nx              The size of a row for each dimension of the matrix
 * \tparam  n_fields        Size of the domain
 *
 */
#pragma once
#ifndef DOMAIN_HXX
#define DOMAIN_HXX

#include <CAM/aggregate.hxx>
#include <CAM/building_units.hxx>
#include <CAM/utils.hxx>
#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
namespace CAM
{
template <auto nx, typename fields_array_t>
class Domain
{
  static const unsigned int dim = nx.size();

 public:
  double jump_parameter_composites;
  static constexpr unsigned int n_fields_ = n_fields<nx>();
  // storing index for new particle
  double indexBU;

  ~Domain()
  {
    std::for_each(buildingUnits.begin(), buildingUnits.end(),
                  [&](CAM::BuildingUnit<nx>* unit) { delete unit; });
    std::for_each(aggregates.begin(), aggregates.end(),
                  [&](CAM::Aggregate<nx>* aggregate) { delete aggregate; });
  }
  /**
   * @brief Construct a new Domain object
   *
   * @param _jump_parameter_composites How far composites particles are allowed to jump.
   */
  Domain(double _jump_parameter_composites = -1)
  : jump_parameter_composites((_jump_parameter_composites != -1.)
                                ? _jump_parameter_composites
                                : 5)  // double _jump_parameter_composites = 1.0
  {
    if constexpr (std::is_same<fields_array_t,
                               std::vector<typename fields_array_t::value_type>>::value)
      domainFields.resize(n_fields_, 0);
    else
    {
      static_assert(
        std::is_same<fields_array_t,
                     std::array<typename fields_array_t::value_type, n_fields<nx>()>>::value,
        "The fields array has incorrect size");
      domainFields.fill(0);
      indexBU = 0;
    }
  }
  /**
   * @brief Place any building unit (BU) into the domain
   * Basis function for placement for all special BUs
   * @param _unit building unit
   * @return true: cells pf BU are marked with individual index, BU is stored in BU vector
   * @return false: BU could not be placed. all its cells are removed
   */
  bool placeBU(CAM::BuildingUnit<nx>* _unit)
  {
    bool success = true;
    indexBU++;
    _unit->number = indexBU;
    // std::cout<<_unit->number<<" ";
    std::vector<unsigned int> fields = _unit->getFieldIndices();
    std::for_each(fields.begin(), fields.end(),
                  [&](unsigned int field)
                  {
                    if (domainFields[field] == 0)
                    {
                      success = success && true;
                      domainFields[field] = _unit->number;
                    }
                    else
                    {
                      success = 0;
                    }
                  });

    if (success == 0)
    {
      std::for_each(fields.begin(), fields.end(),
                    [&](unsigned int field)
                    {
                      if (domainFields[field] == _unit->number)
                        domainFields[field] = 0;
                    });
      indexBU--;
    }
    else
    {
      buildingUnits.push_back(_unit);
    }
    return success;
  }
  /**
   * @brief Places HyperSphere (2D: Circle, 3D: Sphere) into domain
   *
   * @param _position field index of center point; no position is given (-1) -> random position
   * @param _radius all cells are completely within the radius
   * @param _jump_parameter How far sphere is allowed to jump.
   * @return true: sphere is placed
   * @return false: sphere could not be placed. all its cells are removed
   */
  bool placeSphere(int _position = -1, double _radius = 1, double _jump_parameter = 1)
  {
    if (_position == -1)
    {
      rand_seed = std::chrono::system_clock::now().time_since_epoch().count();
      std::srand(rand_seed);
      _position = std::rand() % (n_fields_);
    }
    HyperSphereBU<nx>* sphere = new HyperSphereBU<nx>(1, _jump_parameter, _position, _radius);
    bool success = placeBU(sphere);
    if (!success)
      delete sphere;
    return success;
  }
  /**
   * @brief Place a limited hyper plane (2D: Rectangle, 3D: cuboid) into domain
   *
   * @param _position field index of center point; no position is given (-1) -> random position
   * @param _extent size of plane in each dimension
   * @param _jump_parameter How far hyper plane is allowed to jump.
   * @return true: plane is placed
   * @return false: plane could not be placed. all its cells are removed
   */
  bool placePlane(int _position = -1,
                  std::vector<unsigned int> _extent = std::vector<unsigned int>(nx.size(), 0),
                  double _jump_parameter = 1)
  {
    if (_position == -1)
    {
      rand_seed = std::chrono::system_clock::now().time_since_epoch().count();
      std::srand(rand_seed);
      _position = std::rand() % (n_fields_);
    }
    HyperPlaneBU<nx>* plane = new HyperPlaneBU<nx>(1, _jump_parameter, _position, _extent);
    bool success = placeBU(plane);
    if (!success)
      delete plane;
    return success;
  }
  /**
   * @brief Creates domain with building units containing only one cell
   *
   * @param _porosity The percentage of void space, not occupied by solid/BU.
   * @param _jump_parameter How far individual particles are allowed to jump.
   * @param random_seed If given, sets random seed to given seed.
   */
  void placeSingleCellBURandomly(double _porosity = 0.90,
                                 double _jump_parameter = 1,
                                 unsigned int random_seed = 0)
  {
    if (random_seed == 0)
    {
      rand_seed = std::chrono::system_clock::now().time_since_epoch().count();
    }
    else
      rand_seed = random_seed;

    std::srand(rand_seed);

    unsigned int n_particles = (1. - _porosity) * n_fields_;
    unsigned int position = std::rand() % (n_fields_);
    for (unsigned int i = 0; i < n_particles; ++i)
    {
      while (domainFields[position] != 0)
        position = std::rand() % (n_fields_);
      std::vector<unsigned int> pos(1, position);
      buildingUnits.push_back(new CustomBU<nx>(i + 1, _jump_parameter, pos));
      domainFields[position] = i + 1;
    }
  }
  /*!*
   * @brief  Finds composites (particles containing more then one BU) in domain and stores
   * information in std::vector<CAM::Aggregate<nx>*> aggregates; std::vector<Particle> particles;
   */
  void findAggregates()
  {
    fields_array_t fields = domainFields;
    constexpr unsigned int dim = nx.size();
    unsigned int solids_size, field, neigh_field;
    std::vector<unsigned int> found_solids;
    particles.clear();
    std::for_each(fields.begin(), fields.end(),
                  [](CAM::fieldNumbers_t& field) { field = (field == 0); });

    for (auto first_solid = std::find(fields.begin(), fields.end(), 0); first_solid != fields.end();
         first_solid = std::find(first_solid, fields.end(), 0))
    {
      std::vector<CAM::fieldNumbers_t> aggregateComponents;

      found_solids = std::vector<unsigned int>(1, std::distance(fields.begin(), first_solid));
      fields[found_solids[0]] = uint_max;

      aggregateComponents.push_back(domainFields[found_solids[0]]);
      solids_size = 1;
      for (unsigned int k = 0; k < solids_size; ++k, solids_size = found_solids.size())
      {
        field = found_solids[k];
        for (unsigned int i = 0; i < 2 * dim; ++i)
        {
          neigh_field = aim<nx>(field, direct_neigh<nx>(i));
          if (fields[neigh_field] == 0)
          {
            fields[neigh_field] = uint_max;
            found_solids.push_back(neigh_field);
            unsigned int number = domainFields[neigh_field];
            if (std::find(aggregateComponents.begin(), aggregateComponents.end(), number) ==
                aggregateComponents.end())
            {
              aggregateComponents.push_back(number);
            }
          }
        }
      }
      if (aggregateComponents.size() > 1)
      {
        CAM::Aggregate<nx>* newAggregate = new CAM::Aggregate<nx>();
        for (unsigned int i = 0; i < aggregateComponents.size(); i++)
        {
          typename std::vector<CAM::BuildingUnit<nx>*>::iterator it =
            std::find_if(buildingUnits.begin(), buildingUnits.end(),
                         [&](CAM::BuildingUnit<nx>* unit) -> bool
                         { return unit->number == aggregateComponents[i]; });
          newAggregate->buildingUnits.push_back(*it);
          newAggregate->fieldIndices = found_solids;
          newAggregate->jump_parameter =
            std::pow(newAggregate->fieldIndices.size(), 1.0 / (double)dim);
        }
        aggregates.push_back(newAggregate);
      }
      particles.push_back(Particle(found_solids, aggregateComponents));
    }
  }
  /*!***********************************************************************************************
   * \brief   Returns Domain
   *
   * \retval  domainFields      Information about  status of each cell in domain
   ************************************************************************************************/
  const fields_array_t& fields() const { return domainFields; }
  /*!***********************************************************************************************
   * \brief   Array of particle locations.
   ************************************************************************************************/
  fields_array_t domainFields;
  /*!***********************************************************************************************
   * \brief   Vector of particles.
   ************************************************************************************************/
  std::vector<CAM::BuildingUnit<nx>*> buildingUnits;
  std::vector<CAM::Aggregate<nx>*> aggregates;

  /*!***********************************************************************************************
   * \brief   Describes one connected set of bulk cells/fields.
   * Information gathered only by analysing the domain independ of defintion of BU and composites
   ************************************************************************************************/
  struct Particle
  {
    Particle(std::vector<unsigned int> _fieldIndices, std::vector<CAM::fieldNumbers_t> _numbers)
    {
      fieldIndices = _fieldIndices;
      numbers = _numbers;
    }
    /*!*********************************************************************************************
     * \brief   Indices of building units contained in a particle
     * min_size = 1, max_size = number of building units in domain
     **********************************************************************************************/
    std::vector<CAM::fieldNumbers_t> numbers;
    /*!*********************************************************************************************
     * \brief   Location of the particle.
     **********************************************************************************************/
    std::vector<unsigned int> fieldIndices;
  };
  /**
   * @brief vector of particles (connected components)
   * contains
   */
  std::vector<Particle> particles;
  /*!***********************************************************************************************
   * \brief   Random seed.
   ************************************************************************************************/
  unsigned int rand_seed;

  // -------------------------------------------------------------------------------------------------
  // EVALUATING SECTION STARTS HERE (TODO maybe outsource in different file)
  // -------------------------------------------------------------------------------------------------

  /**
   * @brief Compares bulks of domain A and B
   *
   * @param domain_a
   * @param domain_b
   * @return amount of cells which alternated from or to 0
   */
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
  /**
   * @brief Gives size of each particle in a sorted way
   *
   * @return vector of sizes
   */
  std::vector<unsigned int> particle_size_distribution()
  {
    // Alternative
    //  findAggregates();
    //  std::vector<unsigned int> distribution;
    //  distribution.resize(particles.size());
    //  for(unsigned int i = 0; i < particles.size();i++)
    //  {
    //      distribution.at(i) = particles.at(i).fieldIndices.size();
    //  }

    constexpr unsigned int dim = nx.size();
    unsigned int fluids_size, field, neigh_field;
    std::vector<unsigned int> found_solids, distribution;

    std::for_each(domainFields.begin(), domainFields.end(),
                  [](unsigned int& field) { field = (field == 0); });

    for (auto first_fluid = std::find(domainFields.begin(), domainFields.end(), 0);
         first_fluid != domainFields.end();
         first_fluid = std::find(first_fluid, domainFields.end(), 0))
    {
      found_solids = std::vector<unsigned int>(1, std::distance(domainFields.begin(), first_fluid));
      domainFields[found_solids[0]] = CAM::uint_max;
      fluids_size = 1;
      for (unsigned int k = 0; k < fluids_size; ++k, fluids_size = found_solids.size())
      {
        field = found_solids[k];
        for (unsigned int i = 0; i < 2 * dim; ++i)
        {
          neigh_field = aim<nx>(field, direct_neigh<nx>(i));
          if (domainFields[neigh_field] == 0)
          {
            domainFields[neigh_field] = uint_max;
            found_solids.push_back(neigh_field);
          }
        }
      }
      distribution.push_back(found_solids.size());
    }
    std::sort(distribution.begin(), distribution.end());
    return distribution;
  }
  /**************************************************************************************************
   * @brief Gives the number of particles (connected component) inside the domain
   *
   * @return amount of particles
   ***************************************************************************************************/
  unsigned int n_solid_comp() { return particle_size_distribution().size(); }
  /*!*********************************************************************************************
   * \brief   Counts surfaces of all particles.
   *
   * \retval  n_surfaces      Surfaces of all particles.
   **********************************************************************************************/
  unsigned int getNSurfaces(std::vector<unsigned int> _fields)
  {
    unsigned int n_surfaces = 0;
    std::for_each(_fields.begin(), _fields.end(),
                  [&](const unsigned int field)
                  {
                    for (unsigned int i = 0; i < 2 * dim; ++i)
                      if (domainFields[aim<nx>(field, direct_neigh<nx>(i))] == 0)
                        ++n_surfaces;
                  });
    return n_surfaces;
  }
  /*!*********************************************************************************************
   * \brief   Find the longest shortest path inside domain.
   *
   * \retval  max_distance     Length of the path.
   **********************************************************************************************/
  unsigned int max_min_distance(Particle _particle)
  {
    std::vector<unsigned int> starts(1, _particle.fieldIndices[0]);
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
        for (unsigned int j = 0; j < 2 * dim; ++j)
        {
          neigh_field = aim<nx>(starts[i], direct_neigh<nx>(j));
          if (std::find(_particle.numbers.begin(), _particle.numbers.end(),
                        domainFields[neigh_field]) != _particle.numbers.end() &&
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
  double max_dimension_ratio(Particle _particle)
  {
    std::array<unsigned int, dim> width_dim;
    for (unsigned int i = 0; i < dim; ++i)
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
  unsigned int max_min_distance(unsigned int field, Particle _particle)
  {
    unsigned int max_distance = 0;
    fields_array_t domain_fields = domainFields;
    std::vector<unsigned int> visits(1, field);
    domain_fields[field] = uint_max;

    unsigned int neigh_field, visits_size = visits.size();
    for (unsigned int i = 0; i < visits_size; ++i, visits_size = visits.size())
    {
      field = visits[i];
      for (unsigned int j = 0; j < 2 * dim; ++j)
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
  //     /*!*********************************************************************************************
  //      * \brief   Find the longest amount of moves in certain direction inside particle.
  //      *
  //      * \param   dir_dim          Certain dimension
  //      * \retval  max_distance     Amount of moves
  //      **********************************************************************************************/
  unsigned int directed_max_min_distance(unsigned int dir_dim, Particle _particle)
  {
    unsigned int max_distance = 0, field = _particle.fieldIndices[0];
    fields_array_t domain_fields = domainFields;
    std::vector<unsigned int> visits(1, field);
    domain_fields[field] = uint_max - n_fields_;
    CAM::fieldNumbers_t min_val = domain_fields[field], max_val = domain_fields[field];

    unsigned int neigh_field, visits_size = visits.size();
    for (unsigned int i = 0; i < visits_size; ++i, visits_size = visits.size())
    {
      field = visits[i];
      for (unsigned int j = 0; j < 2 * dim; ++j)
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
  /**
   * @brief Gives the avarage particle size
   *
   * @return avarage particle size
   */
  constexpr double average_particle_size()
  {
    unsigned int n_solids = 0;
    for (unsigned int i = 0; i < domainFields.size(); ++i)
      n_solids += (domainFields[i] != 0);
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
    findAggregates();
    unsigned int n_single_cells =
      std::count_if(particles.begin(), particles.end(),
                    [](Particle particle) -> bool { return particle.fieldIndices.size() == 1; });

    unsigned int n_particles = particles.size();

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

    std::for_each(particles.begin(), particles.end(),
                  [&](Particle particle)
                  {
                    n_solids_part = particle.fieldIndices.size();
                    n_surfaces_part = getNSurfaces(particle.fieldIndices);

                    n_solids += n_solids_part;
                    n_surfaces += n_surfaces_part;

                    local_sphericity =
                      std::pow((double)n_solids_part, (double)(dim - 1) / (double)dim) /
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
      std::pow((double)n_surfaces, ((double)dim / (double)(dim - 1))) / (double)n_solids;

    double variance_particle_sizes = 0;
    for (unsigned int i = 0; i < n_particles; ++i)
    {
      variance_particle_sizes += std::pow(particles[i].fieldIndices.size() - mean_particle_size, 2);
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
  // /*!***********************************************************************************************
  //  * \brief   Computes connected fluid areas.
  //  *
  //  * \retval  n_fluid_comp      1st entry of array is connected fluid areas.
  //  *                            2nd entry of array is periodic connected fluid areas.
  //  ************************************************************************************************/
  std::array<unsigned int, 2> n_fluid_comp()
  {
    unsigned int n_connected_fluids = 0;
    unsigned int n_periodic_fluids = 0;
    unsigned int fluids_size, field, neigh_field;
    std::vector<unsigned int> found_fluids;

    bool periodic;
    std::vector<std::array<int, dim>> ref_dist(n_fields_);
    std::for_each(ref_dist.begin(), ref_dist.end(),
                  [](std::array<int, dim> dist) { dist.fill(0); });

    for (auto first_fluid = std::find(domainFields.begin(), domainFields.end(), 0);
         first_fluid != domainFields.end();
         first_fluid = std::find(first_fluid, domainFields.end(), 0))
    {
      periodic = false;
      found_fluids = std::vector<unsigned int>(1, std::distance(domainFields.begin(), first_fluid));
      domainFields[found_fluids[0]] = uint_max;
      fluids_size = 1;
      for (unsigned int k = 0; k < fluids_size; ++k, fluids_size = found_fluids.size())
      {
        field = found_fluids[k];
        for (unsigned int i = 0; i < 2 * dim; ++i)
        {
          neigh_field = aim<nx>(field, direct_neigh<nx>(i));
          if (domainFields[neigh_field] == 0)
          {
            domainFields[neigh_field] = uint_max;
            found_fluids.push_back(neigh_field);
            ref_dist[neigh_field] = ref_dist[field];
            ref_dist[neigh_field][i / 2] += 2 * (i % 2) - 1;
          }
          else if (domainFields[neigh_field] == uint_max)
            for (unsigned int j = 0; j < dim; ++j)
              if ((unsigned int)std::abs(ref_dist[field][j] - ref_dist[neigh_field][j]) > nx[j] - 2)
                periodic = true;
        }
      }
      ++n_connected_fluids;
      n_periodic_fluids += periodic;
    }

    std::for_each(domainFields.begin(), domainFields.end(),
                  [](CAM::fieldNumbers_t& field)
                  {
                    if (field == uint_max)
                      field = 0;
                  });

    return {n_connected_fluids, n_periodic_fluids};
  }

  // -------------------------------------------------------------------------------------------------
  // PRINTING SECTION STARTS HERE
  // -------------------------------------------------------------------------------------------------

  /*!*************************************************************************************************
   * \brief   Prints array of domain. Only used in print_n_dim().
   *
   * \param   fields      Particle index and size
   * \param   nx             Nx dimension and size
   * \param   init_index     Temporary index in current slice
   **************************************************************************************************/
  void print_1d(const fields_array_t& _fields, unsigned int init_index = 0)
  {
    const unsigned int len_numbers = std::log10(_fields.size() - 1) + 1;
    unsigned int index;

    for (unsigned int x = 0; x < nx[0]; ++x)
    {
      index = init_index + x;
      if (_fields[index] == 0)
        std::cout << std::setw(len_numbers) << std::setfill('0') << _fields[index] << "  ";
      else
        std::cout << "\033[0;31m" << std::setw(len_numbers) << std::setfill('0') << _fields[index]
                  << "\033[0m  ";
    }
    std::cout << std::endl;
  }
  /*!*************************************************************************************************
   * \brief   Prints array of fields in n dimensions. Only used in print_array().
   *
   * \tparam  n_dim          Temporary dimension of current slice
   * \param   fields      Particle index and size
   * \param   nx             Nx dimension and size
   * \param   init_index     Temporary index in current slice
   **************************************************************************************************/
  template <unsigned int n_dim>
  void print_n_dim(const fields_array_t& _fields, unsigned int init_index = 0)
  {
    if constexpr (n_dim == 1)
      print_1d(_fields, init_index);
    else
    {
      if (n_dim == 2 && nx.size() > 2)
      {
        unsigned int coord_i_helper = init_index;
        std::cout << std::endl;
        for (unsigned int i = 3; i < nx.size() + 1; ++i)
        {
          coord_i_helper = coord_i_helper / nx[i - 2];
          std::cout << i;
          if (i % 10 == 1 && i != 11)
            std::cout << "st";
          else if (i % 10 == 2 && i != 12)
            std::cout << "nd";
          else if (i % 10 == 3 && i != 13)
            std::cout << "rd";
          else
            std::cout << "th";
          std::cout << " coord: " << coord_i_helper % nx[i - 1] << "  ";
        }
        std::cout << std::endl;
      }
      for (unsigned int i = 0; i < nx[n_dim - 1]; ++i)
        print_n_dim<n_dim - 1>(_fields, (init_index + i) * nx[n_dim - 2]);
    }
  }
  /*!*************************************************************************************************
   * \brief   Runs print_n_dim() function
   *
   * \param   fields      Particle index and size
   * \param   nx             Nx dimension and size
   **************************************************************************************************/

  void print_array()
  {
    unsigned int sum = 0;
    std::for_each(domainFields.begin(), domainFields.end(),
                  [&](CAM::fieldNumbers_t t) { sum += t; });
    std::cout << "summe " << sum << std::endl;

    print_n_dim<nx.size()>(domainFields);
  }
};
}  // namespace CAM
#endif  // DOMAIN_HXX