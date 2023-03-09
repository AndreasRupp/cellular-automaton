#pragma once

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <chrono>

#include <CAM/building_units.hxx>
namespace CAM
{
    /*!*************************************************************************************************
 * \brief   Calculates the size of the domain.
 *
 * \retval n_field        size of the domain
 **************************************************************************************************/
template <auto nx>
static constexpr unsigned int n_fields()
{
unsigned int n_field = 1;
for (unsigned int i = 0; i < nx.size(); ++i)
    n_field *= nx[i];
return n_field;
}

template <auto nx>
static constexpr int direct_neigh(const unsigned int index)
{
static_assert(nx.size() != 0, "Dimension of zero does not make sense.");
int direct_neigh = (index % 2 == 0) ? -1 : 1;
for (unsigned int i = 0; i < index / 2; ++i)
    direct_neigh *= nx[i];
return direct_neigh;
}
    /*!*************************************************************************************************
 * \brief   Find field if one moves from position to move.
 *
 * \param   position  Current position of field that may move.
 * \param   move      Index shift induced by possible move.
 * \retval  index     Index of move target.
 **************************************************************************************************/
template <auto nx>
static constexpr unsigned int aim(const unsigned int position, const int move)
{
  unsigned int coord, new_pos = 0;
  for (unsigned int i = 0; i < nx.size(); ++i)
  {
    coord = (position / direct_neigh<nx>(2 * i + 1) + move / (int)direct_neigh<nx>(2 * i + 1) +
            n_fields<nx>()) %
            nx[i];
    new_pos += coord * direct_neigh<nx>(2 * i + 1);
  }
  return new_pos;
}



template <auto nx, typename fields_array_t = std::array<unsigned int, n_fields<nx>()>>
class Domain
{
  public:
    static constexpr unsigned int n_fields_ = n_fields<nx>();
    Domain()
    {
     if constexpr (std::is_same<fields_array_t,
                               std::vector<typename fields_array_t::value_type> >::value)
      domainFields.resize(n_fields_, 0);
    else
    {
      static_assert(
        std::is_same<fields_array_t,
                     std::array<typename fields_array_t::value_type, n_fields<nx>()> >::value,
        "The fields array has incorrect size");
      domainFields.fill(0);
    }
    }
  
    bool placeBU()//std::vector<int> _location , CAM::BuildingUnit* _unit
    {
        //_domainFields[aim<nx>(aiming, direct_neigh<nx>(i))] != _BU->number
        std::vector<unsigned int> vect(1, 10);
        buildingUnits.push_back(new ParticleBU(1, vect));
        domainFields[10] = 1;
        std::vector<unsigned int> vect1(1, 14);
        buildingUnits.push_back(new ParticleBU(2, vect1));
        domainFields[14] = 2;

        return true;
        
    }
    // bool placeBURandomly(CAM::BuildingUnit* _unit, unsigned int rand_seed = 0)
    // {
    //     if (random_seed == 0)
    //     {
    //         rand_seed = std::chrono::system_clock::now().time_since_epoch().count();
    //     }
    //     else
    //         rand_seed = random_seed;
        
    //     std::srand(rand_seed);

    //     // unsigned int position = std::rand() % (n_fields_);
    //     // for (unsigned int i = 0; i < n_particles; ++i)
    //     // {
    //     // while (fields_[position] != 0)
    //     //     position = std::rand() % (n_fields_);

    //     // particles_.push_back(particle(position, *this, i + 1));
    // }
    
  // -----------------------------------------------------------------------------------------------
 
 
//  /*!***********************************************************************************************
//    * \brief   Jump parameter.
//    ************************************************************************************************/
//   const double jump_parameter_singles_;

//   const double jump_parameter_composites_;
//  /*!***********************************************************************************************
//    * \brief   Number of particles.
//    ************************************************************************************************/
//   unsigned int n_particles;
  /*!***********************************************************************************************
   * \brief   Array of particle locations.
   ************************************************************************************************/
  fields_array_t domainFields;
  /*!***********************************************************************************************
   * \brief   Vector of particles.
   ************************************************************************************************/
  std::vector<CAM::BuildingUnit*> buildingUnits;
  /*!***********************************************************************************************
   * \brief   Random seed.
   ************************************************************************************************/
  unsigned int rand_seed;

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
    print_n_dim<nx.size()>(domainFields);
    }


};
}