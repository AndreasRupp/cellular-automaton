#include <CAM/building_units.hxx>
#include <CAM/cam_interface.hxx>
#include <iostream>
/*!*************************************************************************************************
 * \brief   Main function.
 *
 * Runs CAM, updates and prints the matrix.
 *
 * Parameters:
 *    nx              the size of a row for each dimension of the domain
 *    n_moves         the number of iterations of the CAM
 *    porosity        the percentage of void space, not occupied by solid
 *    stencil_size      how far individual particles are allowed to jump
 **************************************************************************************************/
int main()
{
  constexpr std::array<unsigned int, 2> nx = {10, 10};
  const unsigned int n_moves = 5;
  const unsigned int stencil_size = 1.;
  CAM::CAMInterface<nx, stencil_size> CAM;
  const double porosity = 0.5;
  const double random_seed = 0;
  CAM.place_single_cell_bu_randomly(stencil_size, porosity, random_seed);

  // std::vector<unsigned int> extent = {3, 1};
  // std::vector<unsigned int> face_charges_pos = {1, 1, 1, 1};
  // std::vector<unsigned int> face_charges_neg = {-1, -1, -1, -1};
  // std::cout << "is placed? " << CAM.place_plane(stencil_size, extent, face_charges_pos) <<

  CAM.print_array();
  std::cout << std::endl << std::endl;
  for (unsigned int i = 0; i < n_moves; ++i)
  {
    CAM.do_cam();
    CAM.print_array();
    std::cout << std::endl << std::endl;
  }

  std::cout << std::endl << "Characteristics / Measures:" << std::endl;
  const std::vector<double> meas = CAM.eval_measures();
  for (unsigned int k = 0; k < 12; ++k)
    std::cout << "Meas[" << k << "] = " << meas[k] << std::endl;
}
