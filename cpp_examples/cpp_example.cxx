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
 *    jump_param      how far individual particles are allowed to jump
 **************************************************************************************************/
int main()
{
  constexpr std::array<unsigned int, 2> nx = {10, 10};
  const unsigned int n_moves = 5;
  const double jump_param = 1.;

  CAM::CAMInterface<nx> CAM;

  const double porosity = 0.5;
  CAM.place_singleCellBU_randomly(porosity, jump_param);

  // std::vector<unsigned int> extent = {3, 1};
  // std::cout << "is placed? " << CAM.place_plane(-1, extent, jump_param) << std::endl;
  CAM.print_array();
  std::cout << std::endl << std::endl;
  for (unsigned int i = 0; i < n_moves; ++i)
  {
    CAM.do_CAM();
    CAM.print_array();
    std::cout << std::endl << std::endl;
  }

  std::cout << std::endl << "Characteristics / Measures:" << std::endl;
  const std::array<double, 12> meas = CAM.eval_measures();
  for (unsigned int k = 0; k < 12; ++k)
    std::cout << "Meas[" << k << "] = " << meas[k] << std::endl;
}
