from __future__ import print_function
from datetime import datetime

import os, sys
import numpy as np
import itertools


def cam_particles(n_steps, jump_parameter, distribution, debug_mode=False):
  start_time = datetime.now()
  print("Starting time is", start_time)

  try:
    import CAM
  except (ImportError, ModuleNotFoundError) as error:
    sys.path.append(os.path.dirname(os.path.abspath(__file__)) + os.sep  + ".." + os.sep + "import")
    import CAM

  try:
    from plot import plot, plot_to_file, plot_to_vtk
  except (ImportError, ModuleNotFoundError) as error:
    sys.path.append(os.path.dirname(os.path.abspath(__file__)) + os.sep  + ".." + os.sep + 
      "python_functions")
    from plot import plot, plot_to_file, plot_to_vtk

  
  const                 = CAM.config()
  const.nx              = [250, 250, 250]
  const.debug_mode      = debug_mode
  PyCAM = CAM.include(const)
  Domain = PyCAM(jump_parameter)

  n_fields = np.prod(const.nx)
  print(n_fields)
  n_solids = 0.1 * n_fields
  each_placement = n_solids/6
  print(each_placement)
  size_g = 2*2*30
  size_i = 2*2*6

  print(each_placement/size_i)
  print(each_placement/size_g)
  n_g = each_placement/size_g
  n_i = each_placement/size_i
   #goethite
  summe = 0
  success = 0
  while success < n_g:
    success = success + Domain.place_plane(jump_parameter, [2,2,30], -1)
  print(success)
  summe = summe + success
  success = 0
  while success < n_g:
    success = success + Domain.place_plane(jump_parameter, [2,30,2], -1)
  print(success)
  summe = summe + success
  success = 0
  while success <n_g:
    success = success + Domain.place_plane(jump_parameter, [30,2,2], -1)
  print(success)
  summe = summe + success
  print(summe)
  #illite
  success = 0
  while success < n_i:
    success = success + Domain.place_plane(jump_parameter, [2,2,6], -1)
  print(success)
  summe = summe + success
  success = 0
  while success < n_i:
    success = success + Domain.place_plane(jump_parameter, [2,6,2], -1)
  print(success)
  summe = summe + success
  success = 0
  while success <  n_i:
    success = success + Domain.place_plane(jump_parameter, [6,2,2], -1)
  print(success)
  summe = summe + success


  # Domain.print_array()

  # Domain.place_particles()
    #   //  constexpr std::array<unsigned int, stencil_size<nx,5>()> possible_moves =
    # //   CAM::get_stencil_a<nx, 5>();
  save_data = np.zeros( (n_steps + 1, np.prod(const.nx)) ) 
  save_data[0] = Domain.fields()
  # plot_to_vtk("output/cam", [save_data[0]], const.nx)
  for step in range(n_steps):
    Domain.do_cam()
    print("done")
    save_data[step+1] = Domain.fields()

  # print(Domain.average_particle_size())
  # print(Domain.bulk_distance(save_data[0],save_data[-1]))
  end_time = datetime.now() 
  print("Program ended at", end_time, "after", end_time-start_time)
  
  # if not os.path.exists('output'):  os.makedirs('output')
  # plot_to_vtk("output/cam", [save_data[0]], const.nx)
  if not os.path.exists('output'):  os.makedirs('output')
  plot_to_vtk("output/cam", save_data, const.nx)
  # if not os.path.exists('output1'):  os.makedirs('output1')
  # plot_to_vtk("output1/cam", [save_data[-1]], const.nx)
  print("done")
  # plot_to_file(const.nx, save_data[-1], 'output/cam.png')
  #plot(const.nx, save_data, 0)
# --------------------------------------------------------------------------------------------------
# Function main.
# --------------------------------------------------------------------------------------------------
def main(debug_mode):
  n_steps =5
  cam_particles(n_steps,5,0.5 ,debug_mode)


# --------------------------------------------------------------------------------------------------
# Define main function.
# -------------------------------------------------------------------------------------------------- 
if __name__ == "__main__":
  main(len(sys.argv) > 1 and sys.argv[1] == "True")
