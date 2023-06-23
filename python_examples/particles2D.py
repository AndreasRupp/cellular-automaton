from __future__ import print_function
from datetime import datetime

import os, sys
import numpy as np
import itertools


def cam_particles(n_steps, jump_parameter, distribution,porosity, debug_mode=False):
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
  const.nx              = [500, 500]
  const.debug_mode      = debug_mode
  PyCAM = CAM.include(const)
  Domain = PyCAM(jump_parameter)

  n_fields = np.prod(const.nx)
  print(n_fields)
  n_solids = (1 - porosity) * n_fields
  
  placement_i = n_solids * distribution
  placement_g = n_solids * (1- distribution)
  size_g = 2*30
  size_i = 2*6

  print(placement_i/size_i)
  print(placement_g/size_g)
  n_g = placement_g/size_g
  n_i = placement_g/size_i
   #goethite
  summe = 0
  success = 0
  while success < n_g:
    success = success + Domain.place_plane(jump_parameter, [2,30], -1, 1)
  print(success)
  summe = summe + success
  # success = 0
  # while success < n_g:
  #   success = success + Domain.place_plane(jump_parameter, [30,2], -1, -1)
  # print(success)
  # summe = summe + success
  #illite
  success = 0
  while success < n_i:
    success = success + Domain.place_plane(jump_parameter, [2,6], -1, -1)
  print(success)
  summe = summe + success
  # success = 0
  # while success < n_i:
  #   success = success + Domain.place_plane(jump_parameter, [6,2], -1, 1)
  # print(success)
  # summe = summe + success


  # Domain.print_array()

  # Domain.place_particles()
    #   //  constexpr std::array<unsigned int, stencil_size<nx,5>()> possible_moves =
    # //   CAM::get_stencil_a<nx, 5>();
  save_data = np.zeros( (n_steps + 1, np.prod(const.nx)) ) 
  save_data[0] = Domain.fields()
  # plot_to_vtk("output/cam", [save_data[0]], const.nx)
  for step in range(n_steps):
    Domain.do_cam()
    save_data[step+1] = Domain.fields()

  # print(Domain.average_particle_size())
  # print(Domain.bulk_distance(save_data[0],save_data[-1]))
  end_time = datetime.now() 
  print("Program ended at", end_time, "after", end_time-start_time)
  save_data[(save_data < n_g) & (save_data > 0)] = 1
  save_data[save_data >= n_g] = 2
  # if not os.path.exists('output'):  os.makedirs('output')
  # plot_to_vtk("output/cam", [save_data[0]], const.nx)
  if not os.path.exists('output2D'):  os.makedirs('output2D')
  plot_to_vtk("output/cam", save_data, const.nx)
  # if not os.path.exists('output1'):  os.makedirs('output1')
  # plot_to_vtk("output1/cam", [save_data[-1]], const.nx)
  print("done")
  # plot_to_file(const.nx, save_data[-1], 'output/cam.png')
  plot(const.nx, save_data, 0)
# --------------------------------------------------------------------------------------------------
# Function main.
# --------------------------------------------------------------------------------------------------
def main(debug_mode):
  n_steps =50
  cam_particles(n_steps,5,0.75, 0.9,debug_mode)


# --------------------------------------------------------------------------------------------------
# Define main function.
# -------------------------------------------------------------------------------------------------- 
if __name__ == "__main__":
  main(len(sys.argv) > 1 and sys.argv[1] == "True")
