from __future__ import print_function
from datetime import datetime

import os, sys
import numpy as np


# --------------------------------------------------------------------------------------------------
# Function bilaplacian_test.
# --------------------------------------------------------------------------------------------------
def cam_test(n_steps, debug_mode=False):
  start_time = datetime.now()
  print("Starting time is", start_time)

  try:
    import CAM
  except (ImportError, ModuleNotFoundError) as error:
    sys.path.append(os.path.dirname(os.path.abspath(__file__)) + os.sep  + ".." + os.sep + "import")
    import CAM

  try:
    from plot import plot, plot_to_file
  except (ImportError, ModuleNotFoundError) as error:
    sys.path.append(os.path.dirname(os.path.abspath(__file__)) + os.sep  + ".." + os.sep + 
      "python_functions")
    from plot import plot, plot_to_file, plot_to_vtk

  
  const                 = CAM.config()
  const.nx              = [50, 50, 50]
  const.debug_mode      = debug_mode
  jump_parameter_composites  = 10
  jump_parameter = 5
  PyCAM = CAM.include(const)
  Domain = PyCAM(jump_parameter_composites)

 # Domain.placeSingleCellBURandomly(0.75, jump_parameter, 0)

  success = 0
  while success < 100:
    success = success + Domain.placeSphere( -1, 3, jump_parameter)
  #print("Nr of spheres " + str(success))
  success = 0
  while success < 100:
    success = success + Domain.placePlane( -1, [1,1,5], jump_parameter)
  #print("Nr of planes " + str(success))

  #Domain.placeParticles()
  
  save_data = np.zeros( (n_steps + 1, np.prod(const.nx)) ) 
  print(save_data.shape)
  save_data[0] = Domain.fields()

  for step in range(n_steps):
    Domain.doCAM()
    save_data[step+1] = Domain.fields()

  # print(Domain.average_particle_size())
  # print(Domain.bulk_distance(save_data[0],save_data[-1]))
  end_time = datetime.now() 
  print("Program ended at", end_time, "after", end_time-start_time)
  
  if not os.path.exists('output'):  os.makedirs('output')
  #plot_to_file(const.nx, save_data[-1], 'output/cam.png')
  filename = "output/cam"
  plot_to_vtk(filename, save_data, const.nx)
  



  #plot(const.nx, save_data, 0)
  
  

# --------------------------------------------------------------------------------------------------
# Function main.
# --------------------------------------------------------------------------------------------------
def main(debug_mode):
  n_steps =5
  cam_test(n_steps, debug_mode)


# --------------------------------------------------------------------------------------------------
# Define main function.
# -------------------------------------------------------------------------------------------------- 
if __name__ == "__main__":
  main(len(sys.argv) > 1 and sys.argv[1] == "True")
