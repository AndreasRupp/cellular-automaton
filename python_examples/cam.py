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
    from plot import plot, plot_to_file

  
  const                 = CAM.config()
  const.nx              = [30, 30, 30]
  const.debug_mode      = debug_mode
  jump_parameter_composites  = 10

  PyCAM = CAM.include(const)
  Domain = PyCAM(jump_parameter_composites)
  #Domain.placeSingleCellBURandomly(0.75, 5, 0)
  # success = 0
  # while success < 10:
  #   success = success + Domain.placeSphere( -1, 3, 5)
  # print(success)
  # success = 0
  # while success < 10:
  #   success = success + Domain.placePlane( -1, [1,1,4], 5)
  # print(success)
  #Domain.placeSphere( 1, 0)
  Domain.placeParticles()
  #Domain.print_array()
  save_data = np.zeros( (n_steps + 1, np.prod(const.nx)) ) 
  save_data[0] = Domain.fields()
  
  for step in range(n_steps):
    Domain.doCAM()
    save_data[step+1] = Domain.fields()
  print(Domain.average_particle_size())
  print(Domain.bulk_distance(save_data[0],save_data[-1]))
  end_time = datetime.now() 
  print("Program ended at", end_time, "after", end_time-start_time)
  
  plot_to_file(const.nx, save_data[-1], 'cam.png')
  plot(const.nx, save_data, 0)
  

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
