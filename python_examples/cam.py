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
  const.nx              = [25, 25,25]
  const.debug_mode      = debug_mode

  PyCAM = CAM.include(const)
  Domain = PyCAM()
  #Domain.placeSingleCellBURandomly(0.75, 5, 0)
  #Domain.placeSphere( 1, 0)
  Domain.placeBU()
  #Domain.print_array()
  save_data = np.zeros( (n_steps + 1, np.prod(const.nx)) ) 
  save_data[0] = Domain.fields()
  
  for step in range(n_steps):
    Domain.doCAM()
    save_data[step+1] = Domain.fields()
  
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
