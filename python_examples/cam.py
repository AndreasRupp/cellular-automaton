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
    from plot import plot, plot_update
  except (ImportError, ModuleNotFoundError) as error:
    sys.path.append(os.path.dirname(os.path.abspath(__file__)) + os.sep  + ".." + os.sep + 
      "python_functions")
    from plot import plot, plot_update
  
  const                 = CAM.config()
  const.nx              = [100, 100]
  const.debug_mode      = debug_mode

  PyCAM = CAM.include(const)
  CAM_wrapper = PyCAM( 0.9, 1 )

  save_data = np.zeros( (n_steps + 1, np.prod(const.nx)) ) 
  save_data[0] = CAM_wrapper.fields()

  for step in range(n_steps):
    data = CAM_wrapper.move_particles()
    save_data[step+1] = data
    # CAM_wrapper.print_array(data)
    # if step != n_steps-1: print()
  
  end_time = datetime.now()
  print("Program ended at", end_time, "after", end_time-start_time)
  
  plt = plot_update(const.nx, save_data[-1])
  plt.savefig('output/cam.png')

  plot(const.nx, save_data, 0)
  

# --------------------------------------------------------------------------------------------------
# Function main.
# --------------------------------------------------------------------------------------------------
def main(debug_mode):
  n_steps = 1000
  cam_test(n_steps, debug_mode)


# --------------------------------------------------------------------------------------------------
# Define main function.
# -------------------------------------------------------------------------------------------------- 
if __name__ == "__main__":
  main(len(sys.argv) > 1 and sys.argv[1] == "True")
