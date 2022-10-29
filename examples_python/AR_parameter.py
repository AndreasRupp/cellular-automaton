from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

import cil_estimator as cil

import os, sys
  

# --------------------------------------------------------------------------------------------------
# Function main.
# --------------------------------------------------------------------------------------------------
def main(debug_mode):
  try:
    import CAM
  except (ImportError, ModuleNotFoundError) as error:
    sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../import")
    import CAM

  # try:
  #   import cil_estimator as cil
  # except (ImportError, ModuleNotFoundError) as error:
  #   sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../../cil_estimator.git")
  #   import cil_estimator as cil

  nx = 50                                 # Size of considered domain!

  const            = CAM.config()
  const.nx         = [nx, nx]
  const.debug_mode = debug_mode
  PyCAM            = CAM.include(const)
  CAM_wrapper      = PyCAM( 0, 0 )

  # ------------------------------------------------------------------------------------------------
  def run_cam(jump_size, nx, debug_mode=False):
    n_steps = 5

    const                 = CAM.config()
    const.nx              = [nx, nx]
    const.debug_mode      = debug_mode

    PyCAM = CAM.include(const)
    CAM_wrapper = PyCAM( 0.7, jump_size )

    for step in range(n_steps):
      CAM_wrapper.move_particles()

    return CAM_wrapper.fields()
  # ------------------------------------------------------------------------------------------------

  start_time = datetime.now()
  print("Starting time is", start_time)

  n_iter = 100
  data   = [[0] * (nx * nx)] * n_iter
  values = [0] * 10
  radii  = range(int(0.4*nx*nx),int(0.5*nx*nx),2)

  for iter in range(n_iter):
    data[iter] = run_cam(5, nx, debug_mode)
  func = cil.estimator( data, radii, CAM_wrapper.bulk_distance, [20] * 5 )

  for jump_size in range(10):
    for iter in range(n_iter):
      data[iter] = run_cam(jump_size, nx, debug_mode)
    values[jump_size] = func.evaluate( data )

  end_time = datetime.now()
  print("Program ended at", end_time, "after", end_time-start_time)

  plt.plot(range(10), values, 'ro')
  plt.show()


# --------------------------------------------------------------------------------------------------
# Define main function.
# -------------------------------------------------------------------------------------------------- 
if __name__ == "__main__":
  main(len(sys.argv) > 1 and sys.argv[1] == "True")
