from __future__ import print_function

from scipy.stats import chi2
import matplotlib.pyplot as plt
from datetime import datetime

# from cil_estimator import estimator as cil

import os, sys
  

# --------------------------------------------------------------------------------------------------
# Function main.
# --------------------------------------------------------------------------------------------------
def main(debug_mode):
  start_time = datetime.now()
  print("Starting time is", start_time)

  try:
    import CAM
  except (ImportError, ModuleNotFoundError) as error:
    sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../import")
    import CAM

  try:
    import cil_estimator as cil
  except (ImportError, ModuleNotFoundError) as error:
    print("No installed cil_estimator package found! Using local cil_estimator.")
    sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../submodules/cil_estimator.git")
    import cil_estimator as cil

  nx = 50                                 # Size of considered domain!

  const            = CAM.config()
  const.nx         = [nx, nx]
  const.debug_mode = debug_mode
  PyCAM            = CAM.include(const)

  # ------------------------------------------------------------------------------------------------
  def run_cam(jump_size, nx, debug_mode=False):
    n_steps = 5

    const                 = CAM.config()
    const.nx              = [nx, nx]
    const.debug_mode      = debug_mode

    PyCAM = CAM.include(const)
    CAM_wrapper = PyCAM( 0.3, jump_size )

    for step in range(n_steps):
      CAM_wrapper.move_particles()

    return CAM_wrapper.fields()
  # ------------------------------------------------------------------------------------------------

  n_iter = 200
  data   = [[0] * (nx * nx)] * n_iter
  values = [0] * 10
  lower_bound = int(0.3*nx*nx)
  upper_bound = int(0.3*nx*nx)
  step        = 1
  radii  = range(int(0.3*nx*nx),int(0.55*nx*nx),1)

  for iter in range(n_iter):
    data[iter] = run_cam(5, nx, debug_mode)

  end_time = datetime.now()
  print("CAM data acquired at", end_time, "after", end_time-start_time)

  func = cil.estimator(data, radii, PyCAM.bulk_distance, [10] * 20)

  fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(18, 5))
  ax[0] = cil.plot_correlation_vectors(func, ax[0])
  ax[0] = cil.plot_mean_vector(func, ax[0])

  func.choose_radii(10, min_value_shift=0.1, max_value_shift=-0.1)
  
  ax[0] = cil.plot_correlation_vectors(func, ax[0], 'r.')
  
  end_time = datetime.now()
  print("Objective function setup at", end_time, "after", end_time-start_time)

  for jump_size in range(len(values)):
    for iter in range(n_iter):
      data[iter] = run_cam(jump_size, nx, debug_mode)
    values[jump_size] = func.evaluate( data )

  
  end_time = datetime.now()
  print("Program ended at", end_time, "after", end_time-start_time)
  

  ax[1].plot(range(len(values)), values, 'ro')
  plt.show()


# --------------------------------------------------------------------------------------------------
# Define main function.
# -------------------------------------------------------------------------------------------------- 
if __name__ == "__main__":
  main(len(sys.argv) > 1 and sys.argv[1] == "True")