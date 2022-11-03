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
  try:
    import CAM
  except (ImportError, ModuleNotFoundError) as error:
    sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../import")
    import CAM

  try:
    import cil_estimator as cil
  except (ImportError, ModuleNotFoundError) as error:
    sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../../cil_estimator.git")
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
  lower_bound = int(0.3*nx*nx)
  upper_bound = int(0.3*nx*nx)
  step        = 1
  radii  = range(int(0.3*nx*nx),int(0.5*nx*nx),1)
  # radii  = range(int(0.35*nx*nx),int(0.475*nx*nx),10)

  for iter in range(n_iter):
    data[iter] = run_cam(5, nx, debug_mode)
  func = cil.estimator( data, radii, PyCAM.bulk_distance, [20] * 5 )

  import numpy as np


  

  # plotter = cil.plot_correlation_vectors(func)
  # plotter = cil.plot_mean_vector(func)
  func.choose_radii(10)

  vec_mat = func.correlation_vector_matrix
  derivative = [sum(vec_mat[i+1]) -  sum(vec_mat[i]) for i in range(len(vec_mat)-1)]
  plt.plot(func.radii[:-1], derivative)
    # plt.hist(derivative, 10, cumulative=False, density=False, facecolor='g', alpha=0.75)
  plt.show()


  # plotter = cil.plot_correlation_vectors(func, plotter, 'r.')
  # plotter.show()

  # for jump_size in range(10):
  #   for iter in range(n_iter):
  #     data[iter] = run_cam(jump_size, nx, debug_mode)
  #   values[jump_size] = func.evaluate( data )

  end_time = datetime.now()
  print("Program ended at", end_time, "after", end_time-start_time)

  # print(stats.chisquare(data))

  # mat_of_correlation_vec = func.correlation_matrix()

  # rv = stats.chi2(mat_of_correlation_vec)
  # plt.plot(func.radii, rv.pdf(func.radii), 'k-', lw=2, label='frozen pdf')
  # r = chi2.rvs(func.radii)
  # plt.hist(r, density=True, histtype='stepfilled', alpha=0.2)
  # plt.legend(loc='best', frameon=False)
  # plt.show()
  

  # plt.plot(range(10), values, 'ro')
  # plt.show()


# --------------------------------------------------------------------------------------------------
# Define main function.
# -------------------------------------------------------------------------------------------------- 
if __name__ == "__main__":
  main(len(sys.argv) > 1 and sys.argv[1] == "True")
