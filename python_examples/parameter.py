from __future__ import print_function

import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime

import os, sys
  

# --------------------------------------------------------------------------------------------------
# Parameter identification test.
# --------------------------------------------------------------------------------------------------
def parameter_identification_test(nx, porosity, n_steps, jump_parameter, n_iter, values, bins,
  n_choose_bins, subset_sizes, min_value_shift, max_value_shift, debug_mode):
  start_time = datetime.now()
  print("Starting time is", start_time)

  n_fields = np.prod(nx)
  data = [[0] * n_fields] * n_iter

  try:
    import CAM
  except (ImportError, ModuleNotFoundError) as error:
    sys.path.append(os.path.dirname(os.path.abspath(__file__)) + os.sep  + ".." + os.sep + "import")
    import CAM

  try:
    import ecdf_estimator as ecdf
  except (ImportError, ModuleNotFoundError) as error:
    print("No installed ecdf_estimator package found! Using local ecdf_estimator.")
    # sys.path.append(os.path.dirname(os.path.abspath(__file__)) + os.sep + ".." + os.sep + ".." +
    #   os.sep + "ecdf_estimator.git")
    sys.path.append(os.path.dirname(os.path.abspath(__file__)) + os.sep + ".." + os.sep +
      "submodules" + os.sep + "ecdf_estimator.git")
    import ecdf_estimator as ecdf

  const            = CAM.config()
  const.nx         = nx
  const.debug_mode = debug_mode
  PyCAM            = CAM.include(const)

  # ------------------------------------------------------------------------------------------------
  def run_cam(jump_parameter, nx, porosity, n_steps, debug_mode=False):
    const            = CAM.config()
    const.nx         = nx
    const.debug_mode = debug_mode

    PyCAM       = CAM.include(const)
    CAM_wrapper = PyCAM(porosity, jump_parameter)

    for step in range(n_steps):
      CAM_wrapper.move_particles()

    return CAM_wrapper.fields()
  # ------------------------------------------------------------------------------------------------

  for iter in range(n_iter):
    data[iter] = run_cam(jump_parameter, nx, porosity, n_steps, debug_mode)

  end_time = datetime.now()
  print("CAM data acquired at", end_time, "after", end_time-start_time)

  func = ecdf.estimator(data, bins, PyCAM.bulk_distance, subset_sizes)

  fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(18, 5))
  ax[0,0] = ecdf.plot_ecdf_vectors(func, ax[0,0])
  ax[0,0] = ecdf.plot_mean_vector(func, ax[0,0])
  ax[1,0] = ecdf.plot_chi2_test(func, ax[1,0], 20)

  func.choose_bins(n_choose_bins, min_value_shift, max_value_shift)
  
  ax[0,0] = ecdf.plot_ecdf_vectors(func, ax[0,0], 'r.')
  ax[1,1] = ecdf.plot_chi2_test(func, ax[1,1])
  
  end_time = datetime.now()
  print("Objective function setup at", end_time, "after", end_time-start_time)

  for jump_size in range(len(values)):
    for iter in range(n_iter):
      data[iter] = run_cam(jump_size, nx, porosity, n_steps, debug_mode)
    values[jump_size] = func.evaluate( data )

  end_time = datetime.now()
  print("Program ended at", end_time, "after", end_time-start_time)
  
  ax[0,1].plot(range(len(values)), values, 'ro')

  if not os.path.exists('output'):  os.makedirs('output')
  plt.savefig('output/prameter.png')
  plt.show()


# --------------------------------------------------------------------------------------------------
# Define main function.
# -------------------------------------------------------------------------------------------------- 
if __name__ == "__main__":

  # Configure the cellular automaton method (CAM).
  nx             = [50, 50]
  porosity       = 0.3
  n_steps        = 5
  jump_parameter = 5

  nf = np.prod(nx)  # This is the number of fields in the CAM and can be derived from nx.

  # Configure the eCDF method.
  n_iter          = 2000                                        # Number of dataset sizes.
  values          = [0] * 10                                    # Length defines tested jump_param.
  bins            = range(int(0.375*nf),int(0.475*nf),5)        # All radii that are checked.
  n_choose_bins   = 20                                          # Number of selected radii for algo.
  subset_sizes    = [40] * 50                                   # Multiplies to n_iter!
  min_value_shift = 0.1                                         # Cutoff value for small values.
  max_value_shift = -0.1                                        # 1 + "cutoff value for large val."

  debug_mode = len(sys.argv) > 1 and sys.argv[1] == "True"
  parameter_identification_test(nx, porosity, n_steps, jump_parameter, n_iter, values, bins,
    n_choose_bins, subset_sizes, min_value_shift, max_value_shift, debug_mode)
