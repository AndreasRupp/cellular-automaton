from __future__ import print_function

import matplotlib.pyplot as plt
import numpy as np
import os, sys
from datetime import datetime
  

# --------------------------------------------------------------------------------------------------
# Parameter identification test.
# --------------------------------------------------------------------------------------------------
def ecdf_identify(nx, porosity, n_steps, jump_parameter, ecdf_type, subset_sizes, n_choose_bins,
  min_value_shift, max_value_shift, jump_params, distance_fct, debug_mode, file_name, is_plot):
  start_time = datetime.now()
  print("Starting time is", start_time)

  try:
    import CAM
  except (ImportError, ModuleNotFoundError) as error:
    sys.path.append(os.path.dirname(os.path.abspath(__file__)) + os.sep  + ".." + os.sep + "import")
    import CAM

  try:
    import ecdf_estimator as ecdf
  except (ImportError, ModuleNotFoundError) as error:
    print("No installed ecdf_estimator package found! Using local ecdf_estimator.")
    sys.path.append(os.path.dirname(os.path.abspath(__file__)) + os.sep + ".." + os.sep + ".." +
      os.sep + "cil_estimator.git")
    # sys.path.append(os.path.dirname(os.path.abspath(__file__)) + os.sep + ".." + os.sep +
    #   "submodules" + os.sep + "ecdf_estimator.git")
    import ecdf_estimator as ecdf

  const            = CAM.config()
  const.nx         = nx
  const.debug_mode = debug_mode
  PyCAM            = CAM.include(const)

  # ------------------------------------------------------------------------------------------------
  def run_cam(jump_parameter, nx, porosity, n_steps, debug_mode=False):
    CAM_wrapper = PyCAM(porosity, jump_parameter)
    for step in range(n_steps):  CAM_wrapper.move_particles()
    return CAM_wrapper.fields()
  # ------------------------------------------------------------------------------------------------

  n_fields, n_iter = np.prod(nx), np.sum(subset_sizes)
  data = [ run_cam(jump_parameter, nx, porosity, n_steps, debug_mode) for _ in range(n_iter) ]

  end_time = datetime.now()
  print("CAM data acquired at", end_time, "after", end_time-start_time)

  fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(18, 5))

  if ecdf_type == "standard":
    min_val, max_val, _ = ecdf.estimate_radii_values(data[0:subset_sizes[0]],
      data[subset_sizes[0]:subset_sizes[0]+subset_sizes[1]], distance_fct )
    bins = np.linspace(min_val, max_val, 50)
    func = ecdf.estimator(data, bins, distance_fct, subset_sizes)
    ax[0,0] = ecdf.plot_ecdf_vectors(func, ax[0,0])
    ax[0,0] = ecdf.plot_mean_vector(func, ax[0,0], 'k.')
    ax[1,0] = ecdf.plot_chi2_test(func, ax[1,0], 20)
    func.choose_bins(n_choose_bins, min_value_shift, max_value_shift)
  elif ecdf_type == "bootstrap":
    min_val, max_val, _ = ecdf.estimate_radii_values(data[0:subset_sizes[0]],
      data[subset_sizes[0]:subset_sizes[0]+subset_sizes[1]], distance_fct )
    bins = np.linspace(min_val, max_val, 50)
    func = ecdf.bootstrap_estimator(data[0:subset_sizes[0]],
            data[subset_sizes[0]:subset_sizes[0]+subset_sizes[1]], bins, distance_fct)
    ax[0,0] = ecdf.plot_ecdf_vectors(func, ax[0,0])
    ax[0,0] = ecdf.plot_mean_vector(func, ax[0,0], 'k.')
    ax[1,0] = ecdf.plot_chi2_test(func, ax[1,0], 20)
    func.choose_bins(n_choose_bins, min_value_shift, max_value_shift)
  elif isinstance(ecdf_type,list):
    if not len(distance_fct) == len(n_choose_bins) and len(n_choose_bins) == len(ecdf_type):
      print("ERROR: Same amount of distance, bin, and type choices needed.")
    func_list = []
    for index in range(len(distance_fct)):
      distance, n_bins = distance_fct[index], n_choose_bins[index]
      min_val, max_val, _ = ecdf.estimate_radii_values(data[0:subset_sizes[0]],
        data[subset_sizes[0]:subset_sizes[0]+subset_sizes[1]], distance )
      bins = np.linspace(min_val, max_val, 50)
      if ecdf_type[index] == "standard":
        aux_func = ecdf.estimator(data, bins, distance, subset_sizes)
      elif ecdf_type[index] == "bootstrap":
        aux_func = ecdf.bootstrap_estimator(data[0:subset_sizes[0]],
                   data[subset_sizes[0]:subset_sizes[0]+subset_sizes[1]], bins, distance)
      aux_func.choose_bins(n_choose_bins[index], min_value_shift, max_value_shift)
      func_list.append(aux_func)
    func = ecdf.multiple_estimator( func_list )
  else:
    print("ERROR: Type of ecdf is invalid.")

  ax[0,0] = ecdf.plot_ecdf_vectors(func, ax[0,0], 'r.')
  ax[0,0] = ecdf.plot_mean_vector(func, ax[0,0], 'k.')
  ax[1,1] = ecdf.plot_chi2_test(func, ax[1,1])

  end_time = datetime.now()
  print("Objective function setup at", end_time, "after", end_time-start_time)

  values = [ func.evaluate( [ run_cam(jump_param, nx, porosity, n_steps, debug_mode) \
             for _ in range(subset_sizes[0]) ] ) for jump_param in jump_params ]

  ax[0,1].plot(jump_params, values, 'ro')

  if not os.path.exists('output'):  os.makedirs('output')
  plt.savefig('output/'+ file_name + '.png')

  end_time = datetime.now()
  print("Program ended at", end_time, "after", end_time-start_time)

  if is_plot:  plt.show()


def run_test_from_class(test_class):
  ecdf_identify(test_class.nx, test_class.porosity, test_class.n_steps, test_class.jump_parameter,
    test_class.ecdf_type, test_class.subset_sizes, test_class.n_choose_bins,
    test_class.min_value_shift, test_class.max_value_shift, test_class.jump_params,
    test_class.distance_fct, test_class.debug_mode, test_class.file_name, test_class.is_plot)
