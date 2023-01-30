from __future__ import print_function

import matplotlib.pyplot as plt
import numpy as np
import os, sys
from datetime import datetime

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
  #   os.sep + "cil_estimator.git")
  sys.path.append(os.path.dirname(os.path.abspath(__file__)) + os.sep + ".." + os.sep +
    "submodules" + os.sep + "ecdf_estimator.git")
  import ecdf_estimator as ecdf
  

# --------------------------------------------------------------------------------------------------
# Parameter identification test.
# --------------------------------------------------------------------------------------------------
def ecdf_identify(nx, porosity, n_steps, jump_parameter, ecdf_type, subset_sizes, n_bins, n_runs,
  min_value_shift, max_value_shift, jump_params, distance_fct, debug_mode, file_name, is_plot):

  def run_cam(jump_parameter, nx, porosity, n_steps, debug_mode=False):
    CAM_wrapper = PyCAM(porosity, jump_parameter)
    for step in range(n_steps):  CAM_wrapper.move_particles()
    return CAM_wrapper.fields()

  def generate_ecdf(data, subset_sizes, distance_fct, n_bins, ecdf_type, ax=None):
    min_val, max_val, distance_data = ecdf.estimate_radii_values(data[0:subset_sizes[0]],
      data[subset_sizes[0]:subset_sizes[0]+subset_sizes[1]], distance_fct )
    bins = np.linspace(min_val, max_val, 50)
    if ax is not None:
      if ecdf_type == "standard":
        aux_func = ecdf.standard(data, bins, distance_fct, subset_sizes)
      elif ecdf_type == "bootstrap":
        aux_func = ecdf.bootstrap(data[0:subset_sizes[0]],
                     data[subset_sizes[0]:subset_sizes[0]+subset_sizes[1]], bins, distance_fct)
      ax = ecdf.plot_ecdf_vectors(aux_func, ax, 'm.')
    bins = ecdf.choose_bins(distance_data, bins, n_bins)
    if ecdf_type == "standard":
      func = ecdf.standard(data, bins, distance_fct, subset_sizes)
    elif ecdf_type == "bootstrap":
      func = ecdf.bootstrap(data[0:subset_sizes[0]],
               data[subset_sizes[0]:subset_sizes[0]+subset_sizes[1]], bins, distance_fct)
    return func, ax


  start_time = datetime.now()
  print("Starting time is", start_time)

  CAM_config            = CAM.config()
  CAM_config.nx         = nx
  CAM_config.debug_mode = debug_mode
  PyCAM                 = CAM.include(CAM_config)

  n_fields, n_iter = np.prod(nx), np.sum(subset_sizes)
  data = [ run_cam(jump_parameter, nx, porosity, n_steps, debug_mode) for _ in range(n_iter) ]

  end_time = datetime.now()
  print("CAM data acquired at", end_time, "after", end_time-start_time)

  fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(18, 5))

  if ecdf_type == "standard" or ecdf_type == "bootstrap":
    func, ax[0,0] = generate_ecdf(data, subset_sizes, distance_fct, n_bins, ecdf_type, ax[0,0])
  elif isinstance(ecdf_type,list):
    if not len(distance_fct) == len(n_bins) and len(n_bins) == len(ecdf_type):
      print("ERROR: Same amount of distance, bin, and type choices needed.")
    func_list = []
    for i in range(len(distance_fct)):
      aux_func, _ = generate_ecdf(data, subset_sizes, distance_fct[i], n_bins[i], ecdf_type[i])
      func_list.append( aux_func )
    func = ecdf.multiple( func_list )
  else:
    print("ERROR: Type of ecdf is invalid.")

  ax[0,0] = ecdf.plot_ecdf_vectors(func, ax[0,0], 'c.')
  ax[0,0] = ecdf.plot_mean_vector(func, ax[0,0], 'k.')
  ax[1,0] = ecdf.plot_chi2_test(func, ax[1,0])

  end_time = datetime.now()
  print("Objective function setup at", end_time, "after", end_time-start_time)

  means_log = [0.] * len(jump_params)
  means_nor = [0.] * len(jump_params)
  for _ in range(n_runs):
    values = [ ecdf.evaluate( func, [ run_cam(jump_param, nx, porosity, n_steps, debug_mode) \
             for _ in range(subset_sizes[0]) ] ) for jump_param in jump_params ]
    ax[0,1].plot(jump_params, values, 'ro')
    means_log = [ means_log[i] + values[i] / n_runs for i in range(len(jump_params))]
    values = [ np.exp(-0.5 * value)   for value in values ]
    values = [ value / np.sum(values) for value in values ]
    ax[1,1].plot(jump_params, values, 'ro')
    means_nor = [ means_nor[i] + values[i] / n_runs for i in range(len(jump_params))]

  ax[0,1].plot(jump_params, means_log, 'bo')
  ax[1,1].plot(jump_params, means_nor, 'bo')

  if not os.path.exists('output'):  os.makedirs('output')
  plt.savefig('output/'+ file_name + '.png')

  end_time = datetime.now()
  print("Program ended at", end_time, "after", end_time-start_time)

  if is_plot:  plt.show()


def run_test_from_class(test_class):
  ecdf_identify(test_class.nx, test_class.porosity, test_class.n_steps, test_class.jump_parameter,
    test_class.ecdf_type, test_class.subset_sizes, test_class.n_bins, test_class.n_runs,
    test_class.min_value_shift, test_class.max_value_shift, test_class.jump_params,
    test_class.distance_fct, test_class.debug_mode, test_class.file_name, test_class.is_plot)
