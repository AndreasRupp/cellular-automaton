from __future__ import print_function

import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
import multiprocessing
import os, sys, time
  

# --------------------------------------------------------------------------------------------------
# Parameter identification test.
# --------------------------------------------------------------------------------------------------
def parameter_identification_test(nx, porosity, n_steps, jump_parameter,
  subset_sizes, n_choose_bins, min_value_shift, max_value_shift, jump_params, bins,
  distance_fct,  debug_mode, file_name, is_plot = 0):
  start_time = datetime.now()
  print("Starting time is", start_time)

  n_fields = np.prod(nx)
  n_iter   = np.sum(subset_sizes)
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

    for step in range(n_steps):
      CAM_wrapper.move_particles()

    return CAM_wrapper.fields()
  # ------------------------------------------------------------------------------------------------

  for iter in range(n_iter):
    data[iter] = run_cam(jump_parameter, nx, porosity, n_steps, debug_mode)

  end_time = datetime.now()
  print("CAM data acquired at", end_time, "after", end_time-start_time)

  min_val, max_val, _ = ecdf.estimate_radii_values(
    data[0:subset_sizes[0]], data[subset_sizes[0]:subset_sizes[0]+subset_sizes[1]], distance_fct)
  step_size = np.floor((max_val - min_val) / 50)
  bins = range(int(min_val), int(max_val), int(step_size))

  func = ecdf.estimator(data, bins, distance_fct, subset_sizes)

  fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(18, 5))
  ax[0,0] = ecdf.plot_ecdf_vectors(func, ax[0,0])
  ax[0,0] = ecdf.plot_mean_vector(func, ax[0,0], 'k.')
  ax[1,0] = ecdf.plot_chi2_test(func, ax[1,0], 20)

  func.choose_bins(n_choose_bins, min_value_shift, max_value_shift)
  
  ax[0,0] = ecdf.plot_ecdf_vectors(func, ax[0,0], 'r.')
  ax[0,0] = ecdf.plot_mean_vector(func, ax[0,0], 'k.')
  ax[1,1] = ecdf.plot_chi2_test(func, ax[1,1])
  
  end_time = datetime.now()
  print("Objective function setup at", end_time, "after", end_time-start_time)

  values = [0] * len(jump_params)
  data   = [0] * subset_sizes[0]
  for jump_index in range(len(values)):
    for iter in range(subset_sizes[0]):
      data[iter] = run_cam(jump_params[jump_index], nx, porosity, n_steps, debug_mode)
    values[jump_index] = func.evaluate( data )

  end_time = datetime.now()
  print("Program ended at", end_time, "after", end_time-start_time)
  
  ax[0,1].plot(range(len(values)), values, 'ro')

  if not os.path.exists('output'):  os.makedirs('output')
  plt.savefig('output/'+ file_name + '_parameter.png')
  if is_plot :
    plt.show()


def run_test_from_class(test_class):
  parameter_identification_test(test_class.nx, test_class.porosity, test_class.n_steps,
    test_class.jump_parameter, test_class.subset_sizes, test_class.n_choose_bins,
    test_class.min_value_shift, test_class.max_value_shift, test_class.jump_params, test_class.bins,
    test_class.distance_fct, test_class.debug_mode, test_class.file_name, test_class.is_plot)  



# --------------------------------------------------------------------------------------------------
# Define main function.
# -------------------------------------------------------------------------------------------------- 
if __name__ == "__main__":

  test_name    = 'basic_test'
  distances    = ["bulk_distance", "average_distance"]
  domain_sizes = [ 5, 10, 25, 50, 100] 
  sigmas       = [ 1,  5, 10, 25,  50]
  time_points  = [ 0,  5, 10, 25,  50]
  dimensions   = [ 1,  2,  3,  4,   5]

  try:
    import ecdf_test
  except (ImportError, ModuleNotFoundError) as error:
    sys.path.append(os.path.dirname(os.path.abspath(__file__)) + os.sep  + "parameters")
    import ecdf_test
  
  debug_mode = len(sys.argv) > 1 and sys.argv[1] == "True"
  fun_args   = []
  base_test = getattr(ecdf_test, test_name)


  for distance in distances:
    
    for size in domain_sizes:
      fun_args.append( base_test(
        nx           = [size, size],
        distance_fct = distance,
        file_name    = distance + 'domain_size_' + str(size)
        ) )

    for sigma in sigmas:
      fun_args.append( base_test(
        jump_parameter = sigma,
        distance_fct   = distance,
        file_name      = distance + 'sigma_' + str(sigma)
        ) )

    for steps in time_points:
      fun_args.append( base_test(
        n_steps      = steps,
        distance_fct = distance,
        file_name    = distance + 'time_steps_' + str(steps)
        ) )

    for dim in dimensions:
      fun_args.append( base_test(
        nx           = [ 50 for _ in range(dim) ],
        distance_fct = distance,
        file_name    = distance + 'dimension_' + str(dim)
        ) )

  processes = []
  for i in range(len(fun_args)):
    t = multiprocessing.Process(target=run_test_from_class, args=(fun_args[i],))
    processes.append(t)
    t.start()

  while multiprocessing.active_children():
    val = input("Enter your value: ")
    if val == "kill_all_children":
      active = multiprocessing.active_children()
      for child in active:
        child.kill()
      time.sleep(2)

  for one_process in processes:
    one_process.join()
