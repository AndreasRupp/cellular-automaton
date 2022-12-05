from __future__ import print_function

import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
from multiprocessing import Pool
import os, sys
  

# --------------------------------------------------------------------------------------------------
# Parameter identification test.
# --------------------------------------------------------------------------------------------------
def parameter_identification_test(nx, porosity, n_steps, jump_parameter, n_iter, values, bins,
  n_choose_bins, subset_sizes, min_value_shift, max_value_shift, debug_mode, name, isPlotted = 0,):
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
  ax[0,0] = ecdf.plot_mean_vector(func, ax[0,0], 'k.')
  ax[1,0] = ecdf.plot_chi2_test(func, ax[1,0], 20)

  func.choose_bins(n_choose_bins, min_value_shift, max_value_shift)
  
  ax[0,0] = ecdf.plot_ecdf_vectors(func, ax[0,0], 'r.')
  ax[0,0] = ecdf.plot_mean_vector(func, ax[0,0], 'k.')
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
  plt.savefig('output/'+ name + '_parameter.png')
  if isPlotted :
    plt.show()


# --------------------------------------------------------------------------------------------------
# Define main function.
# -------------------------------------------------------------------------------------------------- 
if __name__ == "__main__":

    test_name = 'ecdf_parameter'

    try:
        import ecdf_test
    except (ImportError, ModuleNotFoundError) as error:
        sys.path.append(os.path.dirname(os.path.abspath(__file__)) + os.sep  + "parameters")
        import ecdf_test
    

    debug_mode = len(sys.argv) > 1 and sys.argv[1] == "True"


    domain_size = [5, 10 , 25, 50, 100] 
    sigma = [1,5,10,25,50]

    used_test = ecdf_test.ecdf_parameter()
    items = []
    for i in range(len(domain_size)):
      used_test.nx = [domain_size[i],domain_size[i]]
      name = 'domainSize' + str(domain_size[i])
      item = (used_test.nx, used_test.porosity, used_test.n_steps,
          used_test.jump_parameter, used_test.n_iter, used_test.values, used_test.bins,
          used_test.n_choose_bins, used_test.subset_sizes, used_test.min_value_shift,
          used_test.max_value_shift, debug_mode, name)
      items.append(item)

    with Pool() as pool:
      pool.starmap(parameter_identification_test,items)


    used_test = ecdf_test.ecdf_parameter()
    items = []
    for i in range(len(sigma)):
      used_test.jump_parameter = sigma[i]
      name = 'Sigma' + str(sigma[i])
      item = (used_test.nx, used_test.porosity, used_test.n_steps,
          used_test.jump_parameter, used_test.n_iter, used_test.values, used_test.bins,
          used_test.n_choose_bins, used_test.subset_sizes, used_test.min_value_shift,
          used_test.max_value_shift, debug_mode, name)
      items.append(item)

    with Pool() as pool:
      pool.starmap(parameter_identification_test,items)


    used_test = ecdf_test.ecdf_parameter()
    items = []  
    domain = []
    for i in range(5):
      domain.append(50)
      used_test.nx = domain
      name = 'Dim' + str(sigma[i])
      item =  (used_test.nx, used_test.porosity, used_test.n_steps,
          used_test.jump_parameter, used_test.n_iter, used_test.values, used_test.bins,
          used_test.n_choose_bins, used_test.subset_sizes, used_test.min_value_shift,
          used_test.max_value_shift, debug_mode, name)
      items.append(item)
      
    with Pool() as pool:
      pool.starmap(parameter_identification_test,items)