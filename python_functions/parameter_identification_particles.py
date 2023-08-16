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

try:
  from plot import plot, plot_to_file, plot_to_vtk
except (ImportError, ModuleNotFoundError) as error:
  sys.path.append(os.path.dirname(os.path.abspath(__file__)) + os.sep  + ".." + os.sep + 
      "python_functions")
  from plot import plot_to_file


illite_fine = [2, 6]
illite_medium = [2, 30]
goethite_fine = [1, 17]
goethite_coarse = [2, 34]
single_cell = [1,1]
uniform_positive = [1] * 6
uniform_negative = [-1] * 6

face_to_edge = [1,1, 1,1, -1, -1]
face_to_face = [-1,-1, -1, -1, -1,-1]

def split(x, n):
  ret = [0] * n
  #if(x < n):
     #print('error')
  #el
  if (x % n == 0):
      for i in range(n):
          ret[i] = x//n

  else:
      zp = n - (x % n)
      pp = x//n
      for i in range(n):
          if(i>= zp):
              ret[i] = pp + 1
          else:
              ret[i] = pp
  return ret
# --------------------------------------------------------------------------------------------------
# Parameter identification test.
# --------------------------------------------------------------------------------------------------
def ecdf_identify(nx, porosity, n_steps, jump_parameter, subaggregate_threshold, distribution, ecdf_type, subset_sizes, n_bins, n_runs,
  min_value_shift, max_value_shift, sigmas, distance_fct, debug_mode, file_name, is_plot):

  file = "output_mcmc_particles"
  if not os.path.exists(file):  os.makedirs(file)

  def stencil_size(jump_parameter, area):
    return (jump_parameter/(np.sqrt(area)))
  def run_cam(jump_parameter, subaggregate_threshold, distribution, nx, porosity, n_steps, debug_mode=False):
    CAM_wrapper = PyCAM(jump_parameter, subaggregate_threshold)
    n_fields = np.prod(nx)
    n_solids = round((1 - porosity) * n_fields)

    placement_i = n_solids * distribution
    placement_g = n_solids * (1- distribution)
    type_g = single_cell
    type_i = illite_fine

    charges_g = uniform_positive = [1] * 6
    charges_i = uniform_positive = [-1] * 6

    size_g = np.prod(type_g)
    size_i = np.prod(type_i)
    n_dim = len(nx)
    n_g = round(placement_g/size_g)
    n_i = round(placement_i/size_i)
    n_g_d = split(n_g, n_dim)
    n_i_d = split(n_i, n_dim)

    summe = 0

    #typ 1
    st_si_g = stencil_size(jump_parameter,size_g)
    for d in range(n_dim):
      success = 0
      type_g = list(np.roll(type_g,1))
      charges_g = list(np.roll(charges_g,2))
      while success <  n_g_d[d]:
        success = success + CAM_wrapper.place_plane(st_si_g,type_g,  -1, charges_g)
      summe = summe + success
    success = success + CAM_wrapper.place_plane(st_si_g,type_g,  -1, charges_g)
    summe = summe + success


     #type 2
    st_si_i = stencil_size(jump_parameter,size_i)
    for d in range(n_dim):
      success = 0
      type_i = list(np.roll(type_i,1))
      charges_i = list(np.roll(charges_i,2))
      while success < n_i_d[d]:
        success = success + CAM_wrapper.place_plane(st_si_i,type_i ,-1,charges_i)
      summe = summe + success
    success = success + CAM_wrapper.place_plane(st_si_i,type_i ,-1,charges_i)
    summe = summe + success
    # CAM_wrapper.place_single_cell_bu_randomly(jump_parameter, porosity , 0)


    if isinstance(n_steps, list):
      for _ in range(n_steps[0]):  CAM_wrapper.do_cam()
      data = [ CAM_wrapper.fields() ]
      for i in range(len(n_steps) - 1):
        for _ in range(n_steps[i], n_steps[i+1]): CAM_wrapper.do_cam()
        data.append( CAM_wrapper.fields() )
      return data
    else:
      for _ in range(n_steps):  CAM_wrapper.do_cam()
      return CAM_wrapper.fields()

  def generate_ecdf(data, subset_sizes, distance_fct, n_bins, ecdf_type, ax=None):
    min_val, max_val, distance_data = ecdf.estimate_radii_values(data[0:subset_sizes[0]],
      data[subset_sizes[0]:subset_sizes[0]+subset_sizes[1]], distance_fct )
    bins = np.linspace(min_val, max_val, 50)
    if ax is not None:
      if ecdf_type == "standard":
        aux_func = ecdf.standard(data, bins, distance_fct, subset_sizes)
      elif ecdf_type == "bootstrap":
        aux_func = ecdf.bootstrap(data, bins, distance_fct, subset_sizes[0], subset_sizes[1])
      ax = ecdf.plot_ecdf_vectors(aux_func, ax, 'm.')
    bins = ecdf.choose_bins(distance_data, bins, n_bins)
    if ecdf_type == "standard":
      func = ecdf.standard(data, bins, distance_fct, subset_sizes)
    elif ecdf_type == "bootstrap":
      func = ecdf.bootstrap(data, bins, distance_fct, subset_sizes[0], subset_sizes[1])
    return func, ax

  #for [jump_parameter, distribution] in sigmas:
    #print(jump_parameter, distribution)
  # print(len(sigmas))

  start_time = datetime.now()
  print("Starting time is", start_time)

  CAM_config            = CAM.config()
  CAM_config.nx         = nx
  CAM_config.debug_mode = debug_mode
  PyCAM                 = CAM.include(CAM_config)


  #-------plot----------
  # print("plot_start")
  # for  [jump_parameter, distribution] in sigmas:
  #   plot_data = run_cam(jump_parameter, subaggregate_threshold, distribution, nx, porosity, n_steps, debug_mode)
  #   plot_to_file(nx, plot_data, file + "/" + file_name + "_jp_" + str(jump_parameter) + "_distr_" + str(distribution) + '.png')
  # print("plot_end")
  #-------plot------------------

  n_fields, n_iter = np.prod(nx), np.sum(subset_sizes)
  data = [ run_cam(jump_parameter, subaggregate_threshold, distribution, nx, porosity, n_steps, debug_mode) for _ in range(n_iter) ]

  end_time = datetime.now()
  print("CAM data acquired at", end_time, "after", end_time-start_time)


  fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(18, 5))

  if not isinstance(n_steps, list) and (ecdf_type == "standard" or ecdf_type == "bootstrap"):
    func, ax[0,0] = generate_ecdf(data, subset_sizes, distance_fct, n_bins, ecdf_type, ax[0,0])
  elif isinstance(n_steps, list):
    if not all(n_steps[i] <= n_steps[i+1] for i in range(len(n_steps) - 1)):
      raise Exception("Number of stps must be ascending!")
    if not len(n_steps) == len(distance_fct) == len(ecdf_type) == len(n_bins) == len(n_steps):
      raise Exception("Illegal combination of parameter lists.")
    data = np.transpose(data, (1,0,2))
    func_list = []
    for i in range(len(distance_fct)):
      aux_func, _ = generate_ecdf(data[i], subset_sizes, distance_fct[i], n_bins[i], ecdf_type[i])
      func_list.append( aux_func )
    func = ecdf.multiple( func_list )
  else:
    print("ERROR: Type of ecdf is invalid.")

  ax[0,0] = ecdf.plot_ecdf_vectors(func, ax[0,0], 'c.')
  ax[0,0] = ecdf.plot_mean_vector(func, ax[0,0], 'k.')
  ax[1,0] = ecdf.plot_chi2_test(func, ax[1,0])

  end_time = datetime.now()
  print("Objective function setup at", end_time, "after", end_time-start_time)

  means_log = [0.] * len(sigmas)
  means_nor = [0.] * len(sigmas)

  for _ in range(n_runs):
    if not isinstance(n_steps, list):
      values = [ ecdf.evaluate( func, [ run_cam(jump_parameter, subaggregate_threshold, distribution, nx, porosity, n_steps, debug_mode) \
                 for _ in range(subset_sizes[0]) ] ) for [jump_parameter, distribution] in sigmas ]
    else:
      values = [ ecdf.evaluate( func, np.transpose( [ \
        run_cam(jump_parameter,subaggregate_threshold, distribution, nx, porosity, n_steps, debug_mode) for _ in range(subset_sizes[0]) \
        ], (1,0,2) ) ) for  [jump_parameter, distribution] in sigmas ]
    # ax[0,1].plot(sigmas, values, 'ro')
    ax[0,1].plot(range(len(sigmas)), values, 'ro')
    means_log = [ means_log[i] + values[i] / n_runs for i in range(len(sigmas))]
    values = [ np.exp(-0.5 * value)   for value in values ]
    values = [ value / np.sum(values) for value in values ]
    # ax[1,1].plot(sigmas, values, 'ro')
    ax[1,1].plot(range(len(sigmas)), values, 'ro')
    means_nor = [ means_nor[i] + values[i] / n_runs for i in range(len(sigmas))]

  labels = []
  for [jump_parameter, distribution] in sigmas:
    label = "[" + str(jump_parameter) +"," + str(distribution)+ "]"
    labels.append(label)

  #ax.set_xticks(x)
  #ax.set_xticklabels(xlabels, rotation=40, ha=ha[n])
  #ax[0,1].set(xticks=range(len(sigmas)), xticklabels=labels, rotation = 40)


  ax[0,1].plot(range(len(sigmas)), means_log, 'bo')
  ax[0,1].set_xticks(range(len(sigmas)))
  ax[0,1].set_xticklabels(labels, rotation=40, ha = 'right')

  ax[1,1].plot(range(len(sigmas)), means_nor, 'bo')
  ax[1,1].set_xticks(range(len(sigmas)))
  ax[1,1].set_xticklabels(labels, rotation=40, ha = 'right')

  plt.savefig(file + "/" + file_name + '.png')

  end_time = datetime.now()
  print("Program ended at", end_time, "after", end_time-start_time)

  if is_plot:  plt.show()


def run_test_from_class(test_class):
  ecdf_identify(test_class.nx, test_class.porosity, test_class.n_steps, test_class.jump_parameter, test_class.subaggregate_threshold, test_class.distribution, 
    test_class.ecdf_type, test_class.subset_sizes, test_class.n_bins, test_class.n_runs,
    test_class.min_value_shift, test_class.max_value_shift, test_class.sigmas,
    test_class.distance_fct, test_class.debug_mode, test_class.file_name, test_class.is_plot)
