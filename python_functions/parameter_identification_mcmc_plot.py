from __future__ import print_function

import matplotlib.pyplot as plt
import numpy as np
import os, sys
from datetime import datetime
from pymcmcstat.MCMC import MCMC
from pymcmcstat import mcmcplot as mcp
#import xlsxwriter
import math 
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
  from plot import plot_to_file, plot_to_vtk

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

jump_parameter_factor = 1

# --------------------------------------------------------------------------------------------------
# Parameter identification test.
# --------------------------------------------------------------------------------------------------
def mcmc_identify(nx, porosity, n_steps, jump_parameter, subaggregate_threshold, distribution, type_i, type_g, charges_i, charges_g, ecdf_type, subset_sizes, n_bins,
  min_value_shift, max_value_shift, nsimu, qcov, adaptint, parameter_minmax, method,sigma2, N0, S20, distance_fct, debug_mode, filepath, file_name, is_plot):
  print(filepath, flush=True)

  if not os.path.exists(filepath):  os.makedirs(filepath)

  def diameter(PS):
    return 2 * math.pow(PS * (0.025 * 1000) ** 3 / math.pi * 3 / 4, (1/3))

  def stencil_size(jump_parameter, area, nx):
    return (jump_parameter/(area ** (1.0/len(nx))))
  def run_cam(jump_parameter, subaggregate_threshold, distribution, type_i, type_g, charges_i, charges_g, nx, porosity, n_steps, debug_mode=False):
    jump_parameter = jump_parameter * 1/jump_parameter_factor

    CAM_wrapper = PyCAM(jump_parameter, subaggregate_threshold)
    n_fields = np.prod(nx)
    n_solids = round((1 - porosity) * n_fields)

    placement_i = n_solids * distribution
    placement_g = n_solids * (1- distribution)

    size_g = np.prod(type_g)
    size_i = np.prod(type_i)
    n_dim = len(nx)

    n_g = round(placement_g/size_g)
    n_i = round((placement_i )/size_i)

    ##--------------------------------------------------
    n_fit_i = 0
    n_fit_g = 0
    _type_g = type_g
    _type_i = type_i
    for d in range(n_dim):
      fit_i = [i <= j for i,j in zip(_type_i,nx)]
      fit_g = [i <= j for i,j in zip(_type_g,nx)]
      if(not False in fit_i):
        n_fit_i = n_fit_i + 1
      if(not False in fit_g):
        n_fit_g = n_fit_g + 1
      _type_g = list(np.roll(_type_g,1))
      _type_i = list(np.roll(_type_i,1))

    #print(n_fit_g, n_fit_i)
    n_g_d = split(n_g, n_fit_g)
    n_i_d = split(n_i, n_fit_i) 

    _type_g = type_g
    _type_i = type_i
    #print(n_g_d, n_i_d)

    for d in range(n_dim):
      fit_i = [i <= j for i,j in zip(_type_i,nx)]
      fit_g = [i <= j for i,j in zip(_type_g,nx)]
      if(False in fit_i):
        n_i_d.insert(d,0)
      if(False in fit_g):
        n_g_d.insert(d,0)
      _type_g = list(np.roll(_type_g,1))
      _type_i = list(np.roll(_type_i,1))

    ##------------------------------------------------------------

    summe = 0

    #typ 1
    st_si_g = stencil_size(jump_parameter,size_g, nx)
    start_time = datetime.now()
    for d in range(n_dim):
      #print("dim ", d)
      success = 0
      while success <  n_g_d[d]:
        success = success + CAM_wrapper.place_plane(st_si_g,type_g,  -1, charges_g)
      summe = summe + success
      type_g = list(np.roll(type_g,1))
      charges_g = list(np.roll(charges_g,2))

     #type 2
    st_si_i = stencil_size(jump_parameter,size_i, nx)
    #print(jump_parameter, "->", st_si_i,st_si_g,"dist", distribution)
    for d in range(n_dim):
      #print("dim ",d)
      success = 0
      while success < n_i_d[d]:
        success = success + CAM_wrapper.place_plane(st_si_i,type_i ,-1,charges_i)
      summe = summe + success
      type_i = list(np.roll(type_i,1))
      charges_i = list(np.roll(charges_i,2))
   
    end_time = datetime.now()

    if isinstance(n_steps, list):
      for _ in range(n_steps[0]):  CAM_wrapper.do_cam()
      data = [ CAM_wrapper.particle_size_distribution_d() ]
      for i in range(len(n_steps) - 1):
        for _ in range(n_steps[i], n_steps[i+1]): CAM_wrapper.do_cam()
        data.append( CAM_wrapper.particle_size_distribution_d() )
      return data
    else:
      #plot_data = np.zeros( (n_steps + 1, np.prod(nx)) )
      #print("porosity ", sum([j > 0 for j in CAM_wrapper.fields()])/np.prod(nx))
      #data = CAM_wrapper.fields()
      #plot_data[0] = data
      for step in range(n_steps):  
        CAM_wrapper.do_cam()
        #data = CAM_wrapper.fields()
        #print("porosity ", sum([j > 0 for j in CAM_wrapper.fields()])/np.prod(nx))
        #plot_data[step + 1] = data
      #return 0
     
      numbers = CAM_wrapper.particle_size_distribution_d()
      #print(numbers)
      #print("väerndern")
      #for i in list(set(type_g +type_i)): 
        #numbers=[j for j in numbers if j != i] 
      #print(numbers)
      return numbers
      #return plot_data  


  def generate_ecdf(data, subset_sizes, distance_fct, n_bins, ecdf_type, ax=None):
    #print(subset_sizes)
    #print(subset_sizes[0])
    #print(len(data))
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


  start_time = datetime.now()
  print("Starting MCMC Script time is", start_time, flush=True)

  CAM_config            = CAM.config()
  CAM_config.nx         = nx
  CAM_config.debug_mode = debug_mode
  PyCAM                 = CAM.include(CAM_config)

  
  #-------plot----------
  #print("plot_start")
  jump_params = np.round(np.linspace(parameter_minmax[0][0],parameter_minmax[0][1],5, endpoint = True))
  distrib = np.round(np.linspace(parameter_minmax[1][0], parameter_minmax[1][1], 5, endpoint = True),1)
  jump_params = jump_parameter
  distrib = distribution

  test_start_time = datetime.now()
  plot_data = run_cam(jump_params, subaggregate_threshold, distrib,type_i, type_g, charges_i, charges_g, nx, porosity, n_steps, debug_mode)

  test_end_time = datetime.now()
  print("1 CAM nx ",str(nx), " types ",str(type_i)," ",str(type_g)," data duration", test_end_time, "after", test_end_time- test_start_time, flush=True)

  sigmas = []

  n_fields, n_iter = np.prod(nx), np.sum(subset_sizes)

  data = [] 
  file1 = open(filepath + '/CAMdata.txt', 'r')
  Lines = file1.readlines()
  count = 0
  for line in Lines:
    numbers = line.replace('\n', '').split(' ')
    numbers = [int(i) for i in numbers]
    #print(numbers)
    data.append(numbers)
    count = count + 1
    if(count >=n_iter):
       break
  #return
  #print(data)
  #print("----------------------------")
  print(len(data), n_iter)

  end_time = datetime.now()
  print("CAM data acquired at", end_time, "after", end_time-start_time, flush=True)

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
  print("Objective function setup at", end_time, "after", end_time-start_time, flush=True)

  
  def evaluate_objective_function(theta, data=None):
      #print("t0" + str(theta[0]) +"t1" +str(theta[1]))
      theta1 = 0.25
      theta0 = 5
      if not isinstance(n_steps, list):
          data = [ run_cam(theta[0], subaggregate_threshold, theta[1], type_i, type_g, charges_i, charges_g, nx, porosity, n_steps, debug_mode) for _ in range(subset_sizes[0]) ]
      else:
          data = np.transpose([ run_cam(theta[0], subaggregate_threshold, theta[1], type_i, type_g, charges_i, charges_g, nx, porosity, n_steps, debug_mode) for _ in range(subset_sizes[0]) ], (1,0,2))
      return func.evaluate( data )
  print("start sos")
  jump_params = np.round(np.linspace(0,20,20, endpoint = True)) 
  distrib = np.round(np.linspace(0,1, 101, endpoint = True),2)
  sigmas = []
  for jp in jump_params:
    for dis in distrib:
        sigma = [jp,dis]
        print(sigma)
        with open('plotSOS.txt' , 'a') as fp:
           try:
              fp.write('%.3f %.3f %d \n' % (jp, dis, evaluate_objective_function([jp,dis]) ))#
           except:
              fp.write('%.3f %.3f %d \n' % (jp, dis, 0 ))	#


               
            
  print("Program ended at", end_time, "after", end_time-start_time, flush=True)



def run_test_from_class(test_class):
  #print("testclass jp", test_class.jump_parameter, " dst ", test_class.distribution)
  mcmc_identify(test_class.nx, test_class.porosity, test_class.n_steps, test_class.jump_parameter, test_class.subaggregate_threshold, test_class.distribution, 
    test_class.type_i, test_class.type_g, test_class.charges_i, test_class.charges_g, 
    test_class.ecdf_type, test_class.subset_sizes, test_class.n_bins,
    test_class.min_value_shift, test_class.max_value_shift,
    test_class.nsimu, test_class.qcov, test_class.adaptint, test_class.parameter_minmax, test_class.method,test_class.sigma2, test_class.N0, test_class.S20,
    test_class.distance_fct, test_class.debug_mode, test_class.filepath, test_class.file_name, test_class.is_plot)
