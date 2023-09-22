from __future__ import print_function

import matplotlib.pyplot as plt
import numpy as np
import os, sys
from datetime import datetime
from pymcmcstat.MCMC import MCMC
from pymcmcstat import mcmcplot as mcp

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
  print(filepath)
  #print("identifyjp ", jump_parameter, " dist ", distribution)
  if not os.path.exists(filepath):  os.makedirs(filepath)

  def stencil_size(jump_parameter, area, nx):
    return (jump_parameter/(area ** (1.0/len(nx))))
  def run_cam(jump_parameter, subaggregate_threshold, distribution, type_i, type_g, charges_i, charges_g, nx, porosity, n_steps, debug_mode=False):
    jump_parameter = jump_parameter * 1/jump_parameter_factor
    #print("jp ", jump_parameter, " dist ", distribution)
    CAM_wrapper = PyCAM(jump_parameter, subaggregate_threshold)
    n_fields = np.prod(nx)
    n_solids = round((1 - porosity) * n_fields)

    placement_i = n_solids * distribution
    placement_g = n_solids * (1- distribution)

    size_g = np.prod(type_g)
    size_i = np.prod(type_i)
    
    n_dim = len(nx)
    n_g = round(placement_g/size_g)
    #rest = placement_g - n_g * size_g
    #print(rest)
    n_i = round((placement_i )/size_i)
    n_g_d = split(n_g, n_dim)
    n_i_d = split(n_i, n_dim)

    summe = 0

    #typ 1
    st_si_g = stencil_size(jump_parameter,size_g, nx)
    for d in range(n_dim):
      success = 0
      type_g = list(np.roll(type_g,1))
      charges_g = list(np.roll(charges_g,2))
      while success <  n_g_d[d]:
        success = success + CAM_wrapper.place_plane(st_si_g,type_g,  -1, charges_g)
      summe = summe + success
    #success = success + CAM_wrapper.place_plane(st_si_g,type_g,  -1, charges_g)
    #summe = summe + success


     #type 2
    st_si_i = stencil_size(jump_parameter,size_i, nx)
    #print(jump_parameter, "->", st_si_i,st_si_g,"dist", distribution)
    for d in range(n_dim):
      success = 0
      type_i = list(np.roll(type_i,1))
      charges_i = list(np.roll(charges_i,2))
      while success < n_i_d[d]:
        success = success + CAM_wrapper.place_plane(st_si_i,type_i ,-1,charges_i)
      summe = summe + success
    #success = success + CAM_wrapper.place_plane(st_si_i,type_i ,-1,charges_i)
    #summe = summe + success
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

  #for [jump_parameter, distribution] in sigmas:
    #print(jump_parameter, distribution)
  # print(len(sigmas))

  start_time = datetime.now()
  print("Starting MCMC Script time is", start_time)

  CAM_config            = CAM.config()
  CAM_config.nx         = nx
  CAM_config.debug_mode = debug_mode
  PyCAM                 = CAM.include(CAM_config)

  parameter_minmax[0][0] = parameter_minmax[0][0] * jump_parameter_factor
  parameter_minmax[0][1] = parameter_minmax[0][1] * jump_parameter_factor
  
  #-------plot----------
  print("plot_start")
  jump_params = np.round(np.linspace(parameter_minmax[0][0],parameter_minmax[0][1],5, endpoint = True))
  distrib = np.round(np.linspace(parameter_minmax[1][0], parameter_minmax[1][1], 5, endpoint = True),1)
  #print(distrib)
  #print(jump_params)
  sigmas = []
  for jp in jump_params:
	  for dis in distrib:
		  sigmas.append([jp,dis])
  #print(sigmas)  
  #for  [jump_parameter_v, distribution_v] in sigmas:
    #plot_data = run_cam(jump_parameter_v, subaggregate_threshold, distribution_v,type_i, type_g, charges_i, charges_g, nx, porosity, n_steps, debug_mode)
    #print(jump_parameter_v, distribution_v)
    #plot_to_file(nx, plot_data, filepath + "/" + "domain_[" + str(jump_parameter_v) + "," + str(distribution_v) + '].png')
  #return
  print("plot_end")
  #-------plot------------------
  #print("gallo")
  n_fields, n_iter = np.prod(nx), np.sum(subset_sizes)
  #print("datajp ", jump_parameter, " dist ", distribution)
  data = [ run_cam(jump_parameter * jump_parameter_factor, subaggregate_threshold, distribution,type_i, type_g, charges_i, charges_g, nx, porosity, n_steps, debug_mode) for _ in range(n_iter) ]
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

  
  def evaluate_objective_function(theta, data=None):
	  #print("t0" + str(theta[0]) +"t1" +str(theta[1]))
	  if not isinstance(n_steps, list):
		  data = [ run_cam(theta[0], subaggregate_threshold, theta[1], type_i, type_g, charges_i, charges_g, nx, porosity, n_steps, debug_mode) for _ in range(subset_sizes[0]) ]
	  else:
		  data = np.transpose([ run_cam(theta[0], subaggregate_threshold, theta[1], type_i, type_g, charges_i, charges_g, nx, porosity, n_steps, debug_mode) for _ in range(subset_sizes[0]) ], (1,0,2))
	  return func.evaluate( data )
  
  # Initialize the toolbox.
  mcstat = MCMC()

  for i in range(len(subset_sizes)):    
    #print(np.sum(subset_sizes[0:i]),np.sum(subset_sizes[0:i+1]))
    #print(range(subset_sizes[i]))
    # data structure 
    x = range(subset_sizes[i])
    #print(x)
    y = data[int(np.sum(subset_sizes[0:i])):int(np.sum(subset_sizes[0:i+1]))]
    #print(len(y))
    mcstat.data.add_data_set(x , y)

  # Add parameters for MCMC sampling with starting values., theta0=parameter_minmax[0][0],theta0= parameter_minmax[1][0],+ parameter_minmax[0][1])/2  + parameter_minmax[1][1])/2
  #print(parameter_minmax[0][0])theta0= (parameter_minmax[1][0] +parameter_minmax[1][1])/2,
  #print(parameter_minmax[0][1])theta0=(parameter_minmax[0][0] +parameter_minmax[0][1])/2,
  mcstat.parameters.add_model_parameter(name='jump_parameter' , sample=True, minimum = parameter_minmax[0][0])#, maximum = parameter_minmax[0][1], prior_mu = 20
  mcstat.parameters.add_model_parameter(name='distribution',sample=True, minimum = parameter_minmax[1][0], maximum = parameter_minmax[1][1])#

  # Passing to the toolbox our likelihood function
  mcstat.model_settings.define_model_settings(sos_function=evaluate_objective_function,sigma2= sigma2 , N0 = N0, S20 =S20)

  # Defining simulation options
  mcstat.simulation_options.define_simulation_options(
    nsimu=nsimu, #10000,  # number of elements in the chain
    savesize = 1,
    savedir = filepath,
    save_to_txt = True,
    updatesigma=1,
    #burnintime = 50,
    method = method,#'mh',am
    #alphatarget = 0.2,
	  qcov=qcov,  # initial covariance matrix (if we do not know this beforehand, we start from identity matrix)
	  adaptint=adaptint  # adaptation interval of the method (helps MCMC algorithm to converge faster)
  ) 
  mcstat.model_settings.display_model_settings(['sos_function', 'model_function', 'sigma2', 'N', 'N0', 'S20', 'nbatch'])

  # Starting the simulation (the most time-consuming part)
  mcstat.run_simulation()
  mcstat.model_settings.display_model_settings(['sos_function', 'model_function', 'sigma2', 'N', 'N0', 'S20', 'nbatch'])
  # Extracting results
  results = mcstat.simulation_results.results
  chain = results['chain']  # parameter chains
  # s2chain = results['s2chain']  # values of the likelihood for the respective parameter chains
  names = results['names']  # parameter names
  
  burnin = int(chain.shape[0]/2)
  # display chain statistics
  stat = mcstat.chainstats(chain[burnin:, :], results, True)
  mean0 = stat['mean'][0]
  mean1 = stat['mean'][1]
  std0 = stat['std'][0]
  std1 = stat['std'][1]
  # Plotting results
  f1 = mcp.plot_chain_panel(chain, names=names)
  plt.savefig(filepath + "/" + file_name + '_chain.png')
  f2 = mcp.plot_chain_panel(chain[burnin:,:], names=names)
  plt.savefig(filepath + "/" + file_name + '_chain_burned_mean[' +str(round(mean0,2)) +','+str(round(mean1,2))+']std[' + str(round(std0,2)) + ',' + str(round(std1,2)) + '].png')
  f3 = mcp.plot_density_panel(chain[burnin:,:], names)
  plt.savefig(filepath + "/" + file_name + '_densitiy.png')
  #f3 = mcp.plot_chain_panel(s2chain, names=names)
  #plt.savefig(filepath + "/" + file_name + '_s2chain.png')
  f4 = mcp.plot_pairwise_correlation_panel(chain[burnin:,:], names=names)
  plt.savefig(filepath + "/" + file_name + '_pairwise_correlation.png')
  
  print("Program ended at", end_time, "after", end_time-start_time)

  if is_plot:  plt.show()


def run_test_from_class(test_class):
  #print("testclass jp", test_class.jump_parameter, " dst ", test_class.distribution)
  mcmc_identify(test_class.nx, test_class.porosity, test_class.n_steps, test_class.jump_parameter, test_class.subaggregate_threshold, test_class.distribution, 
    test_class.type_i, test_class.type_g, test_class.charges_i, test_class.charges_g, 
    test_class.ecdf_type, test_class.subset_sizes, test_class.n_bins,
    test_class.min_value_shift, test_class.max_value_shift,
    test_class.nsimu, test_class.qcov, test_class.adaptint, test_class.parameter_minmax, test_class.method,test_class.sigma2, test_class.N0, test_class.S20,
    test_class.distance_fct, test_class.debug_mode, test_class.filepath, test_class.file_name, test_class.is_plot)
