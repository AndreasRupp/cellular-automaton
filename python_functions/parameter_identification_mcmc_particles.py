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
  #print("identifyjp ", jump_parameter, " dist ", distribution)
  if not os.path.exists(filepath):  os.makedirs(filepath)
  #workbook = xlsxwriter.Workbook(filepath + "/" + "measure_jp_" + str(jump_parameter) + "_dist_"+ str(distribution)+ '.xlsx', {'nan_inf_to_errors': True})

  #worksheet = workbook.add_worksheet()

  #eval_measures_names = ['n_solids',
        #  'n_particles',
        # 'n_composites',
         # 'n_solids_in_composites',
         # 'n_solids_in_discrete_BUs',
         # 'mean_particle_size',
         # 'specific_surface_area',
         # 'contact_area_per_volume',
         # 'mean_diameter']
  #worksheet.write(0, 0, 'time_step')
  #column = 1
  #for name in eval_measures_names:
   # worksheet.write(0,column, name)
    #column += 1
  def diameter(PS):
    return 2 * math.pow(PS * (0.025 * 1000) ** 3 / math.pi * 3 / 4, (1/3))

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

    #n_g = round(placement_g/size_g)
    #rest = placement_g - n_g * size_g
    #print(rest)
    #n_i = round((placement_i )/size_i)
    #n_g_d = split(n_g, n_dim)
    #n_i_d = split(n_i, n_dim)
    
    #print("distributon of g and i type ",n_g_d, n_i_d)
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
    #success = success + CAM_wrapper.place_plane(st_si_g,type_g,  -1, charges_g)
    #summe = summe + success


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
    #print("porosity ", sum([j > 0 for j in CAM_wrapper.fields()])/np.prod(nx))
    #success = success + CAM_wrapper.place_plane(st_si_i,type_i ,-1,charges_i)
    #summe = summe + success
    # CAM_wrapper.place_single_cell_bu_randomly(jump_parameter, porosity , 0)
    end_time = datetime.now()
    #print("setup CAM", end_time, "after", end_time- test_start_time)
    # particle_size_distribution_d
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





      #data[(data <= n_g) & (data > 0)] = 1
      #data[data > n_g] = 2
      #plot_data[0] = data  
      
      #return CAM_wrapper.particle_size_distribution_d()
      #return CAM_wrapper.fields()
        #data = CAM_wrapper.fields()
        #eval_data = CAM_wrapper.eval_measures()
        #eval_data[len(eval_measures_names)-1] = diameter(eval_data[5])
        #worksheet.write(step + 2, 0, step+1)
        #for i in range(len(eval_measures_names)):
          #worksheet.write(step + 2, i+1, eval_data[i])
        #data[(data <= n_g) & (data > 0)] = 1
        #data[data > n_g] = 2
        #plot_data[step + 1] = data
      #d1 = CAM_wrapper.particle_size_distribution_d()
      #d2 = CAM_wrapper.particle_size_distribution(CAM_wrapper.fields())
      #print("##################################")
      #print(d1)
      #print('-----')
      #print(d2)
      #print("##################################")
      

    
      #plot_data[(plot_data <= n_g) & (plot_data > 0)] = 1
      #plot_data[plot_data > n_g] = 2
      #workbook.close()
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

  #for [jump_parameter, distribution] in sigmas:
    #print(jump_parameter, distribution)
  # print(len(sigmas))

  start_time = datetime.now()
  print("Starting MCMC Script time is", start_time, flush=True)

  CAM_config            = CAM.config()
  CAM_config.nx         = nx
  CAM_config.debug_mode = debug_mode
  PyCAM                 = CAM.include(CAM_config)

  #parameter_minmax[0][0] = parameter_minmax[0][0] * jump_parameter_factor
  #parameter_minmax[0][1] = parameter_minmax[0][1] * jump_parameter_factor
  
  #-------plot----------
  #print("plot_start")
  jump_params = np.round(np.linspace(parameter_minmax[0][0],parameter_minmax[0][1],5, endpoint = True))
  distrib = np.round(np.linspace(parameter_minmax[1][0], parameter_minmax[1][1], 5, endpoint = True),1)
  jump_params = jump_parameter
  distrib = distribution

  test_start_time = datetime.now()
  plot_data = run_cam(jump_params, subaggregate_threshold, distrib,type_i, type_g, charges_i, charges_g, nx, porosity, n_steps, debug_mode)
  #print(plot_data)
  # print(len(plot_data))
  #plot_to_file(nx, plot_data[0], filepath + "/" + "domain")
  #plot_to_vtk( filepath + "/" + "domain",plot_data, nx)
  test_end_time = datetime.now()
  print("1 CAM nx ",str(nx), " types ",str(type_i)," ",str(type_g)," data duration", test_end_time, "after", test_end_time- test_start_time, flush=True)
  #plot_to_file(nx, plot_data, filepath + "/" + "domain")
  #print(distrib)
  #print(jump_params)
  sigmas = []
  # for jp in jump_params:
	 #  for dis in distrib:
		#   sigmas.append([jp,dis])
  #print(sigmas)  
  #for nr in range(2):
    #for  [jump_parameter_v, distribution_v] in sigmas:
      #plot_data = run_cam(jump_parameter_v, subaggregate_threshold, distribution_v,type_i, type_g, charges_i, charges_g, nx, porosity, n_steps, debug_mode)
      #print(plot_data)
      #print(jump_parameter_v, distribution_v)
      #plot_to_file(nx, plot_data, filepath + "/" + "domain_[" + str(jump_parameter_v) + "," + str(distribution_v) + ']' + str(nr)+'.png')
  #return
  #print("plot_end")
  #-------plot------------------
  #print("gallo")
  n_fields, n_iter = np.prod(nx), np.sum(subset_sizes)
  #print("datajp ", jump_parameter, " dist ", distribution)
  #data = [ run_cam(jump_parameter * jump_parameter_factor, subaggregate_threshold, distribution,type_i, type_g, charges_i, charges_g, nx, porosity, n_steps, debug_mode) for _ in range(n_iter) ]
  #with open(filepath + '/CAMdata.txt', 'a') as fp:
    #for item in data:
      #print(item)
      #writeSpace = False
      #for number in item:
        #print(number)
        #if(writeSpace):
          #fp.write(" ")
        #writeSpace=True
        #fp.write("%s" % number)
      #fp.write("\n")
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

    #mcstat.data.add_data_set(x , y)
  # Add parameters for MCMC sampling with starting values., theta0=parameter_minmax[0][0],theta0= parameter_minmax[1][0],+ parameter_minmax[0][1])/2  + parameter_minmax[1][1])/2
  #print(parameter_minmax[0][0])#theta0= (parameter_minmax[1][0] +parameter_minmax[1][1])/2,
  #print(parameter_minmax[0][1])#theta0=(parameter_minmax[0][0] +parameter_minmax[0][1])/2,
  print("theta0 ", parameter_minmax[1][0],parameter_minmax[0][0], flush=True)
  mcstat.parameters.add_model_parameter(name='jump_parameter',theta0=parameter_minmax[0][0], minimum = 0)#, , prior_mu = 20
  mcstat.parameters.add_model_parameter(name='distribution',theta0 = parameter_minmax[1][0], minimum = 0, maximum = 1.0)#

  # Passing to the toolbox our likelihood function
  mcstat.model_settings.define_model_settings(sos_function=evaluate_objective_function)#, sigma2= sigma2 , N0 = N0, S20 =S20

  #file1 = open(filepath + '/chainfile.txt', 'r')
  #Lines = file1.readlines()
  #print("len ", len(Lines))

  # Defining simulation options
  mcstat.simulation_options.define_simulation_options(
    nsimu= 3000,  # number of elements in the chain
    #method = method,
    qcov=np.eye(2),  # initial covariance matrix (if we do not know this beforehand, we start from identity matrix)
    adaptint=500,  # adaptation interval of the method (helps MCMC algorithm to converge faster)
    savedir = filepath,
    save_to_txt = True,
    savesize  =1,
    #updatesigma=1,
    #method = method,
  )
  mcstat.model_settings.display_model_settings(['sos_function', 'model_function', 'sigma2', 'N', 'N0', 'S20', 'nbatch'])

  # Starting the simulation (the most time-consuming part)
  print("start mcmc", flush= True)
  mcstat.run_simulation()

  # while(True):
    # file1 = open(filepath + '/chainfile.txt', 'r')
    # Lines = file1.readlines()
    # if(len(Lines) >= 4000):
       # return
    # # Defining simulation options
    # mcstat1 = MCMC()
    # mcstat1.data.add_data_set(x , y)
    # mcstat1.parameters.add_model_parameter(name='jump_parameter' , sample=True, minimum = parameter_minmax[0][0], maximum = parameter_minmax[0][1])#, , prior_mu = 20
    # mcstat1.parameters.add_model_parameter(name='distribution',sample=True, minimum = parameter_minmax[1][0], maximum = parameter_minmax[1][1])#
    # mcstat1.model_settings.define_model_settings(sos_function=evaluate_objective_function,sigma2= sigma2 , N0 = N0, S20 =S20)
    # mcstat1.simulation_options.define_simulation_options(
      # nsimu=nsimu, #10000,  # number of elements in the chain
      # savesize = 1,
      # savedir = filepath,
      # save_to_txt = True,
      # json_restart_file=filepath +'/results' ,
      # results_filename='results',
      # save_to_json=True,
      # save_lightly=True,
      # method = method,
      # #qcov=qcov,  # initial covariance matrix (if we do not know this beforehand, we start from identity matrix)
	  # adaptint=adaptint  # adaptation interval of the method (helps MCMC algorithm to converge faster)
    # ) 
    # #mcstat1.model_settings.display_model_settings([ 'model_function', 'sigma2', 'N', 'N0', 'S20', 'nbatch'])

    # # Starting the simulation (the most time-consuming part)
    # print("start mcmc")
    # mcstat1.run_simulation()
    # #print("len ", len(Lines))


  mcstat.model_settings.display_model_settings(['sos_function', 'model_function', 'sigma2', 'N', 'N0', 'S20', 'nbatch'])
  # Extracting results
  results = mcstat.simulation_results.results
  chain = results['chain']  # parameter chains
  # s2chain = results['s2chain']  # values of the likelihood for the respective parameter chains
  names = results['names']  # parameter names
  
  burnin = int(chain.shape[0] * 9/10) #100#int(chain.shape[0]/2)
  # display chain statistics
  stat = mcstat.chainstats(chain[burnin:, :], results, True)
  #mean0 = stat['mean'][0]
  #mean1 = stat['mean'][1]
  #std0 = stat['std'][0]
  #std1 = stat['std'][1]
  # Plotting results
  f1 = mcp.plot_chain_panel(chain, names=names)
  plt.savefig(filepath + "/" + file_name + '_chain.png')
  f2 = mcp.plot_chain_panel(chain[burnin:,:], names=names)
  plt.savefig(filepath + "/" + file_name + '_chain_burned.png')#_mean[' +str(round(mean0,2)) +','+str(round(mean1,2))+']std[' + str(round(std0,2)) + ',' + str(round(std1,2)) + ']
  f3 = mcp.plot_density_panel(chain[burnin:,:], names)
  plt.savefig(filepath + "/" + file_name + '_densitiy.png')
  #f3 = mcp.plot_chain_panel(s2chain, names=names)
  #plt.savefig(filepath + "/" + file_name + '_s2chain.png')
  f4 = mcp.plot_pairwise_correlation_panel(chain[burnin:,:], names=names)
  plt.savefig(filepath + "/" + file_name + '_pairwise_correlation.png')
  
  print("Program ended at", end_time, "after", end_time-start_time, flush=True)

  if is_plot:  plt.show()


def run_test_from_class(test_class):
  #print("testclass jp", test_class.jump_parameter, " dst ", test_class.distribution)
  mcmc_identify(test_class.nx, test_class.porosity, test_class.n_steps, test_class.jump_parameter, test_class.subaggregate_threshold, test_class.distribution, 
    test_class.type_i, test_class.type_g, test_class.charges_i, test_class.charges_g, 
    test_class.ecdf_type, test_class.subset_sizes, test_class.n_bins,
    test_class.min_value_shift, test_class.max_value_shift,
    test_class.nsimu, test_class.qcov, test_class.adaptint, test_class.parameter_minmax, test_class.method,test_class.sigma2, test_class.N0, test_class.S20,
    test_class.distance_fct, test_class.debug_mode, test_class.filepath, test_class.file_name, test_class.is_plot)
