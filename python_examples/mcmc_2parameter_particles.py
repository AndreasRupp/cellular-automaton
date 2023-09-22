#  This tutorial was created with pymcmcstat package (https://pymcmcstat.readthedocs.io/en/latest/)

import numpy as np
from pymcmcstat import mcmcplot as mcp
from pymcmcstat.MCMC import MCMC
import matplotlib.pyplot as plt
import sys, os
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
  import ecdf_test
except (ImportError, ModuleNotFoundError) as error:
  sys.path.append(os.path.dirname(os.path.abspath(__file__)) + os.sep  + "parameters")
  import ecdf_test

illite_fine = [ 2, 6]
illite_medium = [ 2, 30]
goethite_fine = [1, 17]
goethite_coarse = [2,34]
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
  
  
  
nx             = [50, 50]
porosity       = 0.9
n_steps        = 10
jump_parameter = 5
subaggregate_threshold = 0.01
distribution           = 0.5

subset_sizes    = [100] * 40
min_value_shift = 0.1
max_value_shift = -0.1
n_choose_bins   = 20
jump_params     = "default"

debug_mode   = False
file_name    = "basic_test"
is_plot      = False

start_time = datetime.now()
print("Starting time is", start_time)

const            = CAM.config()
const.nx         = nx
const.ca_settings      = [False, False, False, False, False]
const.const_jump_parameter = jump_parameter
const.debug_mode      = debug_mode
PyCAM            = CAM.include(const)

CAM_wrapper = PyCAM()
distance_fct = CAM_wrapper.bulk_distance # <--- Needs to be removed

def stencil_size(jump_parameter, area, nx):
	return (jump_parameter/(area ** (1.0/len(nx))))
# --------------------------------------------------------------------------------------------------
# def run_cam(jump_parameter1, jump_parameter2, nx, porosity, n_steps, debug_mode=False):
def run_cam(jump_parameter, distribution, subaggregate_threshold, nx, porosity, n_steps, debug_mode=False):
	CAM_wrapper = PyCAM(jump_parameter, subaggregate_threshold)
	n_fields = np.prod(nx)
	n_solids = round((1 - porosity) * n_fields)

	placement_i = n_solids * distribution
	placement_g = n_solids * (1- distribution)
	type_g = goethite_coarse
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
	st_si_g = stencil_size(jump_parameter,size_g, nx)
	for d in range(n_dim):
		success = 0
		type_g = list(np.roll(type_g,1))
		charges_g = list(np.roll(charges_g,2))
		while success <  n_g_d[d]:
			success = success + CAM_wrapper.place_plane(st_si_g,type_g, -1, charges_g)
		summe = summe + success
	success = success + CAM_wrapper.place_plane(st_si_g,type_g,  -1, charges_g)
	summe = summe + success


	 #type 2
	st_si_i = stencil_size(jump_parameter,size_i, nx)
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
		for _ in range(n_steps[0]):
			CAM_wrapper.do_cam()
		data = [ CAM_wrapper.fields() ]
		for i in range(len(n_steps) - 1):
			for _ in range(n_steps[i], n_steps[i+1]):
				CAM_wrapper.do_cam()
			data.append( CAM_wrapper.fields() )
		return data
	else:
		for _ in range(n_steps):
			CAM_wrapper.do_cam()
		return CAM_wrapper.fields()
# --------------------------------------------------------------------------------------------------

n_fields = np.prod(nx)
n_iter   = np.sum(subset_sizes)
data = [ run_cam(jump_parameter, distribution, subaggregate_threshold, nx, porosity, n_steps, debug_mode) \
        for _ in range(n_iter) ]

end_time = datetime.now()
print("CAM data acquired at", end_time, "after", end_time-start_time)

min_val, max_val, distance_data = ecdf.estimate_radii_values(data[0:subset_sizes[0]],
  data[subset_sizes[0]:subset_sizes[0]+subset_sizes[1]], distance_fct )
bins = ecdf.choose_bins(distance_data, np.linspace(min_val, max_val, 50), n_choose_bins)
func = ecdf.standard(data, bins, distance_fct, subset_sizes)
  
end_time = datetime.now()
print("Objective function setup at", end_time, "after", end_time-start_time)


def evaluate_objective_function(theta, data=None):
  data = [ run_cam(theta[0], theta[1], subaggregate_threshold, nx, porosity, n_steps, debug_mode) for _ in range(subset_sizes[0]) ]
  return func.evaluate( data )


# Initialize the toolbox.
mcstat = MCMC()

# Add parameters for MCMC sampling with starting values.
mcstat.parameters.add_model_parameter(name='jump_parameter', theta0=0.0, minimum = 0, maximum = 10)
mcstat.parameters.add_model_parameter(name='distribution', theta0= 0.0, minimum = 0.0, maximum = 1.0)

# Passing to the toolbox our likelihood function
mcstat.model_settings.define_model_settings(sos_function=evaluate_objective_function)

# Defining simulation options
mcstat.simulation_options.define_simulation_options(
    nsimu=10000, #10000,  # number of elements in the chain
    qcov=np.eye(2),  # initial covariance matrix (if we do not know this beforehand, we start from identity matrix)
    adaptint=500  # adaptation interval of the method (helps MCMC algorithm to converge faster)
)

# Starting the simulation (the most time-consuming part)
mcstat.run_simulation()

# Extracting results
results = mcstat.simulation_results.results
chain = results['chain']  # parameter chains
# s2chain = results['s2chain']  # values of the likelihood for the respective parameter chains
names = results['names']  # parameter names

if not os.path.exists('output'):  os.makedirs('output')
# Plotting results
f1 = mcp.plot_chain_panel(chain, names=names)
plt.savefig('output/mcmc_2parameter_particles1.png')
f2 = mcp.plot_pairwise_correlation_panel(chain, names=names)
plt.savefig('output/mcmc_2parameter_particles2.png')
# plt.show()
#plt.savefig('output/mcmc_2parameter_particles.png')
