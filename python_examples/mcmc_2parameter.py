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


nx             = [50, 50]
porosity       = 0.7
n_steps        = 5
jump_parameter = 5

subset_sizes    = [100] * 2
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
const.debug_mode = debug_mode
PyCAM            = CAM.include(const)

CAM_wrapper = PyCAM()
distance_fct = CAM_wrapper.bulk_distance # <--- Needs to be removed

# --------------------------------------------------------------------------------------------------
def run_cam(jump_parameter1, jump_parameter2, nx, porosity, n_steps, debug_mode=False):
  CAM_wrapper = PyCAM(jump_parameter2)
  CAM_wrapper.place_single_cell_bu_randomly(jump_parameter1, porosity,  0)
  for step in range(n_steps): 
    CAM_wrapper.do_cam()
  return CAM_wrapper.fields()
# --------------------------------------------------------------------------------------------------

n_fields = np.prod(nx)
n_iter   = np.sum(subset_sizes)
data = [ run_cam(jump_parameter, jump_parameter, nx, porosity, n_steps, debug_mode) \
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
  data = [ run_cam(theta[0], theta[1], nx, porosity, n_steps, debug_mode) for _ in range(n_iter) ]
  return func.evaluate( data )


# Initialize the toolbox.
mcstat = MCMC()

# Add parameters for MCMC sampling with starting values.
mcstat.parameters.add_model_parameter(name='theta_1', theta0=5.0)
mcstat.parameters.add_model_parameter(name='theta_2', theta0=5.0)

# Passing to the toolbox our likelihood function
mcstat.model_settings.define_model_settings(sos_function=evaluate_objective_function)

# Defining simulation options
mcstat.simulation_options.define_simulation_options(
    nsimu=10000,  # number of elements in the chain
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

# Plotting results
f1 = mcp.plot_chain_panel(chain, names=names)
f2 = mcp.plot_pairwise_correlation_panel(chain, names=names)
# plt.show()
if not os.path.exists('output'):  os.makedirs('output')
plt.savefig('output/mcmc_2parameter.png')
