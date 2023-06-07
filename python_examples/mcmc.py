#  This tutorial was created with pymcmcstat package (https://pymcmcstat.readthedocs.io/en/latest/)

import numpy as np
from pymcmcstat import mcmcplot as mcp
from pymcmcstat.MCMC import MCMC
import matplotlib.pyplot as plt

import sys, os
from datetime import datetime




nx             = [50, 50]
porosity       = 0.7
n_steps        = 5
jump_parameter = 5

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

try:
  import ecdf_test
except (ImportError, ModuleNotFoundError) as error:
  sys.path.append(os.path.dirname(os.path.abspath(__file__)) + os.sep  + "parameters")
  import ecdf_test


const            = CAM.config()
const.nx         = nx
const.debug_mode = debug_mode
PyCAM            = CAM.include(const)


distance_fct = PyCAM.bulk_distance # <--- Needs to be removed

# --------------------------------------------------------------------------------------------------
def run_cam(jump_parameter1, jump_parameter2, nx, porosity, n_steps, debug_mode=False):
  CAM_wrapper = PyCAM(jump_parameter2)
  CAM_wrapper.place_singleCellBU_randomlyy(porosity, jump_parameter1, 0)
  for step in range(n_steps):  CAM_wrapper.do_CAM()
  return CAM_wrapper.fields()
# --------------------------------------------------------------------------------------------------

n_fields = np.prod(nx)
n_iter   = np.sum(subset_sizes)
data = [[0] * n_fields] * n_iter


for iter in range(n_iter):
  data[iter] = run_cam(jump_parameter, nx, porosity, n_steps, debug_mode)

end_time = datetime.now()
print("CAM data acquired at", end_time, "after", end_time-start_time)

min_val, max_val, _ = ecdf.estimate_radii_values(
data[0:subset_sizes[0]], data[subset_sizes[0]:subset_sizes[0]+subset_sizes[1]], distance_fct )
bins = np.linspace(min_val, max_val, 50)

func = ecdf.estimator(data, bins, distance_fct, subset_sizes)

# fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(18, 5))
# ax[0,0] = ecdf.plot_ecdf_vectors(func, ax[0,0])
# ax[0,0] = ecdf.plot_mean_vector(func, ax[0,0], 'k.')
# ax[1,0] = ecdf.plot_chi2_test(func, ax[1,0], 20)

func.choose_bins(n_choose_bins, min_value_shift, max_value_shift)
  
# ax[0,0] = ecdf.plot_ecdf_vectors(func, ax[0,0], 'r.')
# ax[0,0] = ecdf.plot_mean_vector(func, ax[0,0], 'k.')
# ax[1,1] = ecdf.plot_chi2_test(func, ax[1,1])
  
end_time = datetime.now()
print("Objective function setup at", end_time, "after", end_time-start_time)





# Let us assume that we have only two parameters in our model. Suppose that from the experiments we have estimated
# the mean vector and covariance matrix of our likelihood. In this simple example we just prescribe them manually.
# mu = np.array([[1], [1]])
# Sigma = np.array([[5, 0], [0, 1]])

# pymcmcstat toolbox requires the user to define a function, which evaluates the likelihood at given parameter vector
# theta. If some additional information is needed (which does not depend on parameters), it can be supplied as a second
# argument data.
# In practice, this simple function should be replaced by the call to our 'objective_function' object.
def evaluate_objective_function(theta, data=None):
  data   = [[0] * n_fields] * subset_sizes[0]
  for iter in range(subset_sizes[0]):
    data[iter] = run_cam(theta, nx, porosity, n_steps, debug_mode)
  return func.evaluate( data )

  # theta_vec = np.array(theta)
  # theta_vec = np.reshape(theta, (2, 1))
  # mu_centered = theta_vec - mu
  # result = np.matmul(np.transpose(mu_centered), np.linalg.solve(Sigma, mu_centered))
  # result = float(result)
  # return result


# Initialize the toolbox.
mcstat = MCMC()

# Add parameters for MCMC sampling with starting values.
mcstat.parameters.add_model_parameter(name='theta_1', theta0=5.0)
# mcstat.parameters.add_model_parameter(name='theta_2', theta0=0.0)

# Passing to the toolbox our likelihood function
mcstat.model_settings.define_model_settings(sos_function=evaluate_objective_function)

# Defining simulation options
mcstat.simulation_options.define_simulation_options(
    nsimu=3000,  # number of elements in the chain
    qcov=np.eye(1),  # initial covariance matrix (if we do not know this beforehand, we start from identity matrix)
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
# f2 = mcp.plot_pairwise_correlation_panel(chain, names=names)
# plt.show()
if not os.path.exists('output'):  os.makedirs('output')
plt.savefig('output/mcmc.png')
