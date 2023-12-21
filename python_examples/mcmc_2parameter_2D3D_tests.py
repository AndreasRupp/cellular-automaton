#!/bin/env python

from __future__ import print_function
from datetime import datetime
import multiprocessing
from multiprocessing import cpu_count
#from pathos.pp import ParallelPool as Pool
from pathos.multiprocessing import ProcessingPool as Pool
import os, sys, time
import numpy as np
sys.path.append(os.getcwd())

try:
  from parameter_identification_mcmc_particles import run_test_from_class
except (ImportError, ModuleNotFoundError) as error:
  sys.path.append(os.path.dirname(os.path.abspath(__file__)) + os.sep  + ".." + os.sep +
    "python_functions")
  from parameter_identification_mcmc_particles import run_test_from_class

try:
  import mcmc_test
except (ImportError, ModuleNotFoundError) as error:
  sys.path.append(os.path.dirname(os.path.abspath(__file__)) + os.sep  + "parameters")
  import mcmc_test


illite_fine = [ 2, 6]
illite_medium = [ 2, 30]
goethite_fine = [2, 17]
goethite_coarse = [2,34]

illite_fine_3D = [ 1, 1,3]
illite_medium_3D= [ 1, 1,15]
goethite_fine_3D= [1,1, 8]
goethite_coarse_3D = [1,1,17]

uniform_positive_3D = [1] * 6 
uniform_negative_3D = [-1] * 6 
uniform_positive = [1] * 4
uniform_negative = [-1] * 4 

face_to_edge = [1,1, -1, -1]
face_to_face = [-1,-1, -1, -1]
single_cell = [1,1]

output_folder = "result_2_parameter_compare2D3D/"
output_folder_scenarios = "result_mcmc_testing/"
# --------------------------------------------------------------------------------------------------
# Define main function.
# -------------------------------------------------------------------------------------------------- 
if __name__ == "__main__":
  #debug_mode = len(sys.argv) > 1 and sys.argv[1] == "True"
  #procid=int(os.environ['SLURM_PROCID'])
  #array_task_id=int(os.environ["SLURM_ARRAY_TASK_ID"])
  #ntasks=int(os.environ["SLURM_NTASKS"])
  #counter=array_task_id * ntasks + procid
  #print(f"hello from worker number {counter} with procid: {procid} and array_task_id: {array_task_id}", flush=True)
  
  #number_of_cores = int(os.environ['TASKS_PER_NODE'])#SLURM_CPUS_PER_TASK
  number_of_cores = cpu_count()
  #print(number_of_cores, flush=True)
  #print(int(os.environ['JOB_CPUS_PER_NODE']))
  #exit(1)
  test_name    = 'particles_test_goethite_illite'

  ##distances    = [ "bulk_distance", "average_distance", "particle_sizes" ]

  domain_sizes  = [1]#[5,4,3,2,1]#5,4,3,2,1]#[50,15,10]# ,5,4,3,2,1
  #domain_sizes = [100]
  scenarios = [[[1,1,1], [1,1,1]], [[1,1,1], [1,1,2]], [[1,1,1], [1,1,3]], [[1,1,1], [1,1,5]], [[1,1,1], [1,1,10]], [[1,1,1], [1,1,15]]]#,[illite_fine_3D, goethite_fine_3D], [illite_fine_3D, goethite_coarse_3D], [illite_medium_3D, goethite_fine_3D], [illite_medium_3D, goethite_coarse_3D]]#, 
              #[illite_fine_3D, goethite_fine_3D], [illite_fine_3D, goethite_coarse_3D], [illite_medium_3D, goethite_fine_3D], [illite_medium_3D, goethite_coarse_3D]]#, 
              #[[1,1,1], [1,2,2]], [[1,1,1], [1,3,3]], [[1,1,1], [1,4,4]], [[1,1,1], [1,5,5]],
              #[[1,1,2], [2,2,2]], [[1,1,3], [3,3,3]], [[1,1,4], [4,4,4]], [[1,1,5], [5,5,5]],[[2,2,2], [4,4,4]], [[2,2,2], [4,4,2]], [[2,2,2], [4,4,6]], [[2,2,2], [4,4,8]]]
  #scenarios =  [[illite_fine_3D, goethite_fine_3D], [illite_fine_3D, goethite_coarse_3D], [illite_medium_3D, goethite_fine_3D], [illite_medium_3D, goethite_coarse_3D]]
  #scenarios = [[[2,2,17],[6,6,2]],[[2,2,34],[6,6,2]],[[2,2,17],[30,30,2]],[[2,2,34],[30,30,2]]]
  #scenarios = [[illite_fine_3D, goethite_fine_3D]]
  #[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
  scenarios = [[[1,1,1], [1,1,1]]]
  dist = [0.25,0.5, 0.75]
  dist = [0.25]
  #jp =[10,15]
  #jp = [5,10]
  #jp = [10]
  nsteps = [10] #, 20

  #scenario = [[illite_fine, goethite_fine], [illite_fine, goethite_coarse], [illite_medium, goethite_fine], [illite_medium, goethite_coarse]]
  #scenario_names = ["illite_fine_goethite_fine", "illite_fine_goethite_coarse", "illite_medium_goethite_fine", "illite_medium_goethite_coarse"]
  print(sys.argv)
  n_runs = [int(int(sys.argv[1]) / len(scenarios))]
  n_runs = [60,61,62,63]
  
  #n_runs = range(1,10)
  methods = ['mh', 'am', 'dr', 'dram']
  
  #s_n  = int(sys.argv[1]) % len(scenarios)
  #scenarios = [scenarios[s_n]]
  #if(s_n <= 5):
    #jp = [5,10]
  #else:
    #jp = [10,15]
  #jp = [20,15]
  jp = [5]
  print("int(sys.argv[1]) ", int(sys.argv[1])," n_runs ", n_runs[0], " scenarios[0] ",scenarios[0], " jp ", jp, " domain_sizes ", domain_sizes, flush=True) 
  #parameter_minmax_ = [[j, j], [dist, d]]
  fun_args  = []
  base_test = getattr(mcmc_test, test_name)
  print("start loop", flush=True)
  for n in nsteps:
    for type_i, type_g in scenarios:
      for ds in domain_sizes:
        nx = [50,50,ds]
        #print(type_i, type_g, flush=True)
        output_folder_scenarios = output_folder + "/nsteps_" + str(n) + "/domain_" + str(nx) + "/" + "type_i_" + str(type_i) + "_type_g_" + str(type_g) + "/"
        for d in dist: 
          for j in jp:
            for run in n_runs:
              filepath = output_folder_scenarios + 'distribution_'+ str(d)+"_jump_parameter_"+ str(j) + "_run_" + str(run)
              if(min(type_i) <= min(nx) and min(type_g) <= min(nx)):
                # try:
                  # #print(filepath+ "/chainfile.txt")
                  # with open(filepath+ "/chainfile.txt", 'r') as fp:
                    # for count, line in enumerate(fp):
                        # pass
                    # if(count+1 < 4000):
                        # print("no_count ", count+1, " ", filepath)
                        # try:
                          # os.remove(output_folder_scenarios + 'distribution_'+ str(d)+"_jump_parameter_"+ str(j) + "_run_" + str(run) +"/chainfile.txt")
                          # os.remove(output_folder_scenarios + 'distribution_'+ str(d)+"_jump_parameter_"+ str(j) + "_run_" + str(run) +"/covchainfile.txt")
                          # os.remove(output_folder_scenarios + 'distribution_'+ str(d)+"_jump_parameter_"+ str(j) + "_run_" + str(run) +"/s2chainfile.txt")
                          # os.remove(output_folder_scenarios + 'distribution_'+ str(d)+"_jump_parameter_"+ str(j) + "_run_" + str(run) +"/sschainfile.txt")
                        # except:
                           # print("not remove")
                        # fun_args.append( base_test(nx = [50,50,ds], type_i = type_i  ,charges_i = uniform_negative_3D, type_g = type_g, charges_g = uniform_positive_3D,
  			# distribution = d,jump_parameter = j,  n_steps = n,
  			# distance_fct    = "particle_sizes",
  			# ca_settings = [True, True, False, False, False],
  			# #ca_settings = [True, True, True, False, True],
  			# #subaggregate_threshold = 0.5,
  			# file_name    = "nsteps_" + str(n) + "_domain_" + str([50,50, ds]) + "type_i_" + str(type_i) + "_type_g_" + str(type_g) + 'distribution_'+ str(d)+"_jump_parameter_"+ str(j) +"_run_" + str(run),
  			# filepath = output_folder_scenarios + 'distribution_'+ str(d)+"_jump_parameter_"+ str(j) + "_run_" + str(run))) 


                # except:
                  # print("no ", filepath)
                factor = 2 
                parameter_minmax_ = [[[j, j], [d, d]],[[j*factor, j*factor], [d, d]],[[j, j], [d*factor, d*factor]],[[j*factor, j*factor], [d*factor, d*factor]]] 
                #parameter_minmax_ = [[[j, j], [d, d]], [[j, j], [d*factor, d*factor]]] 
                fun_args.append( base_test(nx = nx, type_i = type_i  ,charges_i = uniform_negative_3D, type_g = type_g, charges_g = uniform_positive_3D,
  			    distribution = d,jump_parameter = j,  n_steps = n,
  			    distance_fct    = "particle_sizes",
                	    nsimu = 4000,
  			    ca_settings = [True, True, False, False, False],
                parameter_minmax = parameter_minmax_[run-n_runs[0]],
                method = 'dram',
  			    #ca_settings = [True, True, True, False, True],
  			    #subaggregate_threshold = 0.5,
  			    file_name    = "nsteps_" + str(n) + "_domain_" + str(nx) + "type_i_" + str(type_i) + "_type_g_" + str(type_g) + 'distribution_'+ str(d)+"_jump_parameter_"+ str(j) +"_run_" + str(run),
  			    filepath = output_folder_scenarios + 'distribution_'+ str(d)+"_jump_parameter_"+ str(j) + "_run_" + str(run))) 
		
   
  
  print(len(fun_args), flush=True)
  #run_test_from_class(fun_args[0])  
  #exit(1)
  # for n in nsteps:
  #   for d in dist: 
  #     for j in jp:
  #       fun_args.append( base_test(nx = [100,100], type_i = illite_medium  ,charges_i = uniform_negative,type_g = goethite_coarse, charges_g = uniform_positive,
  #   distribution = d,jump_parameter = j,  n_steps = n,
  #   distance_fct    = "particle_sizes",
  #   #ca_settings = [True, True, True, False, True],
  #   #subaggregate_threshold = 0.5,
  #   file_name    = 'distribution_'+ str(d)+"_jump_parameter_"+ str(j) +"_n_steps_" + str(n),
  #   filepath = output_folder_scenarios +'result_mcmcm_distribution_'+ str(d)+"_jump_parameter_"+ str(j) + "_n_steps_" + str(n)))                
  
  

  #pool = Pool(nodes=number_of_cores)#nodes=number_of_cores  processes=number_of_cores
  pool = Pool(ncpus=number_of_cores)
  pool.map(run_test_from_class, fun_args)#chunksize=1


  #run_test_from_class(fun_args[0])
  processes = []
  #for fun_arg in fun_args:
    #t = multiprocessing.Process(target=run_test_from_class, args=(fun_arg,))
    #processes.append(t)
    #t.start()


  #while multiprocessing.active_children():
    #val = input("Enter your value: ")
    #if val == "kill_all_children":
     # active = multiprocessing.active_children()
     # for child in active:
      #  child.kill()
     # time.sleep(2)

  #for one_process in processes:
    #one_process.join()
  
  print("End of script", flush=True)
  

  
