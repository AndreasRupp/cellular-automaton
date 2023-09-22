from __future__ import print_function

from datetime import datetime
import multiprocessing
import os, sys, time
import numpy as np
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

illite_fine_3D = [ 2, 2,6]
illite_medium_3D = [ 2, 2,30]
goethite_fine_3D = [2,2, 17]
goethite_coarse_3D = [2,2,34]

illite_fine_3D = [ 1, 1,3]
illite_medium_3D = [ 1, 1,15]
goethite_fine_3D = [1,1, 8]
goethite_coarse_3D = [1,1,17]

single_cell = [1,1]
output_folder_scenarios = "final_results_3D_half_size/"  #_max_accuracy_max_iteration
#version 2 bei mit 0 40 jump
# --------------------------------------------------------------------------------------------------
# Define main function.
# -------------------------------------------------------------------------------------------------- 
if __name__ == "__main__":
  debug_mode = len(sys.argv) > 1 and sys.argv[1] == "True"
  
  test_name    = 'particles_test_goethite_illite'
  #test_name    = 'particles_test'
  distances    = [ "bulk_distance", "average_distance", "particle_sizes" ]
  n_bins       = [ 20,              8,                  10               ]

  n_choose_bins = range(2, 30, 4)
  domain_sizes  = [100 ] #10, 25,  50, 
  time_points   = [  0, 10, 25,  50 ]


  fun_args  = []
  base_test = getattr(mcmc_test, test_name)
  
  # run_test_from_class(base_test())
  # #run_test_from_class(base_test(n_steps = 1, file_name    = 'n_steps_1'))
  # for i in range(10):
  #   #print(i)
  #   fun_args.append( base_test(
  #     #distance_fct = [ "bulk_distance", "average_distance", "particle_sizes" ],
  #     #ecdf_type    = [ "standard",      "standard",         "standard"       ],
  #     #n_bins       = [ 5,               3,                   5               ],
  #     #n_steps      = [ i,               i,                   i               ],
  #     #subset_sizes = [100] * 40,
  #     n_steps = i,
  #     file_name    = 'n_steps_' + str(i),
  #     filepath = output_folder_scenarios+ 'result_mcmcm_n_steps_' + str(i)
  #     ) )
  # for j in range(5):
  #   i = (j+1) * 10
  #   fun_args.append( base_test(
  #     #distance_fct = [ "bulk_distance", "average_distance", "particle_sizes" ],
  #     #ecdf_type    = [ "standard",      "standard",         "standard"       ],
  #     #n_bins       = [ 5,               3,                   5               ],
  #     #n_steps      = [ i,               i,                   i               ],
  #     n_steps = i,
  #     file_name    = 'n_steps_' + str(i),
  #     filepath = output_folder_scenarios +'result_mcmcm_n_steps_' + str(i)
  #     ) )

  # fun_args.append( base_test(n_steps = 0, file_name    = 'n_steps_0', filepath = 'result_mcmc_n_steps_0'))
  # fun_args.append( base_test(n_steps = 1, file_name    = 'n_steps_1', filepath = 'result_mcmc_n_steps_0'))
  # fun_args.append( base_test(n_steps = 2, file_name    = 'n_steps_2', filepath = 'result_mcmc_n_steps_0'))
  # fun_args.append( base_test(n_steps = 3, file_name    = 'n_steps_3', filepath = 'result_mcmc_n_steps_0'))
  # fun_args.append( base_test(n_steps = 4, file_name    = 'n_steps_4', filepath = 'result_mcmc_n_steps_0'))
  # fun_args.append( base_test(n_steps = 5, file_name    = 'n_steps_5', filepath = 'result_mcmc_n_steps_0'))
  # fun_args.append( base_test(n_steps = 6, file_name    = 'n_steps_6', filepath = 'result_mcmc_n_steps_0'))
  # fun_args.append( base_test(n_steps = 7, file_name    = 'n_steps_7', filepath = 'result_mcmc_n_steps_0'))
  # fun_args.append( base_test(n_steps = 8, file_name    = 'n_steps_8', filepath = 'result_mcmc_n_steps_0'))
  # fun_args.append( base_test(n_steps = 9, file_name    = 'n_steps_9', filepath = 'result_mcmc_n_steps_0'))
  # fun_args.append( base_test(n_steps = 10, file_name    = 'n_steps_10', filepath = 'result_mcmc_n_steps_0'))
  # fun_args.append( base_test(n_steps = 20, file_name    = 'n_steps_20', filepath = 'result_mcmc_n_steps_0'))
  # fun_args.append( base_test(n_steps = 30, file_name    = 'n_steps_30', filepath = 'result_mcmc_n_steps_0'))
  # fun_args.append( base_test(n_steps = 40, file_name    = 'n_steps_40', filepath = 'result_mcmc_n_steps_0'))
  # fun_args.append( base_test(n_steps = 50, file_name    = 'n_steps_50', filepath = 'result_mcmc_n_steps_0'))

  #fun_args.append( base_test(nsimu = htop, adaptint= 500, file_name    = 'nsim_10000_adaptint_500'))
  #fun_args.append( base_test(nsimu = 5000, adaptint= 250, file_name    = 'nsim_5000_adaptint_250'))
  #fun_args.append( base_test(nsimu = 2500, adaptint= 125, file_name    = 'nsim_2500_adaptint_125'))
  #fun_args.append( base_test(nsimu = 1000, adaptint= 50, file_name    = 'nsim_1000_adaptint_50'))
  #fun_args.append( base_test(nsimu = 500, adaptint= 25, file_name    = 'nsim_500_adaptint_25'))
  


  # fun_args.append( base_test(subset_sizes = [100] * 40, file_name    = 'subset_sizes_100_40'))
  # fun_args.append( base_test(subset_sizes = [100] * 20, file_name    = 'subset_sizes_100_20'))
  # fun_args.append( base_test(subset_sizes = [100] * 10, file_name    = 'subset_sizes_100_10'))
  # fun_args.append( base_test(subset_sizes = [75] * 40, file_name    = 'subset_sizes_75_40'))
  # fun_args.append( base_test(subset_sizes = [50] * 40, file_name    = 'subset_sizes_50_40'))
  # fun_args.append( base_test(subset_sizes = [25] * 40, file_name    = 'subset_sizes_25_40'))
  # fun_args.append( base_test(subset_sizes = [10] * 40, file_name    = 'subset_sizes_10_40'))
  # fun_args.append( base_test(subset_sizes = [10] * 4, file_name    = 'subset_sizes_10_4'))
  # fun_args.append( base_test(subset_sizes = [20] * 8, file_name    = 'subset_sizes_20_8'))
  # fun_args.append( base_test(subset_sizes = [40] * 16, file_name    = 'subset_sizes_40_16'))
  # fun_args.append( base_test(subset_sizes = [60] * 24, file_name    = 'subset_sizes_60_24'))
  # fun_args.append( base_test(subset_sizes = [80] * 32, file_name    = 'subset_sizes_80_32'))


  # fun_args.append( base_test(type_i = illite_fine ,type_g = goethite_fine, file_name    = 'illite_fine_goethite_fine', filepath    = output_folder_scenarios+'result_illite_fine_goethite_fine'))
  # fun_args.append( base_test(type_i = illite_fine ,type_g = goethite_coarse, file_name    = 'illite_fine_goethite_coarse', filepath    = output_folder_scenarios+'result_illite_fine_goethite_coarse'))
  # fun_args.append( base_test(type_i = illite_medium ,type_g = goethite_fine, file_name    = 'illite_medium_goethite_fine', filepath    = output_folder_scenarios+'result_illite_medium_goethite_fine'))
  # fun_args.append( base_test(type_i = illite_medium ,type_g = goethite_coarse, file_name    = 'illite_medium_goethite_coarse', filepath    = output_folder_scenarios+'result_illite_medium_goethite_coarse'))
#--------3D---------------
  # fun_args.append( base_test(nx = [100, 100, 100], type_i = illite_fine_3D ,type_g = goethite_fine_3D, file_name    = 'illite_fine_goethite_fine_3D', filepath    = 'final_results/result_illite_fine_goethite_fine_3D'))
  # fun_args.append( base_test(nx = [100, 100, 100], type_i = illite_fine_3D ,type_g = goethite_coarse_3D, file_name    = 'illite_fine_goethite_coarse_3D', filepath    = 'final_results/result_illite_fine_goethite_coarse_3D'))
  # fun_args.append( base_test(nx = [100, 100, 100], type_i = illite_medium_3D ,type_g = goethite_fine_3D, file_name    = 'illite_medium_goethite_fine_3D', filepath    = 'final_results/result_illite_medium_goethite_fine_3D'))
  # fun_args.append( base_test(nx = [100, 100, 100], type_i = illite_medium_3D ,type_g = goethite_coarse_3D, file_name    = 'illite_medium_goethite_coarse_3D', filepath    = 'final_results/result_illite_medium_goethite_coarse_3D'))

  fun_args.append( base_test(nx = [50, 50, 50], type_i = illite_fine_3D ,type_g = goethite_fine_3D, file_name    = 'illite_fine_goethite_fine_3D', filepath    = 'final_results/result_illite_fine_goethite_fine_3D'))
  fun_args.append( base_test(nx = [50, 50, 50], type_i = illite_fine_3D ,type_g = goethite_coarse_3D, file_name    = 'illite_fine_goethite_coarse_3D', filepath    = 'final_results/result_illite_fine_goethite_coarse_3D'))
  fun_args.append( base_test(nx = [50, 50, 50], type_i = illite_medium_3D ,type_g = goethite_fine_3D, file_name    = 'illite_medium_goethite_fine_3D', filepath    = 'final_results/result_illite_medium_goethite_fine_3D'))
  fun_args.append( base_test(nx = [50, 50, 50], type_i = illite_medium_3D ,type_g = goethite_coarse_3D, file_name    = 'illite_medium_goethite_coarse_3D', filepath    = 'final_results/result_illite_medium_goethite_coarse_3D'))
  
  for i in range(10):
    #print(i)
    fun_args.append( base_test(
      nx = [50, 50, 50], type_i = illite_medium_3D ,type_g = goethite_coarse_3D,
      #distance_fct = [ "bulk_distance", "average_distance", "particle_sizes" ],
      #ecdf_type    = [ "standard",      "standard",         "standard"       ],
      #n_bins       = [ 5,               3,                   5               ],
      #n_steps      = [ i,               i,                   i               ],
      #subset_sizes = [100] * 40,
      n_steps = i,
      file_name    = 'n_steps_' + str(i),
      filepath = output_folder_scenarios+ 'result_mcmcm_n_steps_' + str(i)
      ) )

  # fun_args.append( base_test(
  #   distance_fct = [ "bulk_distance", "average_distance", "particle_sizes" ],
  #   ecdf_type    = [ "standard",      "standard",         "standard"       ],
  #   n_bins       = [ 5,               3,                   5               ],
  #   n_steps      = [ 5,               5,                   5               ],
  #   subset_sizes = [100] * 40,
  #   file_name    = 'multiple_3_nsteps5_subsetsizes_100_40'
  #   ) )

  # fun_args.append( base_test(
  # distance_fct = [ "bulk_distance", "average_distance", "particle_sizes" ],
  # ecdf_type    = [ "standard",      "standard",         "standard"       ],
  # n_bins       = [ 5,               3,                   5               ],
  # n_steps      = [ 10,               10,                   10               ],
  # subset_sizes = [100] * 40                      ,
  # file_name    = 'multiple_3_nsteps10_subsetsizes_100_40'
  # ) )

  # fun_args.append( base_test(
  # distance_fct = [ "bulk_distance", "average_distance", "particle_sizes" ],
  # ecdf_type    = [ "standard",      "standard",         "standard"       ],
  # n_bins       = [ 5,               3,                   5               ],
  # n_steps      = [ 20,               20,                   20               ],
  # subset_sizes = [100] * 40                      ,
  # file_name    = 'multiple_3_nsteps20_subsetsizes_100_40'
  # ) )

  # fun_args.append( base_test(
  # distance_fct = [ "bulk_distance", "average_distance", "particle_sizes" ],
  # ecdf_type    = [ "standard",      "standard",         "standard"       ],
  # n_bins       = [ 5,               3,                   5               ],
  # n_steps      = [ 5,               5,                   5              ],
  # subset_sizes = [50] * 40                      ,
  # file_name    = 'multiple_3_nsteps5_subsetsizes_50_40'
  # ) )

  # fun_args.append( base_test(
  # distance_fct = [ "bulk_distance", "average_distance", "particle_sizes" ],
  # ecdf_type    = [ "standard",      "standard",         "standard"       ],
  # n_bins       = [ 5,               3,                   5               ],
  # n_steps      = [ 5,               5,                   5              ],
  # subset_sizes = [25] * 40                      ,
  # file_name    = 'multiple_3_nsteps5_subsetsizes_25_40'
  # ) )
  # fun_args.append( base_test(
  # distance_fct = [ "bulk_distance", "average_distance", "particle_sizes" ],
  # ecdf_type    = [ "standard",      "standard",         "standard"       ],
  # n_bins       = [ 5,               3,                   5               ],
  # n_steps      = [ 5,               5,                   5              ],
  # subset_sizes = [10] * 40                      ,
  # file_name    = 'multiple_3_nsteps5_subsetsizes_10_40'
  # ) )



  # fun_args.append( base_test(
  # distance_fct = [ "bulk_distance", "average_distance", "particle_sizes" ],
  # ecdf_type    = [ "standard",      "standard",         "standard"       ],
  # n_bins       = [ 5,               3,                   5               ],
  # n_steps      = [ 5,               5,                   5              ],
  # subset_sizes = [50] * 20                      ,
  # file_name    = 'multiple_3_nsteps5_subsetsizes_50_20'
  # ) )

  # fun_args.append( base_test(
  # distance_fct = [ "bulk_distance", "average_distance", "particle_sizes" ],
  # ecdf_type    = [ "standard",      "standard",         "standard"       ],
  # n_bins       = [ 5,               3,                   5               ],
  # n_steps      = [ 5,               5,                   5              ],
  # subset_sizes = [50] * 10                      ,
  # file_name    = 'multiple_3_nsteps5_subsetsizes_50_10'
  # ) )
  # fun_args.append( base_test(
  # distance_fct = [ "bulk_distance", "average_distance", "particle_sizes" ],
  # ecdf_type    = [ "standard",      "standard",         "standard"       ],
  # n_bins       = [ 5,               3,                   5               ],
  # n_steps      = [ 5,               5,                   5              ],
  # subset_sizes = [50] * 5                      ,
  # file_name    = 'multiple_3_nsteps5_subsetsizes_50_5'
  # ) )
  # fun_args.append( base_test(
  # distance_fct = [ "bulk_distance", "average_distance", "particle_sizes" ],
  # ecdf_type    = [ "standard",      "standard",         "standard"       ],
  # n_bins       = [ 5,               3,                   5               ],
  # n_steps      = [ 5,               5,                   5              ],
  # subset_sizes = [25] * 10                      ,
  # file_name    = 'multiple_3_nsteps5_subsetsizes_25_10'
  # ) )


  #find best---------------------
  # fun_args.append( base_test(
  #   distance_fct = [ "bulk_distance", "average_distance", "particle_sizes" ],
  #   ecdf_type    = [ "standard",      "standard",         "standard"       ],
  #   n_bins       = [ 5,               3,                   5               ],
  #   n_steps      = [ 10,               10,                   10               ],
  #   subset_sizes = [25] * 10,
  #   file_name    = 'multiple_3_nsteps10_subsetsizes_25_10'
  #   ) )
  # fun_args.append( base_test(
  #   distance_fct = [ "bulk_distance", "average_distance", "particle_sizes" ],
  #   ecdf_type    = [ "standard",      "standard",         "standard"       ],
  #   n_bins       = [ 5,               3,                   5               ],
  #   n_steps      = [ 10,               10,                   10               ],
  #   subset_sizes = [50] * 20,
  #   file_name    = 'multiple_3_nsteps10_subsetsizes_50_20'
  #   ) )
  # fun_args.append( base_test(
  #   distance_fct = [ "bulk_distance", "average_distance", "particle_sizes" ],
  #   ecdf_type    = [ "standard",      "standard",         "standard"       ],
  #   n_bins       = [ 5,               3,                   5               ],
  #   n_steps      = [ 10,               10,                   10               ],
  #   subset_sizes = [10] * 40,
  #   file_name    = 'multiple_3_nsteps10_subsetsizes_10_40'
  #   ) )  

# distance_fct = [ "bulk_distance", "average_distance", "particle_sizes" ],
#               ecdf_type    = [ "standard",      "standard",         "standard"       ],
#               n_bins       = [ 5,               3,                   5               ],
#               n_steps      = [ 10,               10,                   10               ],
  #vary mcmc ---------------
  # method= ['mh',  'dram']
  # cov = [ 0.1 ,  1,  10, 100 ]
  # sigma_2 = [10 ** 2, 1, 0.1 **2] 
  # subsetsizes = [[50] * 40]#, [50]*40,
  # adaptinterval = [10,50]
  # for sub in subsetsizes:
  #   for m in method:
  #     for c in cov:
  #       for s in sigma_2:    
  #         for ad in adaptinterval:
  #           fun_args.append( base_test(
  #             subset_sizes = sub,
  #             qcov = np.eye(2)*c,
  #             method = m,
  #             sigma2 = s,
  #             adaptint = ad,
  #             file_name    = 'particle_size_subset_sizes_' + str(sub[0])+ '_method_' + m + "_cov_" + str(c) + "_sigma_" +  str(s)+ "_ad_" + str(ad),
  #             filepath = 'vary_mcmc_v2/result_mcmcm_particle_size_subset_sizes_' + str(sub[0]) + "_" + str(len(sub)) + '_method_' + m + "_cov_" + str(c) + "_sigma_" + str(s) +  "_ad_" + str(ad),
  #             ) )


  processes = []
  for fun_arg in fun_args:
    t = multiprocessing.Process(target=run_test_from_class, args=(fun_arg,))
    processes.append(t)
    t.start()





  while multiprocessing.active_children():
    val = input("Enter your value: ")
    if val == "kill_all_children":
      active = multiprocessing.active_children()
      for child in active:
        child.kill()
      time.sleep(2)
  #run_test_from_class(fun_args[3])
  for one_process in processes:
    one_process.join()
  print(len(fun_args))

  
