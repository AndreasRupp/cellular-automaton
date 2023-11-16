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

# goethite_fine_smaller1x = [ 2, 16]
# goethite_fine_smaller2x = [ 2,15]
# goethite_fine_bigger1x = [ 2, 18]
# goethite_fine_bigger2x = [ 2, 19]

# illite_fine_3D = [ 2, 2,6]
# illite_medium_3D = [ 2, 2,30]
# goethite_fine_3D = [2,2, 17]
# goethite_coarse_3D = [2,2,34]

illite_fine_3D = [ 1, 1,3]
illite_medium_3D = [ 1, 1,15]
goethite_fine_3D = [1,1, 8]
goethite_coarse_3D = [1,1,17]

# illite_fine_3D = [ 1, 1,2]
# illite_medium_3D = [ 1, 1,8]
# goethite_fine_3D = [1,1, 4]
# goethite_coarse_3D = [1,1,9]
my_nx = [50,50,50]
# test_i_3D = [ 1, 1,2]
# test_g_3D = [1,1,8]

# illite_fine_3D = [ 1, 1,1]
# illite_medium_3D = [ 1, 1,3]
# goethite_fine_3D = [1,1, 2]
# goethite_coarse_3D = [1,1,4]
uniform_positive = [1] * 4 
uniform_negative = [-1] * 4


uniform_positive_3D = [1] * 6 
uniform_negative_3D = [-1] * 6 
face_to_edge = [1,1, -1, -1]
face_to_face = [-1,-1, -1, -1]
single_cell = [1,1]
#output_folder_scenarios = "final_results_max_acc_3D_quarter_size/"#"final_results_max_accuracy_max_iteration/"#"final_results_max_acc_3D_half_size_v2/"##  ##"final_results_3D_100_100_mcmc_4000_100_/"  #
#output_folder_scenarios = "result_mcmc_subaggregates_goethite_fine_face_to_face_v1"
#output_folder_scenarios = "result_mcmc_porosities"
#output_folder_scenarios = "result_mcmc_3D_nx_50_v1_jp9"
#output_folder_scenarios = "result_mcmc_2D_nx_100_particle_size_illite_medium_goethite_coarse__100_40_v2/"#_illite_medium_goethite_fine
output_folder_scenarios = "result_mcmc_2D_nx_100_particle_size_illite_medium_goethite_coarse__100_40_v2/"#_illite_medium_goethite_fine
output_folder_scenarios = "result_mcmc_testing/"#_illite_medium_goethite_fine
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
  #subaggregateThresholds = [100, 0.01]#1, 5,0.5, 10, 0.1, 
  subaggregateThresholds = [1.999,2.001]
  porosities = [0.8,0.85, 0.9,0.95]
  n_choose_bins = range(2, 30, 4)
  domain_sizes  = [100 ] #10, 25,  50, 
  time_points   = [  0, 10, 25,  50 ]
  dist = [0,0.25,0.5,0.75]
  jp =[0, 10,15, 20,25, 30]
  #jp = [0,5,7.5,10,12.5,15]
  nsteps = [10]#10
  scenario = [[illite_fine, goethite_fine], [illite_fine, goethite_coarse], [illite_medium, goethite_fine], [illite_medium, goethite_coarse]]
  scenario_names = ["illite_fine_goethite_fine", "illite_fine_goethite_coarse", "illite_medium_goethite_fine", "illite_medium_goethite_coarse"]
  fun_args  = []
  base_test = getattr(mcmc_test, test_name)
  #--------------------------------------------------
  # for j in range(len(scenario)):# TODO das hier mal mit 50x50 machen
  #   i = j
  #   s = scenario[i]
  #   s_name = scenario_names[i]
  #   print(s_name)
  #   print(s[0])
  #   print(s[1])
  #   for d in dist:
  #     for j in jp:
  #       for n in nsteps:
  #         fun_args.append( base_test(n_steps = n,distribution = d,jump_parameter = j,type_i = s[0], type_g = s[1],
  #          file_name    = s_name +'_distribution_'+ str(d)+"_jump_parameter_"+ str(j)  +"_n_steps_" + str(n),
  #          filepath =  "result_mcmc_particle_size_via_bu_"+ s_name + "/" + s_name +'_distribution_'+ str(d)+"_jump_parameter_"+ str(j)  +"_n_steps_" + str(n)))





#different porosities
# for i in range(len(scenario)): 
#   s = scenario[i]
#   s_name = scenario_names[i]
#   for p in porosities:
#     fun_args.append( base_test(n_steps = nsteps[0],type_i = s[0], type_g = s[1], porosity = p, 
#            file_name    = s_name +'_porosity_'+ str(p),
#             filepath =  output_folder_scenarios+ "/" + s_name +'_porosity_'+ str(p)))




#different subaggreagte thresholds
  # for sub_thr in  subaggregateThresholds: 
  #   for i in range(10):
  #     #print(i)
  #     fun_args.append( base_test(
  #     #distance_fct = [ "bulk_distance", "average_distance", "particle_sizes" ],
  #     #ecdf_type    = [ "standard",      "standard",         "standard"       ],
  #     #n_bins       = [ 5,               3,                   5               ],
  #     #n_steps      = [ i,               i,                   i               ],
  #     #subset_sizes = [100] * 40,
  #     type_i = goethite_fine, 
  #     type_g = goethite_fine,
  #     charges_i = face_to_face,
  #     charges_g = face_to_face * (-1),
  #     subaggregate_threshold = sub_thr,
  #     ca_settings = [True, True, False, False, True],
  #     n_steps = i,
  #     file_name    = 'n_steps_' + str(i),
  #     filepath = output_folder_scenarios+ "/"+'result_mcmcm_thresh_'+ str(sub_thr) + '_n_steps_' + str(i)
  #     ) )
  #   for j in range(10):
  #     i = (j+1) * 10
  #     fun_args.append( base_test(
  #     # distance_fct = [ "bulk_distance", "average_distance", "particle_sizes" ],
  #     # ecdf_type    = [ "standard",      "standard",         "standard"       ],
  #     # n_bins       = [ 5,               3,                   5               ],
  #     # n_steps      = [ i,               i,                   i               ],
  #     type_i = goethite_fine, 
  #     type_g = goethite_fine,
  #     charges_i = face_to_face,
  #     charges_g = face_to_face * (-1),
  #     subaggregate_threshold = sub_thr,
  #     n_steps = i,
  #     ca_settings = [True, True, False, False, True],
  #     file_name    = 'n_steps_' + str(i),
  #     filepath = output_folder_scenarios +"/"+'result_mcmcm_thresh_'+ str(sub_thr) + '_n_steps_' + str(i)
  #     ) )
  # print("ads----------------------------------------------------------")










  #-------------------------------------------------
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
  # for j in range(30):
  #   i = (j+26) * 10
  #   fun_args.append( base_test(
  #     # distance_fct = [ "bulk_distance", "average_distance", "particle_sizes" ],
  #     # ecdf_type    = [ "standard",      "standard",         "standard"       ],
  #     # n_bins       = [ 5,               3,                   5               ],
  #     # n_steps      = [ i,               i,                   i               ],
  #    n_steps = i,
  #    file_name    = 'n_steps_' + str(i),
  #    filepath = output_folder_scenarios +'result_mcmcm_n_steps_' + str(i)
  #    ) )
  # for d in dist:
  #   for n in nsteps:
  #     for j in jp:
  #       fun_args.append( base_test(n_steps = n,distribution = d,jump_parameter = j,
  #         file_name    = 'distribution_'+ str(d)+"_jump_parameter_"+ str(j)  +"_n_steps_" + str(n),
  #         filepath = output_folder_scenarios +'result_mcmcm_distribution_'+ str(d)+"_jump_parameter_"+ str(j)  +"_n_steps_" + str(n)))
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

  # fun_args.append( base_test(type_i = goethite_fine ,type_g = goethite_fine, file_name    = 'goethite_fine_2x', filepath    = output_folder_scenarios+'goethite_fine_2x'))
  # fun_args.append( base_test(type_i = illite_fine ,type_g = illite_fine, file_name    = 'illite_fine_2x', filepath    = output_folder_scenarios+'illite_fine_2x'))

  # fun_args.append( base_test(type_i = goethite_fine_smaller1x ,type_g = goethite_fine, file_name    = 'goethite_fine_smaller1x', filepath    = output_folder_scenarios+'goethite_fine_smaller1x'))
  # fun_args.append( base_test(type_i = goethite_fine_smaller2x ,type_g = goethite_fine, file_name    = 'goethite_fine_smaller2x', filepath    = output_folder_scenarios+'goethite_fine_smaller2x'))

  # fun_args.append( base_test(type_i = goethite_fine_bigger1x ,type_g = goethite_fine, file_name    = 'goethite_fine_bigger1x', filepath    = output_folder_scenarios+'goethite_fine_bigger1x'))
  # fun_args.append( base_test(type_i = goethite_fine_bigger2x ,type_g = goethite_fine, file_name    = 'goethite_fine_bigger2x', filepath    = output_folder_scenarios+'goethite_fine_bigger2x'))
  # fun_args.append( base_test(type_i = illite_fine ,type_g = goethite_fine, file_name    = 'illite_fine_goethite_fine', filepath    = output_folder_scenarios+'result_illite_fine_goethite_fine'))
  # fun_args.append( base_test(type_i = illite_fine ,type_g = goethite_coarse, file_name    = 'illite_fine_goethite_coarse', filepath    = output_folder_scenarios+'result_illite_fine_goethite_coarse'))
  # fun_args.append( base_test(type_i = illite_medium ,type_g = goethite_fine, file_name    = 'illite_medium_goethite_fine', filepath    = output_folder_scenarios+'result_illite_medium_goethite_fine'))
  # fun_args.append( base_test(type_i = illite_medium ,type_g = goethite_coarse, file_name    = 'illite_medium_goethite_coarse', filepath    = output_folder_scenarios+'result_illite_medium_goethite_coarse'))
#--------3D---------------
  # for i in range(5):
  #   fun_args.append( base_test(jump_parameter   = 10, subaggregate_threshold = 0.99,n_steps =10, ca_settings = [True, True, True, False, False], 
  #     nx = my_nx, type_i = illite_medium_3D ,charges_i = uniform_negative_3D,type_g = goethite_coarse_3D, charges_g = uniform_positive_3D,
  #     file_name    = 'illite_medium_goethite_coarse_3D', filepath    =output_folder_scenarios+"/"+'result_illite_medium_goethite_coarse_3D_' + str(i)))
  #   fun_args.append( base_test(jump_parameter   = 10, subaggregate_threshold = 0.99,n_steps = 10, ca_settings = [True, True, True, False, False], 
  #     nx = my_nx, type_i = illite_medium_3D ,charges_i = uniform_negative_3D,type_g = goethite_fine_3D, charges_g = uniform_positive_3D,
  #     file_name    = 'illite_medium_goethite_fine_3D', filepath    =output_folder_scenarios+"/"+'result_illite_medium_goethite_fine_3D_' + str(i)))


  # fun_args.append( base_test(nx = [100, 100, 100], type_i = illite_fine_3D ,charges_i = uniform_negative_3D, type_g = goethite_fine_3D,charges_g = uniform_positive_3D, file_name    = 'illite_fine_goethite_fine_3D', filepath    = output_folder_scenarios + "/"+'result_illite_fine_goethite_fine_3D'))
  # fun_args.append( base_test(nx = [100, 100, 100], type_i = illite_fine_3D ,charges_i = uniform_negative_3D,type_g = goethite_coarse_3D, charges_g = uniform_positive_3D,file_name    = 'illite_fine_goethite_coarse_3D', filepath    = output_folder_scenarios +"/"+'result_illite_fine_goethite_coarse_3D'))
  # fun_args.append( base_test(nx = [100, 100, 100], type_i = illite_medium_3D ,charges_i = uniform_negative_3D,type_g = goethite_fine_3D, charges_g = uniform_positive_3D,file_name    = 'illite_medium_goethite_fine_3D', filepath    = output_folder_scenarios +"/"+'result_illite_medium_goethite_fine_3D'))
  # fun_args.append( base_test(nx = [100, 100, 100], type_i = illite_medium_3D ,charges_i = uniform_negative_3D,type_g = goethite_coarse_3D, charges_g = uniform_positive_3D,file_name    = 'illite_medium_goethite_coarse_3D', filepath    = output_folder_scenarios +"/"+'result_illite_medium_goethite_coarse_3D'))

  # fun_args.append( base_test(nx = my_nx, type_i = illite_fine_3D ,charges_i = uniform_negative_3D,type_g = goethite_fine_3D, charges_g = uniform_positive_3D, file_name    = 'illite_fine_goethite_fine_3D', filepath    = output_folder_scenarios+ "/"+'result_illite_fine_goethite_fine_3D'))
  # fun_args.append( base_test(nx = my_nx, type_i = illite_fine_3D ,charges_i = uniform_negative_3D,type_g = goethite_coarse_3D, charges_g = uniform_positive_3D,file_name    = 'illite_fine_goethite_coarse_3D', filepath    = output_folder_scenarios+"/"+'result_illite_fine_goethite_coarse_3D'))
  # fun_args.append( base_test(nx = my_nx, type_i = illite_medium_3D ,charges_i = uniform_negative_3D,type_g = goethite_fine_3D,charges_g = uniform_positive_3D, file_name    = 'illite_medium_goethite_fine_3D', filepath    = output_folder_scenarios+"/"+'result_illite_medium_goethite_fine_3D'))
  # fun_args.append( base_test(nx = my_nx, type_i = illite_medium_3D ,charges_i = uniform_negative_3D,type_g = goethite_coarse_3D, charges_g = uniform_positive_3D,file_name    = 'illite_medium_goethite_coarse_3D', filepath    =output_folder_scenarios+"/"+'result_illite_medium_goethite_coarse_3D'))
# for n in nsteps:
  # fun_args.append( base_test(n_steps = n,nx = [25, 25, 25], type_i = illite_fine_3D ,charges_i = uniform_negative_3D,type_g = goethite_fine_3D, 
  #   charges_g = uniform_positive_3D, file_name    = 'illite_fine_goethite_fine_3D_n_steps_' + str(n), filepath    = output_folder_scenarios+ 'result_illite_fine_goethite_fine_3D_n_steps_' + str(n)))
  # fun_args.append( base_test(n_steps = n,nx = [25, 25, 25], type_i = illite_fine_3D ,charges_i = uniform_negative_3D,type_g = goethite_coarse_3D, 
  #   charges_g = uniform_positive_3D,file_name    = 'illite_fine_goethite_coarse_3D_n_steps_' + str(n), filepath    = output_folder_scenarios+'result_illite_fine_goethite_coarse_3D_n_steps_' + str(n)))
  # fun_args.append( base_test(n_steps = n,nx = [25, 25, 25], type_i = illite_medium_3D ,charges_i = uniform_negative_3D,type_g = goethite_fine_3D,
  #     charges_g = uniform_positive_3D, file_name    = 'illite_medium_goethite_fine_3D_n_steps_' + str(n), filepath    = output_folder_scenarios+'result_illite_medium_goethite_fine_3D_n_steps_' + str(n)))
  # fun_args.append( base_test(n_steps = n,nx = [25, 25, 25], type_i = illite_medium_3D ,charges_i = uniform_negative_3D,type_g = goethite_coarse_3D, 
  #     charges_g = uniform_positive_3D,file_name    = 'illite_medium_goethite_coarse_3D_n_steps_' + str(n), filepath    =output_folder_scenarios+'result_illite_medium_goethite_coarse_3D_n_steps_' + str(n)))
  #for n in nsteps:
  #for n in nsteps:
	#  for d in dist:
	#	  for j in jp:
	#		  fun_args.append( base_test(nx = my_nx, type_i = illite_medium_3D ,charges_i = uniform_negative_3D,type_g = goethite_coarse_3D, charges_g = uniform_positive_3D, 
	#		distribution = d,jump_parameter = j,  n_steps = n,
	#		distance_fct    = "particle_sizes",
	#		ca_settings = [True, True, True, False, True],
	#		subaggregate_threshold = 0.5,
	#		file_name    = 'distribution_'+ str(d)+"_jump_parameter_"+ str(j) +"_n_steps_" + str(n),
	#		filepath = output_folder_scenarios +'result_mcmcm_distribution_'+ str(d)+"_jump_parameter_"+ str(j) + "_n_steps_" + str(n)))    
  for n in nsteps:
	  for d in dist: 
		  for j in jp:
			  fun_args.append( base_test(nx = [100,100], type_i = illite_medium  ,charges_i = uniform_negative,type_g = goethite_coarse, charges_g = uniform_positive,
			distribution = d,jump_parameter = j,  n_steps = n,
			distance_fct    = "particle_sizes",
			ca_settings = [True, True, True, False, True],
			subaggregate_threshold = 0.5,
			file_name    = 'distribution_'+ str(d)+"_jump_parameter_"+ str(j) +"_n_steps_" + str(n),
			filepath = output_folder_scenarios +'result_mcmcm_distribution_'+ str(d)+"_jump_parameter_"+ str(j) + "_n_steps_" + str(n)))    

  # for i in range(10):
  #   #print(i)
  #   fun_args.append( base_test(
  #     nx = my_nx, type_i = illite_medium_3D ,charges_i = uniform_negative_3D,type_g = goethite_coarse_3D,charges_g = uniform_positive_3D,
  #     #distance_fct = [ "bulk_distance", "average_distance", "particle_sizes" ],
  #     #ecdf_type    = [ "standard",      "standard",         "standard"       ],
  #     #n_bins       = [ 5,               3,                   5               ],
  #     #n_steps      = [ i,               i,                   i               ],
  #     #subset_sizes = [75] * 40,
  #     n_steps = i,
  #     file_name    = 'illite_medium_goethite_coarse_n_steps_' + str(i),
  #     filepath = output_folder_scenarios+ "/"+ 'result_mcmcm_illite_medium_goethite_coarse_n_steps_' + str(i)
  #    ) )

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
  #run_test_from_class(fun_args[0])
  for one_process in processes:
    one_process.join()
  print(len(fun_args))

  
