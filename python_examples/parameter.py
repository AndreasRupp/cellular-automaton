from __future__ import print_function

from datetime import datetime
import multiprocessing
import os, sys, time

try:
  from parameter_identification import run_test_from_class
except (ImportError, ModuleNotFoundError) as error:
  sys.path.append(os.path.dirname(os.path.abspath(__file__)) + os.sep  + ".." + os.sep +
    "python_functions")
  from parameter_identification import run_test_from_class

try:
  import ecdf_test
except (ImportError, ModuleNotFoundError) as error:
  sys.path.append(os.path.dirname(os.path.abspath(__file__)) + os.sep  + "parameters")
  import ecdf_test


# --------------------------------------------------------------------------------------------------
# Define main function.
# -------------------------------------------------------------------------------------------------- 
if __name__ == "__main__":
  debug_mode = len(sys.argv) > 1 and sys.argv[1] == "True"

  test_name    = 'basic_test'
  
  distances    = [ "bulk_distance", "average_distance", "particle_sizes" ]
  n_bins       = [ 20,              8,                  10               ]
  ecdf_types   = [ "standard" ]
  subset_sizes = [ [100] * 40 ]
  domain_sizes = [ 50 ] 
  sigmas       = [ 5  ]
  time_points  = [ 0,  5, 10, 25,  50]
  dimensions   = [ 1,  2,  3,  4,   5]

  mult_ecdf_types_2 = [ "standard",      "standard",        ]
  mult_distances_2  = [ "bulk_distance", "average_distance" ]
  mult_n_bins_2     = [ 8,               4                  ]

  mult_ecdf_types_3 = [ "standard",      "standard",         "standard"       ]
  mult_distances_3  = [ "bulk_distance", "average_distance", "particle_sizes" ]
  mult_n_bins_3     = [ 5,               2,                   5               ]
  
  
  fun_args  = []
  base_test = getattr(ecdf_test, test_name)

  fun_args.append( base_test(
    distance_fct  = mult_distances_2,
    ecdf_type     = mult_ecdf_types_2,
    n_choose_bins = mult_n_bins_2,
    file_name     = 'multiple_2'
    ) )

  fun_args.append( base_test(
    distance_fct  = mult_distances_3,
    ecdf_type     = mult_ecdf_types_3,
    n_choose_bins = mult_n_bins_3,
    file_name     = 'multiple_3'
    ) )

  for type_index in range(len(ecdf_types)):
    ecdf_type = ecdf_types[type_index]
    subsets   = subset_sizes[type_index]

    for distance_index in range(len(distances)):
      distance      = distances[distance_index]
      n_choose_bins = n_bins[distance_index]
    
      for size in domain_sizes:
        for sigma in sigmas:
          fun_args.append( base_test(
            nx             = [size, size],
            n_choose_bins  = n_choose_bins,
            jump_parameter = sigma,
            distance_fct   = distance,
            ecdf_type      = ecdf_type,
            subset_sizes   = subsets,
            file_name      = ecdf_type + '_' + distance + '_jump-param_' + str(sigma) + \
                             '_size_' + str(size)
            ) )

      # for steps in time_points:
      #   fun_args.append( base_test(
      #     n_steps      = steps,
      #     distance_fct = distance,
      #     ecdf_type    = ecdf_type,
      #     subset_sizes = subsets,
      #     file_name    = ecdf_type + '_' + distance + '_time-steps_' + str(steps)
      #     ) )

      # for dim in dimensions:
      #   fun_args.append( base_test(
      #     nx           = [ 50 for _ in range(dim) ],
      #     distance_fct = distance,
      #     ecdf_type    = ecdf_type,
      #     subset_sizes = subsets,
      #     file_name    = ecdf_type + '_' + distance + '_dimension_' + str(dim)
      #     ) )

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

  for one_process in processes:
    one_process.join()
