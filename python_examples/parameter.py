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

  n_choose_bins = range(2, 30, 4)
  domain_sizes  = [ 10, 25, 50, 100 ] 
  time_points   = [  0, 10, 25,  50 ]


  fun_args  = []
  base_test = getattr(ecdf_test, test_name)

  # fun_args.append( base_test(
  #   distance_fct = [ "bulk_distance", "average_distance" ],
  #   ecdf_type    = [ "standard",      "standard",        ],
  #   n_bins       = [ 8,               4                  ],
  #   n_runs       = 100,
  #   file_name    = 'multiple_2'
  #   ) )

  fun_args.append( base_test(
    distance_fct = [ "bulk_distance", "bulk_distance"    ],
    ecdf_type    = [ "standard",      "standard",        ],
    n_bins       = [ 20,               4                 ],
    n_steps      = [ 5,                10                ],
    file_name    = 'multiple_times'
    ) )

  # fun_args.append( base_test(
  #   distance_fct = [ "bulk_distance", "average_distance", "particle_sizes" ],
  #   ecdf_type    = [ "standard",      "standard",         "standard"       ],
  #   n_bins       = [ 5,               3,                   5               ],
  #   n_runs       = 100,
  #   file_name    = 'multiple_3'
  #   ) )

  # for size in domain_sizes:
  #   fun_args.append( base_test(
  #     nx             = [size, size],
  #     file_name      = 'size_' + str(size)
  #     ) )

  # for steps in time_points:
  #   fun_args.append( base_test(
  #     n_steps      = steps,
  #     file_name    = 'time-steps_' + str(steps)
  #     ) )

  # for n_choose_bin in n_bins:
  #   fun_args.append( base_test(
  #     n_bins    = n_choose_bin,
  #     file_name = 'n-bins_' + str(n_choose_bin)
  #     ) )

  # for distance_index in range(len(distances)):
  #   distance      = distances[distance_index]
  #   n_choose_bins = n_bins[distance_index]
  #   fun_args.append( base_test(
  #     n_bins       = n_choose_bins,
  #     distance_fct = distance,
  #     file_name    = "distance_" + distance
  #   ) )

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
