import numpy as np
import os, sys

try:
  import CAM
except (ImportError, ModuleNotFoundError) as error:
  sys.path.append(os.path.dirname(os.path.abspath(__file__)) + os.sep + ".."  + os.sep +
    ".." + os.sep + "import")
  import CAM


class basic_test:
  def __init__( self,
    nx             = [50, 50],
    porosity       = 0.7,
    n_steps        = 5,
    jump_parameter = 5,

    ecdf_type       = "standard",
    subset_sizes    = [100] * 40,
    min_value_shift = 0.1,
    max_value_shift = -0.1,
    n_bins          = None,
    n_runs          = 10,
    jump_params     = None,
    distance_fct    = "bulk_distance",

    debug_mode   = False,
    file_name    = "basic_test",
    is_plot      = False
    ):
    # Configure the cellular automaton method (CAM).
    self.nx             = nx
    self.porosity       = porosity
    self.n_steps        = n_steps
    self.jump_parameter = jump_parameter

    n_fields = np.prod(nx)

    # Configure the eCDF method.
    self.ecdf_type       = ecdf_type
    self.subset_sizes    = subset_sizes
    self.n_bins          = n_bins
    self.n_runs          = n_runs
    self.max_value_shift = max_value_shift
    self.min_value_shift = min_value_shift

    if n_bins is None:
      if distance_fct == "bulk_distance" or distance_fct == "particle_sizes":  self.n_bins = 20
      elif distance_fct == "average_distance":                                 self.n_bins = 8

    if jump_params is None:  self.jump_params = range(jump_parameter-5, jump_parameter+6)
    else:                    self.jump_params = jump_params

    self.debug_mode = debug_mode
    self.file_name  = file_name
    self.is_plot    = is_plot

    const             = CAM.config()
    const.nx          = self.nx
    const.debug_mode  = self.debug_mode
    self.PyCAM        = CAM.include(const)

    def return_distance_fct( distance_fct ):
      if distance_fct == "bulk_distance":
        return self.PyCAM.bulk_distance
      elif distance_fct == "average_distance":
        def my_distance(a, b):
          avg_part_a = self.PyCAM.average_particle_size(a)
          avg_part_b = self.PyCAM.average_particle_size(b)
          return np.abs(avg_part_a - avg_part_b)
        return my_distance
      elif distance_fct == "particle_sizes":
        def my_distance(a, b):
          part_sizes_a = self.PyCAM.particle_size_distribution(a)
          part_sizes_b = self.PyCAM.particle_size_distribution(b)
          for item in part_sizes_b:
            part_sizes_a.append(item)
          return part_sizes_a
        return my_distance

    if isinstance(distance_fct, list):
      self.distance_fct = [return_distance_fct(entry) for entry in distance_fct]
    else:
      self.distance_fct = return_distance_fct(distance_fct)
