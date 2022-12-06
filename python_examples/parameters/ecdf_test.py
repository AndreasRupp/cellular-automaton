import numpy as np
import os, sys


class basic_test:
  def __init__( self,
    nx             = [50, 50],
    porosity       = 0.7,
    n_steps        = 5,
    jump_parameter = 5,

    subset_sizes    = [40] * 50,
    min_value_shift = 0.1,
    max_value_shift = -0.1,
    n_choose_bins   = "default",
    jump_params     = "default",
    bins            = "default",

    distance_fct = "default",
    debug_mode   = False,
    file_name    = "basic_test",
    is_plot      = 0
    ):
    # Configure the cellular automaton method (CAM).
    self.nx             = nx
    self.porosity       = porosity
    self.n_steps        = n_steps
    self.jump_parameter = jump_parameter

    n_fields = np.prod(nx)

    # Configure the eCDF method.
    self.subset_sizes    = subset_sizes
    self.n_choose_bins   = n_choose_bins
    self.max_value_shift = max_value_shift
    self.min_value_shift = min_value_shift

    if jump_params == "default":  self.jump_params = range(jump_parameter-5, jump_parameter+6)
    else:                         self.jump_params = jump_params
    if bins == "default":  self.bins = range(int(0.375*n_fields),int(0.475*n_fields),5)
    else:                  self.bins = bins
    if n_choose_bins == "default":  self.n_choose_bins = np.min([20, len(bins)])
    else:                           self.n_choose_bins = n_choose_bins 

    self.debug_mode = debug_mode
    self.file_name  = file_name
    self.is_plot    = is_plot

    if distance_fct == "default" or distance_fct == "bulk_distance":
      try:
        import CAM
      except (ImportError, ModuleNotFoundError) as error:
        sys.path.append(os.path.dirname(os.path.abspath(__file__)) + os.sep + ".."  + os.sep + 
          ".." + os.sep + "import")
        import CAM
      const             = CAM.config()
      const.nx          = self.nx
      const.debug_mode  = self.debug_mode
      self.PyCAM        = CAM.include(const)
      self.distance_fct = self.PyCAM.bulk_distance
    elif distance_fct == "average_distance":
      try:
        import CAM
      except (ImportError, ModuleNotFoundError) as error:
        sys.path.append(os.path.dirname(os.path.abspath(__file__)) + os.sep + ".."  + os.sep + 
          ".." + os.sep + "import")
        import CAM
      const             = CAM.config()
      const.nx          = self.nx
      const.debug_mode  = self.debug_mode
      self.PyCAM        = CAM.include(const)
      def my_distance(a, b):
        avg_part_a = self.PyCAM.average_particle_size(a)
        avg_part_b = self.PyCAM.average_particle_size(a)
        return np.abs(avg_part_a - avg_part_b)
      self.distance_fct = my_distance
