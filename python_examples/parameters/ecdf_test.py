import numpy as np
import os, sys


class basic_test:
  def __init__( self,
    nx             = [50, 50],
    porosity       = 0.7,
    n_steps        = 5,
    jump_parameter = 5,

    # subset_sizes    = [100] * 40,
    subset_sizes    = [50] * 10,
    min_value_shift = 0.1,
    max_value_shift = -0.1,
    n_choose_bins   = "default",
    jump_params     = "default",

    distance_fct = "default",
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
    self.subset_sizes    = subset_sizes
    self.n_choose_bins   = n_choose_bins
    self.max_value_shift = max_value_shift
    self.min_value_shift = min_value_shift

    if jump_params == "default":  self.jump_params = range(jump_parameter-5, jump_parameter+6)
    else:                         self.jump_params = jump_params

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
        avg_part_b = self.PyCAM.average_particle_size(b)
        return np.abs(avg_part_a - avg_part_b)
      self.distance_fct = my_distance
    elif distance_fct == "particle_sizes":
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
        part_sizes_a = self.PyCAM.particle_size_distribution(a)
        part_sizes_b = self.PyCAM.particle_size_distribution(b)
        for item in part_sizes_b:
          part_sizes_a.append(item)
        return part_sizes_a
      self.distance_fct = my_distance


    if n_choose_bins == "default":
      if distance_fct == "default" or distance_fct == "bulk_distance":
        self.n_choose_bins = 20
      elif distance_fct == "average_distance":
        self.n_choose_bins = 8
    else:
      self.n_choose_bins = n_choose_bins
