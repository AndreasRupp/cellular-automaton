import numpy as np
import os, sys

try:
  import CAM
except (ImportError, ModuleNotFoundError) as error:
  sys.path.append(os.path.dirname(os.path.abspath(__file__)) + os.sep + ".."  + os.sep +
    ".." + os.sep + "import")
  import CAM

illite_fine = [ 2, 6]
illite_medium = [ 2, 30]
goethite_fine = [2, 17]
goethite_coarse = [2,34]
single_cell = [1,1]

uniform_positive = [1] * 4 
uniform_negative = [-1] * 4 

face_to_edge = [1,1, -1, -1]
face_to_face = [-1,-1, -1, -1]

class basic_test:
  def __init__( self,
    nx             = [50, 50],
    porosity       = 0.7,
    n_steps        = 5,
    jump_parameter = 5,
    ecdf_type       = "standard",
    subset_sizes    = [20] * 40,
    min_value_shift = 0.1,
    max_value_shift = -0.1,
    n_bins          = None,
    distance_fct    = "bulk_distance",
    nsimu			= 10000,
    qcov			= np.eye(2),
    adaptint		= 500,
    


    # ["DFACE_ATTRACTIVITY", "DROTATION", "DROTATION_COMPOSITES", "DSTENCIL_4_ALL_BUS", "DSUB_COMPOSITES"]
    ca_settings      = [False, False, False, False, False],
    const_stencil_size = 5, # only used when DSTENCIL_4_ALL_BUS is set
    subaggregate_threshold = 0.01,
    debug_mode   = False,
    file_name    = "basic_test",
    filepath 	 = "basic_test",
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
    self.max_value_shift = max_value_shift
    self.min_value_shift = min_value_shift
    
	# Configure the MCMC method.
    self.nsimu = nsimu
    self.qcov = qcov
    self.adaptint = adaptint 
	
    if n_bins is None:
      if distance_fct == "bulk_distance" or distance_fct == "particle_sizes":  self.n_bins = 20
      elif distance_fct == "average_distance":                                 self.n_bins = 8

    self.ca_settings   = ca_settings
    self.const_stencil_size = const_stencil_size # only used when DSTENCIL_4_ALL_BUS is set
    self.subaggregate_threshold = subaggregate_threshold
    self.debug_mode = debug_mode
    self.file_name  = file_name
    self.filepath	= filepath
    self.is_plot    = is_plot

    const             = CAM.config()
    const.nx          = self.nx
    const.debug_mode  = self.debug_mode
    self.PyCAM        = CAM.include(const)
    self.CAM_wrapper = self.PyCAM(jump_parameter, subaggregate_threshold)

    def return_distance_fct( distance_fct ):
      if distance_fct == "bulk_distance":
        return self.CAM_wrapper.bulk_distance
      elif distance_fct == "average_distance":
        def my_distance(a, b):
          avg_part_a = self.CAM_wrapper.average_particle_size(a)
          avg_part_b = self.CAM_wrapper.average_particle_size(b)
          return np.abs(avg_part_a - avg_part_b)
        return my_distance
      elif distance_fct == "particle_sizes":
        def my_distance(a, b):
          part_sizes_a = self.CAM_wrapper.particle_size_distribution(a)
          part_sizes_b = self.CAM_wrapper.particle_size_distribution(b)
          for item in part_sizes_b:
            part_sizes_a.append(item)
          return part_sizes_a
        return my_distance

    if isinstance(distance_fct, list):
      self.distance_fct = [return_distance_fct(entry) for entry in distance_fct]
    else:
      self.distance_fct = return_distance_fct(distance_fct)


class particles_test:
  def __init__( self,
    nx               = [50, 50],
    porosity         = 0.7,
    n_steps          = 0,
    jump_parameter   = 0,
    
    ecdf_type       = "standard",
    subset_sizes    = [50] * 40,
    min_value_shift = 0.1,
    max_value_shift = -0.1,
    n_bins          = None,
    distance_fct    = "bulk_distance",
    #distance_fct    = "particle_sizes",
    nsimu			= 10000,
    qcov			= np.eye(2),
    adaptint		= 500,
    parameter_minmax = [[0, 10], [0, 1]],

	type_i			= illite_fine, 
	type_g 			= single_cell,
	charges_i 		= uniform_negative,
	charges_g		= uniform_positive,
    # ["DFACE_ATTRACTIVITY", "DROTATION", "DROTATION_COMPOSITES", "DSTENCIL_4_ALL_BUS", "DSUB_COMPOSITES"]
    ca_settings            = [False, False, False, False, False],
    const_stencil_size     = 5, # only used when DSTENCIL_4_ALL_BUS is set
    subaggregate_threshold = 0.01,
    distribution           = 0.5, 
    debug_mode             = False,
    file_name              = "particles_test",
    #filepath 			   = "output_mcmc_particles_basic_10steps",
    filepath 			   = "AHHHH",
    is_plot                = False
    ):
    # Configure the cellular automaton method (CAM).
    self.nx               = nx
    self.porosity         = porosity
    self.n_steps          = n_steps
    self.jump_parameter   = jump_parameter

    n_fields = np.prod(nx)

    # Configure the eCDF method.
    self.ecdf_type       = ecdf_type
    self.subset_sizes    = subset_sizes
    self.n_bins          = n_bins
    self.max_value_shift = max_value_shift
    self.min_value_shift = min_value_shift
    
    # Configure the MCMC method.
    self.nsimu = nsimu
    self.qcov = qcov
    self.adaptint = adaptint 
    self.parameter_minmax = parameter_minmax
	
    if n_bins is None:
      if distance_fct == "bulk_distance" or distance_fct == "particle_sizes":  self.n_bins = 20
      elif distance_fct == "average_distance":                                 self.n_bins = 8
    
    self.ca_settings             = ca_settings
    self.const_stencil_size      = const_stencil_size # only used when DSTENCIL_4_ALL_BUS is set
    self.subaggregate_threshold  = subaggregate_threshold
    self.distribution            = distribution
    self.debug_mode              = debug_mode
    self.file_name               = file_name
    self.filepath				 = filepath
    self.is_plot                 = is_plot
    
    self.type_i		= type_i
    self.type_g		= type_g
    self.charges_i	= charges_i
    self.charges_g	= charges_g
		
    const             = CAM.config()
    const.nx          = self.nx
    const.debug_mode  = self.debug_mode
    const.ca_settings = ca_settings
    self.PyCAM        = CAM.include(const)
    self.CAM_wrapper = self.PyCAM(jump_parameter, subaggregate_threshold)

    def return_distance_fct( distance_fct ):
      if distance_fct == "bulk_distance":
        return self.CAM_wrapper.bulk_distance
      elif distance_fct == "average_distance":
        def my_distance(a, b):
          avg_part_a = self.CAM_wrapper.average_particle_size(a)
          avg_part_b = self.CAM_wrapper.average_particle_size(b)
          return np.abs(avg_part_a - avg_part_b)
        return my_distance
      elif distance_fct == "particle_sizes":
        def my_distance(a, b):
          part_sizes_a = self.CAM_wrapper.particle_size_distribution(a)
          part_sizes_b = self.CAM_wrapper.particle_size_distribution(b)
          for item in part_sizes_b:
            part_sizes_a.append(item)
          return part_sizes_a
        return my_distance

    if isinstance(distance_fct, list):
      self.distance_fct = [return_distance_fct(entry) for entry in distance_fct]
    else:
      self.distance_fct = return_distance_fct(distance_fct)



class particles_test_goethite_illite:
  def __init__( self,
    nx               = [100, 100],
    porosity         = 0.9,
    n_steps          = 20,
    jump_parameter   = 10,

    ecdf_type       = "standard",
    subset_sizes    = [75] * 40,##[40] * 15,#[100] * 40,
    min_value_shift = 0.1,
    max_value_shift = -0.1,
    n_bins          = None,
    #distance_fct    = "bulk_distance",
    distance_fct    = "particle_sizes",
    #distance_fct    = "average_distance",
    nsimu			= 4000,#10000
    qcov			= np.eye(2) * 100,
    adaptint		= 100,#100
    parameter_minmax = [[0, 25], [0, 1]],
    method = 'mh',
    sigma2 = 10** 2,
    N0 = 1,
    S20 = 0,

	  type_i			= illite_fine, 
	  type_g 			= goethite_coarse,
	  charges_i 		= uniform_negative,
	  charges_g		= uniform_positive,
    # ["DFACE_ATTRACTIVITY", "DROTATION", "DROTATION_COMPOSITES", "DSTENCIL_4_ALL_BUS", "DSUB_COMPOSITES"]
    ca_settings            = [True, True, False, False, False],
    const_stencil_size     = 20, # only used when DSTENCIL_4_ALL_BUS is set
    subaggregate_threshold = 0.01,
    distribution           = 0.5, 
    debug_mode             = False,
    file_name              = "particles_test_goethite_illite",
    filepath 			         = "findbest_output_mcmc_particles_coarse_goethite_fine_illite_mixture",
    #filepath = "aaaaah",
    is_plot                = False
    ):
    # Configure the cellular automaton method (CAM).
    self.nx               = nx
    self.porosity         = porosity
    self.n_steps          = n_steps
    self.jump_parameter   = jump_parameter

    n_fields = np.prod(nx)

    # Configure the eCDF method.
    self.ecdf_type       = ecdf_type
    self.subset_sizes    = subset_sizes
    self.n_bins          = n_bins
    self.max_value_shift = max_value_shift
    self.min_value_shift = min_value_shift
	
	# Configure the MCMC method.
    self.nsimu = nsimu
    self.qcov = qcov
    self.adaptint = adaptint
    self.parameter_minmax = parameter_minmax 
    self.method = method
    self.sigma2 = sigma2
    self.N0 = N0
    self.S20 = S20
	
    if n_bins is None:
      if distance_fct == "bulk_distance" or distance_fct == "particle_sizes":  self.n_bins = 20
      elif distance_fct == "average_distance":                                 self.n_bins = 8
	
    self.type_i		= type_i
    self.type_g		= type_g
    self.charges_i	= charges_i
    self.charges_g	= charges_g
	
    self.ca_settings             = ca_settings
    self.const_stencil_size      = const_stencil_size # only used when DSTENCIL_4_ALL_BUS is set
    self.subaggregate_threshold  = subaggregate_threshold
    self.distribution            = distribution
    self.debug_mode              = debug_mode
    self.file_name               = file_name
    self.filepath				 = filepath
    self.is_plot                 = is_plot

    const             = CAM.config()
    const.nx          = self.nx
    const.debug_mode  = self.debug_mode
    const.ca_settings = ca_settings
    self.PyCAM        = CAM.include(const)
    self.CAM_wrapper = self.PyCAM(jump_parameter, subaggregate_threshold)

    def return_distance_fct( distance_fct ):
      if distance_fct == "bulk_distance":
        return self.CAM_wrapper.bulk_distance
      elif distance_fct == "average_distance":
        def my_distance(a, b):
          avg_part_a = self.CAM_wrapper.average_particle_size(a)
          avg_part_b = self.CAM_wrapper.average_particle_size(b)
          return np.abs(avg_part_a - avg_part_b)
        return my_distance
      elif distance_fct == "particle_sizes":
        def my_distance(a, b):          
          part_sizes_a = a[:]#self.CAM_wrapper.particle_size_distribution(a)
          #part_sizes_b = b[:]#self.CAM_wrapper.particle_size_distribution(b)
          for item in b:#part_sizes_b:
            part_sizes_a.append(item)
          #print("##")
          #print(part_sizes_a)
          #print("---")
          return part_sizes_a
        return my_distance

    if isinstance(distance_fct, list):
      self.distance_fct = [return_distance_fct(entry) for entry in distance_fct]
    else:
      self.distance_fct = return_distance_fct(distance_fct)

