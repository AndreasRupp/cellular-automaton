from dataclasses import dataclass
import numpy as np


@dataclass
class basic_test:
  # Configure the cellular automaton method (CAM).
  nx             = [50, 50]
  porosity       = 0.3
  n_steps        = 5
  jump_parameter = 5

  nf = np.prod(nx)  # This is the number of fields in the CAM and can be derived from nx.

  # Configure the eCDF method.
  n_iter          = 2000                                        # Number of dataset sizes.
  jump_params     = range(jump_parameter-5, jump_parameter+6)   # The checked jump parameters.
  bins            = range(int(0.375*nf),int(0.475*nf),5)        # All radii that are checked.
  n_choose_bins   = np.min([20, len(bins)])                     # Number of selected radii for algo.
  subset_sizes    = [40] * 50                                   # Multiplies to n_iter!
  min_value_shift = 0.1                                         # Cutoff value for small values.
  max_value_shift = -0.1                                        # 1 + "cutoff value for large val."
