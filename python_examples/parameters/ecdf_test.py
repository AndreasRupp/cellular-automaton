from dataclasses import dataclass
import numpy as np

@dataclass
class ecdf_parameter:
  # Configure the cellular automaton method (CAM).
  _nx = [50, 50]
   
  _nf : int = np.prod(_nx)  # This is the number of fields in the CAM and can be derived from nx.
  porosity      = 0.3
  n_steps        = 5
  _jump_parameter = 5

  # Configure the eCDF method.
  n_iter          = 2000                                        # Number of dataset sizes.
  values          = [0] * 10                                    # Length defines tested jump_param.
  _bins            = range(int(0.375*_nf),int(0.475*_nf),5)     # All radii that are checked.
  n_choose_bins   = 20                                          # Number of selected radii for algo.
  subset_sizes    = [40] * 50                                   # Multiplies to n_iter!
  min_value_shift = 0.1                                         # Cutoff value for small values.
  max_value_shift = -0.1            
  #TODO consistency: make all variables defined like that          
  @property
  def nx(self):
    return self._nx
  @property
  def jump_parameter(self):
    return self._jump_parameter
  @property
  def nf(self):
    self._nf = np.prod(self._nx)
    return  self._nf
  @property
  def bins(self):
    self._bins = range(int(0.375*self.nf),int(0.475*self.nf),5) 
    return self._bins
  @nx.setter
  def nx(self,value):
    self._nx = value
  @jump_parameter.setter
  def jump_parameter(self, value):
    self._jump_parameter = value