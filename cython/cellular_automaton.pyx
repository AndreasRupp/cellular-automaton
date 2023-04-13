#distutils : language = c++

cdef class PythonClassName: 
  cdef CythonClassName* thisptr #hold a C++ instance which we're wrapping 
  def __cinit__(self, jump_param_composites = "default" ):
    if jump_param_composites   == "default" : jump_param_composites = 5
    self.thisptr = new CythonClassName(jump_param_composites) 
  def __dealloc__(self): 
    del self.thisptr 
  def print_array(self): 
    self.thisptr.print_array() 
  def fields(self): 
    return self.thisptr.fields() 
  def placeSingleCellBURandomly(self, _porosity, _jump_parameter, random_seed): 
    if _porosity == "default" : _porosity = 0.5 
    if _jump_parameter == "default": _jump_parameter = 1.0 
    if random_seed == "default" : random_seed = 0 
    self.thisptr.placeSingleCellBURandomly(_porosity, _jump_parameter, random_seed) 
  def placeSphere(self, _position, _radius, _jump_parameter): 
    if _jump_parameter == "default" : _jump_parameter = 1.0 
    if _radius == "default" : _radius = 1.0 
    if _position == "default" : _position = -1 
    return self.thisptr.placeSphere(_position, _radius, _jump_parameter) 
  def placePlane(self, _position, _extent, _jump_parameter): 
    if _jump_parameter == "default" : _jump_parameter = 1.0 
    if _extent == "default" : _extent = [ 0 ]* 5 
    if _position == "default": _position = -1 
    return self.thisptr.placePlane(_position, _extent, _jump_parameter) 
  def placeParticles(self): 
    self.thisptr.placeParticles() 
  def doCAM(self) : 
    self.thisptr.doCAM()
  def average_particle_size(self):
    return self.thisptr.average_particle_size()
  def particle_size_distribution(self) : 
    return self.thisptr.particle_size_distribution()
  def bulk_distance(self,vec_a, vec_b):
    return self.thisptr.bulk_distance (vec_a, vec_b)
#  @staticmethod
#cdef class CellularAutomaton:
#cdef CellularAutomatonCP* ptr #hold a C++ instance which we're wrapping
#def __cinit__(self):
#self.ptr = new CellularAutomatonCP()
#def __dealloc__(self):
#del self.ptr
