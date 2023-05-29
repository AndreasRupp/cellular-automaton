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
  def place_single_cell_bu_randomly(self, _porosity, _jump_parameter, random_seed): 
    if _porosity == "default" : _porosity = 0.5 
    if _jump_parameter == "default": _jump_parameter = 1.0 
    if random_seed == "default" : random_seed = 0 
    self.thisptr.place_single_cell_bu_randomly(_porosity, _jump_parameter, random_seed) 
  def place_sphere(self, _position, _radius, _jump_parameter): 
    if _jump_parameter == "default" : _jump_parameter = 1.0 
    if _radius == "default" : _radius = 1.0 
    if _position == "default" : _position = -1 
    return self.thisptr.place_sphere(_position, _radius, _jump_parameter) 
  def place_plane(self, _position, _jump_parameter, _extent): 
    if _jump_parameter == "default" : _jump_parameter = 1.0 
    if _extent == "default" : _extent = [ 0 ]* 3 
    if _position == "default": _position = -1 
    return self.thisptr.place_plane(_position, _jump_parameter, _extent) 
  def place_particles(self): 
    self.thisptr.place_particles() 
  def do_cam(self) : 
    self.thisptr.do_cam()
  def average_particle_size(self):
    return self.thisptr.average_particle_size()
  def particle_size_distribution(self) : 
    return self.thisptr.particle_size_distribution()
  def bulk_distance(self,vec_a, vec_b):
    return self.thisptr.bulk_distance (vec_a, vec_b)

