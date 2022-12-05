# distutils: language=c++

cdef class PythonClassName :
  cdef CythonClassName *thisptr # hold a C++ instance which we're wrapping
  def __cinit__(self, porosity, jump_param_singles, jump_param_composites="default"):
    if jump_param_composites == "default": jump_param_composites = jump_param_singles
    self.thisptr = new CythonClassName (porosity, jump_param_singles, jump_param_composites)
  def __dealloc__(self):
    del self.thisptr
  def move_particles(self):
    return self.thisptr.move_particles ()
  def fields(self):
    return self.thisptr.fields ()
  @staticmethod
  def print_array(vec):
    print_array (vec)
  @staticmethod
  def bulk_distance(vec_a, vec_b):
    return bulk_distance(vec_a, vec_b)
