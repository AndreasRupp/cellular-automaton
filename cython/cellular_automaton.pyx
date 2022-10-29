# distutils: language=c++

cdef class PythonClassName :
  cdef CythonClassName *thisptr # hold a C++ instance which we're wrapping
  def __cinit__(self, porosity, jump_param):
    self.thisptr = new CythonClassName (porosity, jump_param)
  def __dealloc__(self):
    del self.thisptr
  def move_particles(self):
    return self.thisptr.move_particles ()
  def fields(self):
    return self.thisptr.fields ()
  def print_array(self, vec):
    print_array (vec)
  def bulk_distance(self, vec_a, vec_b):
    return bulk_distance(vec_a, vec_b)
