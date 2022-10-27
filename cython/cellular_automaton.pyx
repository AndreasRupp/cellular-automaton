# distutils: language=c++

cdef class PythonClassName :
  cdef CythonClassName *thisptr # hold a C++ instance which we're wrapping
  def __cinit__(self, porosity, jump_param):
    self.thisptr = new CythonClassName (porosity, jump_param)
  def __dealloc__(self):
    del self.thisptr
  def move_particles(self):
    return self.thisptr.move_particles ()
  def print_array(self, vec):
    print_array (vec)
