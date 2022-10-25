from libcpp.vector cimport vector
from libcpp.string cimport string

IncludeFiles

cdef extern from "<CAM/cellular_automaton.hxx>" :
  cdef cppclass CythonClassName C++ClassName :
    CythonClassName ( porosity , jump_param ) except +
    const vector[ CyReplace03 ]& move_particles ()
    void plot_solution (const vector[ unsigned int ]&)
