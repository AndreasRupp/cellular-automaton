from libcpp.vector cimport vector
from libcpp.string cimport string

IncludeFiles

cdef extern from "<CAM/cellular_automaton.hxx>" :
  cdef cppclass CythonClassName C++ClassName :
    CythonClassName ( const double , const double ) except +
    const vector[ unsigned int ]& move_particles ()

cdef extern from "<CAM/domain.hxx>" :
  void print_array "CAM::print_array<std::array<unsigned int, 2>({10, 10})>" (const vector[ unsigned int ]&)
