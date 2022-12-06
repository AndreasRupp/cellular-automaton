from libcpp.vector cimport vector
from libcpp.string cimport string

IncludeFiles

cdef extern from "<CAM/cellular_automaton.hxx>" :
  cdef cppclass CythonClassName C++ClassName :
    CythonClassName ( const double , const double ) except +
    CythonClassName ( const double , const double, const double ) except +
    const vector[ unsigned int ]& move_particles ()
    const vector[ unsigned int ]& fields ()

cdef extern from "<CAM/domain.hxx>" :
  void print_array "CAM::print_array<NX_string>" (const vector[ unsigned int ]&)
  unsigned int bulk_distance "CAM::bulk_distance" (const vector[ unsigned int]&, const vector[ unsigned int ]&)
  double average_particle_size "CAM::average_particle_size<NX_string>" (const vector[ unsigned int ]&)
