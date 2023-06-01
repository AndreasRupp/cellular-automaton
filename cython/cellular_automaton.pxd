from libcpp.vector cimport vector 
from libcpp.string cimport string 
from libcpp cimport bool

IncludeFiles


cdef extern from "<CAM/cam_interface.hxx>" : 
  cdef cppclass CythonClassName C++ClassName : 
    CythonClassName(const double jump_param_composite) except + 
    void print_array() 
    const vector[unsigned int] &fields() 
    void do_cam() 
    void place_single_cell_bu_randomly( double _jump_parameter, double _porosity, unsigned int random_seed)
    bool place_sphere(double _jump_parameter, double _radius, int _position) 
    bool place_plane(double _jump_parameter, vector[unsigned int] _extent, int _position) 
    void place_particles()
    double average_particle_size()
    vector[unsigned int] particle_size_distribution()
    unsigned int bulk_distance (const vector[ unsigned int]&, const vector[ unsigned int ]&)
