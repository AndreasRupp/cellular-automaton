from libcpp.vector cimport vector 
from libcpp.string cimport string 
from libcpp cimport bool

IncludeFiles


cdef extern from "<CAM/cam_interface.hxx>" : 
  cdef cppclass CythonClassName C++ClassName : 
    CythonClassName(const double jump_param_composite) except + 
    void print_array() 
    const vector[unsigned int] &fields() 
    void do_CAM() 
    void place_singleCellBU_randomly( double _porosity,double _jump_parameter, unsigned int random_seed)
    bool place_sphere(int _position, double _radius, double _jump_parameter) 
    bool place_plane(int _position,vector[unsigned int] _extent, double _jump_parameter) 
    void place_particles()
    double average_particle_size()
    vector[unsigned int] particle_size_distribution()
    unsigned int bulk_distance (const vector[ unsigned int]&, const vector[ unsigned int ]&)
