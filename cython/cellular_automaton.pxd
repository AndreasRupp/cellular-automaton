from libcpp.vector cimport vector 
from libcpp.string cimport string 
from libcpp cimport bool

IncludeFiles


cdef extern from "<CAM/cam_interface.hxx>" : 
  cdef cppclass CythonClassName C++ClassName : 
    CythonClassName() except + 
    void print_array() 
    void placeSingleCellBURandomly( double _porosity,double _jump_parameter, unsigned int random_seed)
    void doCAM() 
    const vector[unsigned int] &fields() 
    bool placeSphere(int _position, double _radius, double _jump_parameter) 
    bool placePlane(int _position,vector[unsigned int] _extent, double _jump_parameter) 
    void placeParticles()

#cdef extern from "<CAM/cellular_automaton.hxx>":
#cdef cppclass CellularAutomatonCP \ "CAM::CellularAutomaton< std::array<unsigned int, 2>({10,10}), std::vector<unsigned int> >":
#CellularAutomatonCP() except +
#void apply "CAM::CellularAutomaton<NX_string>::apply"(CythonClassName);
