#include<domain.hxx>

using namespace CAM;

int main(){
  domain<10,10> grid(0.1, 4.);
  grid.move_particles();
}