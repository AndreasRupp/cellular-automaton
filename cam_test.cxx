#include <iomanip>
#include <iostream>

#include <CAM/cellular_automaton.hxx>

template <std::size_t size>
void print_array(const std::array<unsigned int, size> particles)
{
  const unsigned int nx = std::sqrt(size);
  for (unsigned int y = 0; y < nx; ++y)
  {
    for (unsigned int x = 0; x < nx; ++x)
      if (particles[y * nx + x] == 0)
        std::cout << std::setw(2) << std::setfill('0') << particles[y * nx + x] << "  ";
      else
        std::cout << "\033[0;31m" << std::setw(2) << std::setfill('0') << particles[y * nx + x]
                  << "\033[0m"
                  << "  ";
    std::cout << std::endl;
  }
}

int main()
{
  constexpr unsigned int nx = 10;
  constexpr unsigned int ny = 10;
  const unsigned int n_moves = 3;
  const double porosity = 0.96;
  const double jump_param = 4.;
  cellular_automaton<nx, ny> domain(porosity, jump_param);

  print_array(domain.fields());

  for (unsigned int i = 0; i < n_moves; ++i)
  {
    std::cout << std::endl;
    print_array(domain.move_particles());
  }
}