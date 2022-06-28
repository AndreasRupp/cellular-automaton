#include "mex.hpp"
#include "mexAdapter.hpp"

#include <CAM/cellular_automaton.hxx>

static constexpr unsigned int nx = 10;
static constexpr unsigned int ny = nx;

template <unsigned int nx, unsigned int ny>
constexpr std::array<int, 4> cellular_automaton<nx, ny>::particle::direct_neigh_;

class MexFunction : public matlab::mex::Function
{
 public:
  void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs)
  {
    checkArguments(outputs, inputs);

    const unsigned int n_moves = inputs[0][0];
    const double porosity = inputs[1][0];
    const double jump_param = inputs[2][0];
    const unsigned int output_rate = inputs[3][0];
    unsigned int random_seed = 0;
    matlab::data::TypedArray<double> results = std::move(inputs[4]);

    cellular_automaton<nx, ny> domain(porosity, jump_param, inputs.size() == 6 ? inputs[5][0] : 0);

    for (unsigned int k = 0; k < nx * ny; ++k)
      results[k][0] = (domain.fields())[k];

    for (unsigned int i = 0; i < n_moves; ++i)
    {
      domain.move_particles();

      if (output_rate != 0 && (i + 1) % output_rate == 0)
        for (unsigned int k = 0; k < nx * ny; ++k)
          results[k][i + 1] = (domain.fields())[k];
    }

    outputs[0] = std::move(results);
  }

  void checkArguments(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs)
  {
    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
    matlab::data::ArrayFactory factory;

    if (inputs.size() != 5 && inputs.size() != 6)
    {
      matlabPtr->feval(
        u"error", 0,
        std::vector<matlab::data::Array>({factory.createScalar("Four inputs required")}));
    }

    // if (inputs[0].getNumberOfElements() != 1) {
    //     matlabPtr->feval(u"error",
    //         0, std::vector<matlab::data::Array>({ factory.createScalar("Input multiplier must be
    //         a scalar") }));
    // }

    // if (inputs[0].getType() != matlab::data::ArrayType::DOUBLE ||
    //     inputs[0].getType() == matlab::data::ArrayType::COMPLEX_DOUBLE) {
    //     matlabPtr->feval(u"error",
    //         0, std::vector<matlab::data::Array>({ factory.createScalar("Input multiplier must be
    //         a noncomplex scalar double") }));
    // }

    // if (inputs[1].getType() != matlab::data::ArrayType::DOUBLE ||
    //     inputs[1].getType() == matlab::data::ArrayType::COMPLEX_DOUBLE) {
    //     matlabPtr->feval(u"error",
    //         0, std::vector<matlab::data::Array>({ factory.createScalar("Input matrix must be type
    //         double") }));
    // }

    // if (inputs[1].getDimensions().size() != 2) {
    //     matlabPtr->feval(u"error",
    //         0, std::vector<matlab::data::Array>({ factory.createScalar("Input must be m-by-n
    //         dimension") }));
    // }
  }
};
