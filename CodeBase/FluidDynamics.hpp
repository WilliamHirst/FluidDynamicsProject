#ifndef SOLVER_HPP
#define SOLVER_HPP
#define _USE_MATH_DEFINES
#include <cmath>
#include "armadillo"
#include <string>

class FluidDynamics
{
    private:
      double nu;
      double omega;
    public:
      void Initialize(double width, double height, double omega);
};

#endif
