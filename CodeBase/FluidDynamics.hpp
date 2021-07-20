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
      double **nZero;
      double **nOne;
      double **nTwo;
      double **nThree;
      double **nFour;
      double **nFive;
      double **nSix;
      double **nSeven;
      double **nEight;
      double **nNine;
      double **rho;
      double **u_x;
      double **u_y;
      double *equiDens;
      double *Omega;
      double *unit_x;
      double *unit_y;
      double three_Half;
      double nine_Half;
    public:
      void Initialize(double width, double height, double omega);
      double ** createMatrix(int height,int width);
      void find_density(int width, int height);
};

#endif
