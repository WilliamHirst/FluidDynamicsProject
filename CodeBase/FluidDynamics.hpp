#ifndef SOLVER_HPP
#define SOLVER_HPP
#define _USE_MATH_DEFINES
#include <cmath>
#include <armadillo>
#include <string>

class FluidDynamics
{
    private:
      double nu;
      double c;
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
      double finalTime;
      double dt;
    public:
      void Initialize(int height, int width, double omega, double initial_vel, double fTime);
      double ** createMatrix(int height,int width);
      void find_density(int height, int width);
      void deleteMatrix(double **matrix, int height);
      void Lattice_Boltzmann(int height, int width);
      void printDensToFile(int height, int width, std::ofstream &ofile);
};

#endif
