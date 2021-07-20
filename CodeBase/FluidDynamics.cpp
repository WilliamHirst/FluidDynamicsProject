#include "FluidDynamics.hpp"
#include <cmath>
#include <stdlib.h>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <fstream>
using namespace std;

void FluidDynamics::Initialize(int width, int height, double omega)
{
  const double boltzmann_constant = 1.3806503e-23;
  const double nu = 1./3. * (1./omega - 0.5);
  nOne = nTwo = nThree = nFoure = nFive = nSix = nSeven = nEight = nNine = rho = u_x = u_y = createMatrix(width, height);
  omega[0] = (double) 4. / 9.;
  omega[1] = omega[2] = omega[3] = omega[4] = (double) 1. / 9.;
  omega[5] = omega[6] = omega[7] = omega[8] = (double) 1. / 36.;
  unit_x[0] = unit_x[2] = unit_x[4] = unit_y[0] = unit_y[1] = unit_y[3] = 0;
  unit_x[1] = unit_x[5] = unit_x[6] = unit_y[2] = unit_y[5] = unit_y[6] = 1;
  unit_x[3] = unit_x[6] = unit_x[7] = unit_y[4] = unit_y[7] = unit_y[8] = -1;
}
void FluidDynamics::find_density(int width, int height)
{
  // Compute what the nine densities would be if the molecules at this site were in thermal equilibrium
  for(int i = 0; i<9; i++){
    double e_dot_u = unit_x[i] * u_x[width][height] + unit_y[i] * u_y[width][height];
    double u_squared = u_x[width][height] * u_x[width][height] + u_y[width][height] * u_y[width][height];
    equiDens[i] = rho[width][height] * Omega[i] * (1 + 3 * e_dot_u + (double) 9 / 2 * e_dot_u * e_dot_u - 3 / 2 * u_squared);
  }
  //Compute density in next time-step

}

double ** FluidDynamics::createMatrix(double width, double height) {
  // Set up matrix
  double **matrix;
  matrix = new double*[height];
  // Allocate memory
  for(int i=0;i<height;i++)
      matrix[i] = new double[width];
  // Set values to zero
  for(int i = 0; i < height; i++){
      for(int j = 0; j < width; j++){
          matrix[i][j] = 0.0;
      }
  }
  return matrix;
}
