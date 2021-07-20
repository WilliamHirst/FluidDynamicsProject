#include "FluidDynamics.hpp"
#include <cmath>
#include <stdlib.h>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <fstream>
using namespace std;

void FluidDynamics::Initialize(double width, double height, double omega)
{
  const double boltzmann_constant = 1.3806503e-23;
  const double nu = 1./3. * (1./omega - 0.5);
  double **nZero = nOne = nTwo = nThree = nFoure = nFive = nSix = nSeven = nEight = nNine = rho = u_x = u_y = createMatrix(width, height);
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
