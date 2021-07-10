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
}

double ** FluidDynamics::createMatrix(double width, double height) {
  // Set up matrix
  double **matrix;
  matrix = new double*[height];
  // Allocate memory
  for(int i=0;i<height;i++)
      matrix[i] = new double[width];

}
