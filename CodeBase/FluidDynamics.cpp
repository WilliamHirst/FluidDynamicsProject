#include "FluidDynamics.hpp"
#include <cmath>
#include <stdlib.h>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <fstream>
using namespace std;

void FluidDynamics::Initialize(int height, int width, double omega, double initial_vel, double fTime)
{
  nu = 1./3. * (1./omega - 0.5);
  c = 1;
  finalTime = fTime;
  dt = (double) 1 / 30;
  three_Half = (double) 3 / 2;
  nine_Half = (double) 9 / 2;
  nOne = nTwo = nThree = nFour = nFive = nSix = nSeven = nEight = nNine = rho = u_x = u_y = createMatrix(height, width);
  Omega[0] = (double) 4. / 9.;
  Omega[1] = Omega[2] = Omega[3] = Omega[4] = (double) 1. / 9.;
  Omega[5] = Omega[6] = Omega[7] = Omega[8] = (double) 1. / 36.;
  unit_x[0] = unit_x[2] = unit_x[4] = unit_y[0] = unit_y[1] = unit_y[3] = 0;
  unit_x[1] = unit_x[5] = unit_x[6] = unit_y[2] = unit_y[5] = unit_y[6] = 1;
  unit_x[3] = unit_x[6] = unit_x[7] = unit_y[4] = unit_y[7] = unit_y[8] = -1;
  //Initialize rho, u_x and u_y.
  for(int i = 0; i < height; i++){
    for(int j = 0; j < width; j++) {
      rho[i][j] = rho[i][j] = 1;
      u_x[i][j] = initial_vel;
      u_y[i][j] = 0;
    }
  }
}
void FluidDynamics::find_density(int height, int width)
{
  // Compute what the nine densities would be if the molecules at this site were in thermal equilibrium
  for(int i = 0; i<9; i++){
    double e_dot_u = unit_x[i] * u_x[height][width] + unit_y[i] * u_y[height][width];
    double u_squared = u_x[height][width] * u_x[height][width] + u_y[height][width] * u_y[height][width];
    equiDens[i] = rho[height][width] * Omega[i] * (1 + 3 * e_dot_u +  nine_Half * e_dot_u * e_dot_u - three_Half * u_squared);
  }
  //Compute density in next time-step for all 9 lattices.
  nZero[height][width] = nZero[height][width] +  omega * (equiDens[0] - nZero[height][width]);
  nOne[height][width] = nOne[height][width] +  omega * (equiDens[1] - nOne[height][width]);
  nTwo[height][width] = nTwo[height][width] +  omega * (equiDens[2] - nTwo[height][width]);
  nThree[height][width] = nThree[height][width] +  omega * (equiDens[3] - nThree[height][width]);
  nFour[height][width] = nFour[height][width] +  omega * (equiDens[4] - nFour[height][width]);
  nFive[height][width] = nFive[height][width] +  omega * (equiDens[5] - nFive[height][width]);
  nSix[height][width] = nSix[height][width] +  omega * (equiDens[6] - nSix[height][width]);
  nSeven[height][width] = nSeven[height][width] +  omega * (equiDens[7] - nSeven[height][width]);
  nEight[height][width] = nEight[height][width] +  omega * (equiDens[8] - nEight[height][width]);
  nNine[height][width] = nNine[height][width] +  omega * (equiDens[9] - nNine[height][width]);
}
void FluidDynamics::Lattice_Boltzmann(int height, int width)
{
  std::ofstream ofile;
  std::string outfilename = "particleDensity.txt";
  ofile.open(outfilename);
  ofile << std::setprecision(16)<< width << " "<< height <<endl;
  int counter = 0;
  while(counter<finalTime) {
    for(int i = 0; i<height; i++){
      for(int j=0; j<width; j++){
        find_density(i,j);
      }
    }
    printDensToFile(height, width, ofile);
  }
  ofile.close();
}

double ** FluidDynamics::createMatrix(int height, int width) {
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

void FluidDynamics::deleteMatrix(double **matrix, int height){
  for (int i = 0; i < height; i++)
    delete[] matrix[i];
    delete[] matrix;
}
void FluidDynamics::printDensToFile(int height, int width, std::ofstream &ofile){
  for(int i = 0; i<height; i++){
    for(int j=0; j<width; j++){
      ofile << std::setprecision(5) << rho[i][j] << endl;
    }
  }
}
