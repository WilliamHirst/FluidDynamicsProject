#include "FluidDynamics.hpp"
#include <cmath>
#include <stdlib.h>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <fstream>
using namespace std;

void FluidDynamics::Initialize(int height, int width, double inputOmega, double inputInitialVel, double fTime)
{
  nu = 1./3. * (1./omega - 0.5);
  c = 1;
  finalTime = fTime;
  dt = (double) 1. / 30.;
  three_Half = (double) 3 / 2;
  nine_Half = (double) 9 / 2;
  nZero = createMatrix(height, width);
  nOne = createMatrixOne(height, width);
  nTwo = createMatrix(height, width);
  nThree = createMatrix(height, width);
  nFour = createMatrix(height, width);
  nFive = createMatrix(height, width);
  nSix = createMatrix(height, width);
  nSeven = createMatrix(height, width);
  nEight = createMatrix(height, width);
  nNine = createMatrix(height, width);
  u_x = createMatrix(height, width);
  u_y = createMatrix(height, width);
  rho = createMatrix(height, width);
  Omega = new double[10];
  unit_x = new double[10];
  unit_y = new double[10];
  equiDens = new double[10];
  initial_vel = inputInitialVel;
  initialSquared = initial_vel*initial_vel;
  omega = inputOmega;
  Omega[0] = (double) 4. / 9.;
  Omega[1] = Omega[2] = Omega[3] = Omega[4] = (double) 1. / 9.;
  Omega[5] = Omega[6] = Omega[7] = Omega[8] = (double) 1. / 36.;
  unit_x[0] = unit_x[2] = unit_x[4] = unit_y[0] = unit_y[1] = unit_y[3] = 0;
  unit_x[1] = unit_x[5] = unit_x[6] = unit_y[2] = unit_y[5] = unit_y[6] = 1;
  unit_x[3] = unit_x[6] = unit_x[7] = unit_y[4] = unit_y[7] = unit_y[8] = -1;
  //Initialize rho, u_x and u_y.
  for(int i = 0; i < height; i++){
    for(int j = 0; j < width; j++) {
      nZero[i][j] = Omega[0]*(1. - 1.5*initialSquared);
      nOne[i][j] = Omega[1]*(1. - 1.5*initialSquared);
      nTwo[i][j] = Omega[2]*(1. - 1.5*initialSquared);
      nThree[i][j] = Omega[3]*(1. +  3*initial_vel + 4.5*initialSquared - 1.5*initialSquared);
      nFour[i][j] = Omega[4]*(1. -  3*initial_vel + 4.5*initialSquared - 1.5*initialSquared);
      nFive[i][j] = Omega[5]*(1. +  3*initial_vel + 4.5*initialSquared - 1.5*initialSquared);
      nSix[i][j] = Omega[5]*(1. +  3*initial_vel + 4.5*initialSquared - 1.5*initialSquared);
      nSeven[i][j] = Omega[6]*(1. +  3*initial_vel + 4.5*initialSquared - 1.5*initialSquared);
      nEight[i][j] = Omega[7]*(1. +  3*initial_vel + 4.5*initialSquared - 1.5*initialSquared);
      rho[i][j]  = (nZero[i][j] + nOne[i][j] + nTwo[i][j] + nThree[i][j] + nFour[i][j]\
      + nFive[i][j] + nSix[i][j] + nSeven[i][j] + nEight[i][j]);
      u_x[i][j] = (nOne[i][j] + nThree[i][j] + nFive[i][j] + nSix[i][j] + nSeven[i][j] + nEight[i][j])/(rho[i][j]);
      u_y[i][j] = (nTwo[i][j] + nFour[i][j] + nFive[i][j] + nSix[i][j] + nSeven[i][j] + nEight[i][j])/(rho[i][j]);
    }
  }
}
void FluidDynamics::find_density(int height, int width)
{
  rho[height][width]  = (nZero[height][width] + nOne[height][width] + nTwo[height][width] + nThree[height][width] + nFour[height][width]\
  + nFive[height][width] + nSix[height][width] + nSeven[height][width] + nEight[height][width]);
  u_x[height][width] = (nOne[height][width] + nThree[height][width] + nFive[height][width] + nSix[height][width] + nSeven[height][width] + nEight[height][width])/(rho[height][width]);
  u_y[height][width] = (nTwo[height][width] + nFour[height][width] + nFive[height][width] + nSix[height][width] + nSeven[height][width] + nEight[height][width])/(rho[height][width]);
  double ux2 = u_x[height][width]*u_x[height][width];
	double uy2 = u_y[height][width]*u_y[height][width];
  double u2 = ux2 + uy2;
  double uxuy = u_x[height][width]*u_y[height][width];
  double omu215 = 1 - 1.5*u2;
  nZero[height][width] = (1-omega)*nZero[height][width] + omega * Omega[0] * rho[height][width] * omu215;
	nOne[height][width] = (1-omega)*nOne[height][width] + omega * Omega[1] * rho[height][width] * (omu215 + 3*u_y[height][width] + 4.5*uy2);
	nTwo[height][width] = (1-omega)*nTwo[height][width] + omega * Omega[2] * rho[height][width] * (omu215 - 3*u_y[height][width] + 4.5*uy2);
	nThree[height][width] = (1-omega)*nThree[height][width] + omega * Omega[3] * rho[height][width] * (omu215 + 3*u_x[height][width] + 4.5*ux2);
	nFour[height][width] = (1-omega)*nFour[height][width] + omega * Omega[4] * rho[height][width] * (omu215 - 3*u_x[height][width] + 4.5*ux2);
	nFive[height][width] = (1-omega)*nFive[height][width] + omega * Omega[5] * rho[height][width] * (omu215 + 3*(u_x[height][width]+u_y[height][width]) + 4.5*(u2+2*uxuy));
	nSix[height][width] = (1-omega)*nSix[height][width] + omega * Omega[6] * rho[height][width] * (omu215 + 3*(-u_x[height][width]+u_y[height][width]) + 4.5*(u2-2*uxuy));
	nSeven[height][width] = (1-omega)*nSeven[height][width] + omega * Omega[7] * rho[height][width] * (omu215 + 3*(u_x[height][width]-u_y[height][width]) + 4.5*(u2-2*uxuy));
	nEight[height][width] = (1-omega)*nEight[height][width] + omega * Omega[8] * rho[height][width] * (omu215 + 3*(-u_x[height][width]-u_y[height][width]) + 4.5*(u2+2*uxuy));

  if (height == 0) {
    nZero[height][width] = Omega[0]*(1. - 1.5*initialSquared);
    nOne[height][width] = Omega[1]*(1. - 1.5*initialSquared);
    nTwo[height][width] = Omega[2]*(1. - 1.5*initialSquared);
    nThree[height][width] = Omega[3]*(1. +  3*initial_vel + 4.5*initialSquared - 1.5*initialSquared);
    nFour[height][width] = Omega[4]*(1. -  3*initial_vel + 4.5*initialSquared - 1.5*initialSquared);
    nFive[height][width] = Omega[5]*(1. +  3*initial_vel + 4.5*initialSquared - 1.5*initialSquared);
    nSix[height][width] = Omega[5]*(1. +  3*initial_vel + 4.5*initialSquared - 1.5*initialSquared);
    nSeven[height][width] = Omega[6]*(1. +  3*initial_vel + 4.5*initialSquared - 1.5*initialSquared);
    nEight[height][width] = Omega[7]*(1. +  3*initial_vel + 4.5*initialSquared - 1.5*initialSquared);
  }
}
void FluidDynamics::Barrier(int height, int width){
  int j = width/4;
  for(int i = height/4; i< height*3/4; i++){
    //Bounce particles back.
    nOne[i][j-1] += nOne[i][j];
    nTwo[i-1][j] += nTwo[i][j];
    nThree[i-1][j] += nThree[i][j];
    nFour[i+1][j] += nFour[i][j] ;
    nFive[i-1][j+1] += nFive[i][j];
    nSix[i+1][j+1] += nSix[i][j];
    nSeven[i+1][j-1] = nSeven[i][j];
    nEight[i-1][j-1] = nEight[i][j];
    //Set all particles in barrier to 0.
    nOne[i][j] = 0;
    nTwo[i][j] = 0;
    nThree[i][j] = 0;
    nFour[i][j] = 0;
    nFive[i][j] = 0;
    nSix[i][j] = 0;
    nSeven[i][j] = 0;
    nEight[i][j] = 0;
  }
}
void FluidDynamics::Roll(int height, int width) {

  for(int i = 1; i<height-1; i++){
    for(int j=1; j<width-1; j++){
      nOne[i][j] = nOne[i][j+1];
      nTwo[i][j] = nTwo[i+1][j];
      nThree[i][j] = nThree[i+1][j];
      nFour[i][j] = nFour[i-1][j];
      nFive[i][j] = nFive[i+1][j-1];
      nSix[i][j] = nSix[i-1][j-1];
      nSeven[i][j] = nSeven[i-1][j+1];
      nEight[i][j] = nEight[i+1][j+1];
    }
  }
  for(int i = 0; i<height; i++) {
    nOne[i][width-1] = nOne[i][0];
    nFive[i][width-1] = nFive[i][0];
    nEight[i][width-1] = nEight[i][0];
  }
}
void FluidDynamics::Lattice_Boltzmann(int height, int width)
{
  ofstream ofile;
  string outfilename = "particleDensity.txt";
  ofile.open(outfilename);
  ofile << setprecision(16)<< width << " "<< height <<endl;
  double counter = 0.;
  while(counter<finalTime) {
    for(int i = 0; i<height; i++){
      for(int j=0; j<width; j++){
        Roll(height, width);
        find_density(i,j);
        Barrier(height, width);
        printDensToFile(height, width, ofile);
      }
    }
    counter += dt;
  }
  ofile.close();
}

double ** FluidDynamics::createMatrix(int height, int width) {
  // Set up matrix
  double **matrix;
  matrix = new double*[height];
  // Allocate memory
  for(int i=0; i < height ;i++)
      matrix[i] = new double[width];
  // Set values to zero
  return matrix;
}
double ** FluidDynamics::createMatrixOne(int height, int width) {
  // Set up matrix
  double **matrix;
  matrix = new double*[height];
  // Allocate memory
  for(int i=0; i < height ;i++)
      matrix[i] = new double[width];
  // Set values to zero
  for(int i = 0; i < height; i++){
      for(int j = 0; j < width; j++){
          matrix[i][j] = 1.0;
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
      ofile << setprecision(5) << rho[i][j] << " ";
    }
    ofile << setprecision(5) << " " << endl;
  }

}
