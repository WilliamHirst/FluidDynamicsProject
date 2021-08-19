#include "FluidDynamics.hpp"
using namespace arma;
using namespace std;
//This main is used to calculate option values for spesific S and \tau values
//and compare them to analytical results.
int main(int argc, char* argv[])
{
   FluidDynamics FD;
   FD.Initialize(200, 100, 0.5, 0.4, 5);
   FD.Lattice_Boltzmann(200,100);
   return 0;
}
