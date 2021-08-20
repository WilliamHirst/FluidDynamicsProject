#include "FluidDynamics.hpp"
using namespace std;
//This main is used to calculate option values for spesific S and \tau values
//and compare them to analytical results.
int main(int argc, char* argv[])
{
   FluidDynamics FD;
   FD.Initialize(20, 20, 0.5, 0.1, 1.);
   FD.Lattice_Boltzmann(20,20);
   return 0;
}
