#include <precice/SolverInterface.hpp>
#include <iostream>

int main(int argc, char** argv)
{
   /* Initialize the preCICE participant for the fluid solver */
   precice::SolverInterface participant("FluidSolver", "precice-config.xml", 0, 1);

   /* Get a reference to the coupling mesh */
   int meshID = participant.getMeshID("MeshFSI");


  return 0;
}