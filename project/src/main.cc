#include "Array.hh"
#include "FileReader.hh"
#include "Debug.hh"
#include <iostream>
#include <string>
#include "StaggeredGrid.hh"
#include "SORSolver.hh"
#include <cmath>
#include <functional>
#include "FluidSimulator.hh"



int main( int argc, char** argv )
{  
    std::string filename ;
   if(argc < 2)
   {
	   ERROR("Please pass a parameter file for SOR Solver");
   }
   else if(argc > 2)
   {
	   ERROR("Current implementation requires only one Parameter file for SOR Solver");
   }
   else
	   filename = argv[1];


   PRG_LEVEL("Start");
   FileReader fr;
   fr.readFile(filename);
   fr.printParameters();
   //StaggeredGrid grid(fr);
   FluidSimulator fs(fr);
   fs.simulateTimeStepCount();
   //fs.simulate(5);
   
   //grid.k_old().print();


   PRG_LEVEL("Finalize");
   return 0;
}
