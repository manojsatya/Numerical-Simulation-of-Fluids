#include "Array.hh"
#include "FileReader.hh"
#include "Debug.hh"
//#include "StaggeredGrid.hh"
#include "SORSolver.hh"
#include "FluidSimulator.hh"
#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>

#define PI 3.1415926535897932384626433832795028841971

 
int main( int argc, char** argv )
{
	//Check to give only 2 argument in the command
	CHECK_MSG(argc == 2 ,"Please provide only one configuration file"); 
	std::string ifile(argv[1]);
	
	
	PROGRESS("******************NUSIF************");
	
	PROGRESS("******************Reading Input File************");
	std::cout << "Input File given is :"<<"	"<< ifile << std::endl;
	std::cout << "Reading parameters from file :"<<"	"<< ifile << std::endl;
	
	FileReader ifileread;


	/*************Reading parameters from file and outputting*******************/
	ifileread.readFile(ifile);
	ifileread.printParameters();

	size_t timestepcount = ifileread.getParameter<size_t>("timesteps");
	//real duration = 1.0 ; // input in seconds

	FluidSimulator simulation(ifileread);
	

	simulation.simulateTimeStepCount(timestepcount);
	//simulation.simulate(duration);
	
	// To get values from specific point in the domain . Inputs are in (x,y) co-ordinates
	simulation.getU(0.7,0.7); // U velocity value
	simulation.getV(0.7,0.7); // V velocity value
	simulation.getP(0.7,0.7); // Pressure value
	
	
	
   		
	
   return 0;
}
