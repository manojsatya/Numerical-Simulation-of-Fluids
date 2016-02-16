#include "Array.hh"
#include "FileReader.hh"
#include "Debug.hh"
#include "StaggeredGrid.hh"
#include "SORSolver.hh"


#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>

#define PI 3.1415926535897932384626433832795028841971

 // to generate random init guess

void initGridSetup1(StaggeredGrid &grid){

	int imax = (grid.p()).getSize(0);
	int jmax = (grid.p()).getSize(1);

	for (int i = 0 ; i < imax ; ++i)
		for (int j = 0 ; j < jmax ; ++j)
			grid.p()(i,j) = rand() % 10 + 1 ;

	// RHS is already initialzed to zero. 
}

void initGridSetup2(StaggeredGrid &grid){

	int imax = (grid.p()).getSize(0);
	int jmax = (grid.p()).getSize(1);

	for (int i = 0 ; i < imax ; ++i)
		for (int j = 0 ; j < jmax ; ++j)
			grid.p()(i,j) = rand() % 10 + 1;

	real dx = grid.dx(); // Get the value of dx from staggered grid class

 	for (int i = 0 ; i < imax-1 ; ++i)
		for (int j = 0 ; j < jmax-1 ; ++j)

		grid.rhs()(i,j) = std::sin ( 2.0 * (i) * dx * PI ) ; // f(x,y) = sin(2xpi)

}


int main( int argc, char** argv )
{
	//Check to give only 2 argument in the command
	CHECK_MSG(argc == 2 ,"Please provide only one configuration file"); 
	std::string ifile(argv[1]);
	//int xSize,ySize,xcoord,ycoord;
	//real init;
	
	PROGRESS("******************NUSIF************");
	
	PROGRESS("******************Reading Input File************");
	std::cout << "Input File given is :"<<"	"<< ifile << std::endl;
	std::cout << "Reading parameters from file :"<<"	"<< ifile << std::endl;
	FileReader ifileread;
	/*************Reading parameters from file and outputting*******************/
	ifileread.readFile(ifile);

	
	
	StaggeredGrid grid(ifileread); // Grid is object created by Staggered Grid
	srand (time(NULL)); // Random generator
	SORSolver solver(ifileread); // SORSolver reads parameters from par file to solve SOR

	initGridSetup1 (grid); // Sending grid as argument for initGridSetup 	
	bool GridSetup1 = solver.solve(grid);// Solve the grid of initGridSetup1
	CHECK_MSG(GridSetup1 ,"Pressure poisson equation for Grid 1 is not solved");
	PROGRESS("******************Grid 1 Solved************");
	
	
	initGridSetup2 (grid); // Sending grid as argument for initGridSetup 	
	bool GridSetup2 = solver.solve(grid);// Solve the grid of initGridSetup2
	CHECK_MSG(GridSetup2 ,"Pressure poisson equation for Grid 2 is not solved");
	PROGRESS("******************Grid 2 Solved************");
	 	
   		
	
   return 0;
}
