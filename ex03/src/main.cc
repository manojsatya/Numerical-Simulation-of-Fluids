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

 void testCase1( const FileReader & ifileread) {

	Array U(3,4), V(4,3), F(3,4), G(4,3);

	U(0,3) =  0.0;	U(1,3) =  1.0;	U(2,3) =  0.0;
	U(0,2) = -3.0;	U(1,2) = -1.0;	U(2,2) =  1.0;
	U(0,1) =  2.0;	U(1,1) =  0.0;	U(2,1) =  0.0;
	U(0,0) = -3.0;	U(1,0) =  0.0;	U(2,0) = -2.0;

	V(0,2) = -2.0;	V(1,2) =  1.0;	V(2,2) =  1.0;	V(3,2) = -1.0;
	V(0,1) = -1.0;	V(1,1) =  1.0;	V(2,1) =  0.0;	V(3,1) =  1.0;
	V(0,0) = -1.0;	V(1,0) =  0.0;	V(2,0) =  0.0;	V(3,0) =  2.0;

	F(0,3) =  0.0;	F(1,3) =  0.0;			F(2,3) =  0.0;
	F(0,2) = -3.0;	F(1,2) =  -0.99853;		F(2,2) =  1.0;
	F(0,1) =  2.0;	F(1,1) =  0.000853;		F(2,1) =  0.0;
	F(0,0) =  0.0;	F(1,0) =  0.0;			F(2,0) =  0.0;

	G(0,2) =  0.0;	G(1,2) =  1.0;		G(2,2) =  1.0;		G(3,2) =  0.0;
	G(0,1) =  0.0;	G(1,1) =  0.99912;	G(2,1) =  0.00021;	G(3,1) =  0.0;
	G(0,0) =  0.0;	G(1,0) =  0.0;		G(2,0) =  0.0;		G(3,0) =  0.0;

	std::cout << "Taken U: " << std::endl;
	U.print();
	std::cout << "Taken V: " << std::endl;
	V.print();

	FluidSimulator testSimulation ( ifileread);
	
	testSimulation.grid().u().setArray(U);
		
	testSimulation.grid().v().setArray(V);
	
	int im,jm;
	real a,b,tol;
	
	im = ifileread.getParameter<int>("imax");
	jm = ifileread.getParameter<int>("jmax");
	tol =0.001;	
	
	PROGRESS("************************** Test case************************************");
	testSimulation.testSimulate();
	
	std::cout << "Tolerance value taken for checking hand and computed values is :" << tol << std::endl ;
	
	for (int i = 1 ; i < im ; ++i)
		for( int j = 1 ; j <= jm ; ++j){
	
		  a = F(i,j);
		  b = testSimulation.grid().f()(i,j);
		   //std::cout << "a : " << a << " / b : " << b << std::endl;
		  //std ::cout << "absolute Value: "<<std::fabs(a-b) << std::endl ; 
		  CHECK_MSG( std::fabs(a-b) < tol ," F Values not matching . Please check ");
		  //std::cout << "F values are checked" << std::endl;
		}
		PROGRESS("*********F Values are matching ************");
	for (int i = 1 ; i <= im ; ++i)
		for( int j = 1 ; j < jm ; ++j){
	
		  a = G(i,j);
		  b = testSimulation.grid().g()(i,j);
		  //std::cout << "G(" << i << "," << j << ")" << std::endl;
		   //std::cout << "a : " << a << " / b : " << b << std::endl;
		  //std ::cout << "absolute Value: "<<std::fabs(a-b) << std::endl ; 
		  CHECK_MSG( std::fabs(a-b) < tol ," G Values not matching . Please check ");
		  //std::cout << "G values are checked" << std::endl;
		}
	
		PROGRESS("*********G Values are matching ************");
	PROGRESS("************************** End of Test case************************************");
	std::cout << "Hand Calculated F: " << std::endl;
	F.print();
	std::cout << "Computed F: " << std::endl;
	testSimulation.grid().f().print();
	std::cout << "\n Hand Calculated G: " << std::endl;
	G.print();
	std::cout << "Computed G:" << std::endl;
	testSimulation.grid().g().print();
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
	//ifileread.printParameters();

	FluidSimulator simulation(ifileread);
	//simulation.testSimulate();
//	PROGRESS("*********test simulate over************"); 
	
	testCase1(ifileread);
	
	
	
   		
	
   return 0;
}
