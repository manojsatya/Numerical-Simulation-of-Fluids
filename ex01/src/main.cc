#include "Array.hh"
#include "FileReader.hh"
#include "Debug.hh"

#include <iostream>


int main( int argc, char** argv )
{
	//Check to give only 2 argument in the command
	CHECK_MSG(argc == 2 ,"Please provide only one configuration file"); 
	std::string ifile(argv[1]);
	int xSize,ySize,xcoord,ycoord;
	real init;
	
	PROGRESS("******************NUSIF************");
	
	PROGRESS("******************Reading Input File************");
	std::cout << "Input File given is :"<<"	"<< ifile << std::endl;
	std::cout << "Reading parameters from file :"<<"	"<< ifile << std::endl;
	FileReader ifileread;
	ifileread.readFile(ifile);
	//ifileread.printParameters();
	 	
   	xSize = ifileread.getParameter<int>("width");
	std::cout << "Value of width :" << xSize << std::endl ;
	ySize = ifileread.getParameter<int>("height");
	std::cout << "Value of height :" << ySize << std::endl ;
	xcoord = ifileread.getParameter<int>("x");
	std::cout << "Value of x-cord :" << xcoord << std::endl ;
	ycoord = ifileread.getParameter<int>("y");
	std::cout << "Value of y-cord :" << ycoord << std::endl ;
	init = ifileread.getParameter<real>("initial");	
	std::cout << "Value of initial :" << init << std::endl ;
	
	PROGRESS("******************Creating Array ************");
	Array dummy(xSize,ySize);
	dummy.fill(init);
	dummy.print();
	dummy(xcoord,ycoord) *= 2 ;
	std::cout << " After doubling the said x,y coordinates " << std::endl;
	dummy.print();	
	
   return 0;
}
