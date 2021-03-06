#include "Array.hh"
#include <vector>



//===================================================================================================================
//
//  Constructors
//
//===================================================================================================================


Array::Array( int xSize )
{
	CHECK_MSG(xSize > 1, "Width of simulation data has to be more than 1");
	matrix.erase(matrix.begin(),matrix.end()); // Erase vector 
	matrix.resize(xSize,0.0); // Resizing vector and with initializing with 0.0
	imax = xSize ; jmax = 1 ; kmax = 1 ;
   // TODO construct 1D array here
}

Array::Array( int xSize, int ySize )
{
	CHECK_MSG(xSize > 1 && ySize > 1, "Width and height of 2-D Simulation data has to be more than 1 * 1");
	matrix.erase(matrix.begin(),matrix.end()); // Erase vector
	matrix.resize(xSize * ySize,0.0);// Resizing vector and with initializing with 0.0
	imax = xSize ; jmax = ySize ; kmax = 1 ;
   // TODO construct 2D array here
}

Array::Array( int xSize, int ySize, int zSize )
{
	CHECK_MSG(xSize > 1 && ySize > 1 && zSize, "Width,height and depth of 3-D Simulation data has to be more than 1 * 1 * 1");
	matrix.erase(matrix.begin(),matrix.end()); // Erase vector
	matrix.resize(xSize * ySize * zSize ,0.0); // Resizing vector and with initializing with 0.0
	imax = xSize ; jmax = ySize ; kmax = zSize ;
   // TODO construct 3D array here
}




//===================================================================================================================
//
//  Convenience Functions
//
//===================================================================================================================


//initialize the whole array with a constant value
void Array::fill( real value )
{
   // TODO
	std::fill(matrix.begin(),matrix.end(),value);
   // you might want to use std::fill() here
}


// Print the whole array (for debugging purposes)
void Array::print()
{
   // TODO
   // For 2D Arrays the positive x-coordinate goes to the right
   //                   positive y-coordinate goes upwards
   //      -> the line with highest y-value should be printed first
	if ( imax != 1 && jmax != 1 && kmax == 1){
		for(int i = jmax-1 ; i>=0 ; --i){
			for(int j = 0 ; j < imax ; ++j )
				std::cout << matrix[j + imax*i] << " \t";
			std::cout << std::endl;
	} //end of for loops
 } // end of if 

	// For loop to print 3D Array 
	else for (std::vector<real>::iterator it = matrix.begin() ; it!=matrix.end() ; ++it )
		std::cout << *it << "\t " ;

}

int Array::getSize( int dimension ) const
{
ASSERT_MSG(dimension !=0 || dimension != 1 || dimension != 2, "****Invalid dimension Entry******");
	if(dimension == 0) return imax;
	else if (dimension == 1) return jmax;
	else return kmax ;
   //TODO
   //return 0;
}

//return total size of the array
int Array::getSize() const
{
   //TODO
   return imax*jmax*kmax;
}

//bool Array::maxelem(const real& a , const real& b){return std::fabs(a) < std::fabs(b);}	

real Array::get_abs_maxelement(){

	//bool maxelem(const real& a , const real& b){return std::fabs(a) < std::fabs(b);}
	//real abs_Maxelement = std::fabs(*std::max_element(matrix.begin(),matrix.end(),myelem);	
	real abs_Maxelement = std::fabs(*std::max_element(matrix.begin(),matrix.end(),[](const real& a , const real& b){return std::fabs(a) < std::fabs(b);}));
	return abs_Maxelement;
}