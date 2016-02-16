#include "SORSolver.hh"
#include "FileReader.hh"
#include <cmath>

SORSolver::SORSolver (size_t itermax ,const size_t & chkfreq, real omg, real eps ){
	itermax_ = itermax;
	chkfreq_ = chkfreq;
	omg_ = omg;
	eps_ = eps;
}

//Getting the parameter from the parsed file
SORSolver::SORSolver ( const FileReader & configuration ){

	itermax_ = configuration.getParameter<size_t>("itermax");
	chkfreq_ = configuration.getParameter<size_t>("checkfrequency");
	eps_ = configuration.getParameter<real>("eps");
	omg_ = configuration.getParameter<real>("omg");
	CHECK_MSG(omg_ < 2.0 ,"Omega value for SOR Solver should be less than 2"); 
	
}

bool SORSolver::solve( StaggeredGrid & grid ){

	CHECK_MSG(itermax_ > chkfreq_,"Maximum number of frequency should be more than check frequency"); 

// Initialgrid values are given in main function

//grid is object of Staggered Grid. Can access p() and rhs()  variables

int xbound{grid.p().getSize(0)-1};
int ybound{grid.p().getSize(1)-1};

	//std::cout << "Value of xbound:" << xbound << std::endl;
	//std::cout << "Value of ybound:" << ybound << std::endl;
	//std::cout << "Value of itermax:" << itermax_ << std::endl;
	//std::cout << "Value of eps:" << eps_ << std::endl;
	//std::cout << "Value of omg:" << omg_ << std::endl;
/************************************************************************************************/
/********************STEP 1 - Copying the values of the inner points to the boundary points*/
// Clockwise
//Left Boundary Treatment 


for(int j = 0;j <= ybound; ++j)
grid.p()(0,j) = grid.p()(1,j);//Left -> right ( j varies or y axis varies)
//PROGRESS("I am here ")
//Top Boundary Treatment 

for(int i = 0;i <= xbound; ++i)
grid.p()(i,ybound)=grid.p()(i,ybound-1);//down -> up ( i varies or x axis  )

//Right Boundary Treatment

for(int j = 0;j <= ybound; ++j)
grid.p()(xbound,j)=grid.p()(xbound-1,j); // Right -> left (j varies or y axis varies)

//Bottom boundary Treatment

for(int i = 0;i <= xbound; ++i)
grid.p()(i,0) = grid.p()(i,1); // up -> down ( i varies or x axis )

/**************************************************************************************************/ 

/********************STEP 2 - Performing one iteration SOR on the inner points*********/

	real dx = grid.dx();	real dx2 = dx * dx ; real invdx2 = 1.0 / dx2 ; 
	real dy = grid.dy();	real dy2 = dy * dy ; real invdy2 = 1.0 / dy2 ;
	real temp = omg_ * (1.0 / ((2.0 * invdx2) + (2.0 * invdy2)));
	real resNorm = 100.0 ; // To initialize residual much higher to check
    Array<real> res(xbound+1,ybound+1); // Res vector should have same size as p and rhs
	//Avoiding Division as much as possible because it takes too many cycles. PTFS Lessons !!!!
	// SOR Solver
	size_t SORiter = 0;
	//PROGRESS("Step 1");


	while ((resNorm >= eps_) && (SORiter < itermax_) ) {
	//for ( size_t SORiter = 1; SORiter <= itermax_ && resNorm >= eps_; ++SORiter) {
//	PROGRESS("Step 2");


	for(int i = 1 ; i < xbound ; ++i ) // Inner points only 
		for(int j = 1 ; j < ybound ; ++j) {
            if(grid.isFluid(i,j))
	// Equation 3 from Sheet02 . 
	grid.p()(i,j) = (1.0 - omg_)*grid.p()(i,j) + 
    temp * (invdx2 * (grid.p(i,j,EAST) + grid.p(i,j,WEST)) +
        invdy2 * (grid.p(i,j,NORTH) + grid.p(i,j,SOUTH)) -
		grid.rhs()(i,j)) ; }

/*********************************************************************************************/

/*****************************STEP 3 = REPEAT STEP 1 *****************************************/
//Left Boundary Treatment 
//PROGRESS("Step 3");
for(int j = 0;j <= ybound; ++j)
grid.p()(0,j) = grid.p()(1,j);//Left -> right ( j varies or y axis varies)

//Top Boundary Treatment 

for(int i = 0;i <= xbound; ++i)
grid.p()(i,ybound)=grid.p()(i,ybound-1);//down -> up ( i varies or x axis  )

//Right Boundary Treatment

for(int j = 0;j <= ybound; ++j)
grid.p()(xbound,j)=grid.p()(xbound-1,j); // Right -> left (j varies or y axis varies)

//Bottom boundary Treatment

for(int i = 0;i <= xbound; ++i)
grid.p()(i,0) = grid.p()(i,1); // up -> down ( i varies or x axis )

/*********************************************************************************************/
	if (SORiter%chkfreq_ == 0) {
/*****************************STEP 4 = Calculate Residual *****************************************/
//PROGRESS("Step 4");
for(int i = 1 ; i < xbound ; ++i ) // Inner points only 
        for(int j = 1 ; j < ybound ; ++j)
            if(grid.isFluid(i,j))
	// Residual = f{h} - A * p{h}
    res(i,j) = grid.rhs()(i,j) - invdx2 * (grid.p(i,j,EAST) - 2.0 * grid.p(i,j,CENTER) + grid.p(i,j,WEST))
                - invdy2 * (grid.p(i,j,NORTH) - 2.0 * grid.p(i,j,CENTER) + grid.p(i,j,SOUTH));

	// Calculating Norm of Residual 

	resNorm = 0.0 ; // Before calculating norm at every iteration,initializing to 0.0

	for ( int i= 0 ; i <= xbound ; ++i)
		for(int j = 0 ; j <= ybound ; ++j)
		  if(grid.isFluid(i,j))
			resNorm += res(i,j) * res(i,j); // square of res from Equation 4 from sheet02  
			
	
	//resNorm = sqrt(resNorm/((ybound+1)*(xbound+1)));
			resNorm = sqrt(resNorm/(grid.getNumFluid()));
	//std::cout << "Iteration:" << SORiter << "->	" << resNorm << std::endl;
	  
	}  // Square root as mentioned in eq 4
	
	 
	 
/**************************************************************************************************/

/************************STEP 5 = Residual check and repeat from STEP 2 ********************/
		++ SORiter ; 
	
	}// Close while loop
	std::cout << "Residual - > "  << resNorm << std::endl;
	 //std::cout << "Number of iterations it took: " << SORiter << " with residual norm ->	" << resNorm << std::endl;
	return resNorm  ; //<= eps_; // If residual is good enough  

} // End of bool solver 
