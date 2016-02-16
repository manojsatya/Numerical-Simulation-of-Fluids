#include "FluidSimulator.hh"
#include "Debug.hh"


FluidSimulator::FluidSimulator( const FileReader & conf ):grid_(conf),solver_(conf),
	gx_{conf.getParameter<real>("gx")},
	gy_{conf.getParameter<real>("gy")},
	Re_{conf.getParameter<real>("Re")},
	gamma_{conf.getParameter<real>("gamma")},
	dt_{conf.getParameter<real>("dt")},
	timesteps_{conf.getParameter<size_t>("timesteps")}{}


StaggeredGrid &FluidSimulator::grid() {return grid_;}

const StaggeredGrid &FluidSimulator::grid() const{return grid_;}

void FluidSimulator::computeFG(){

	
	real idx = 1/real(grid_.dx()); 
	real idx2 = idx * idx ;
	real idy = 1/real(grid_.dy()); 
	real idy2 = idy * idy ;
 	real iRe = 1.0 / Re_ ;

	int imax = grid_.p().getSize(0)-2 ;
	int jmax = grid_.p().getSize(1)-2 ;


	Array u2_dx(imax+2,jmax+2);
	Array uv_dy(imax+2,jmax+2);
	Array u2_x2(imax+2,jmax+2);
	Array u2_y2(imax+2,jmax+2);
	Array p_dx(imax+2,jmax+2);

	Array uv_dx(imax+2,jmax+2);
	Array v2_dy(imax+2,jmax+2);
	Array v2_x2(imax+2,jmax+2);
	Array v2_y2(imax+2,jmax+2);
	Array p_dy(imax+2,jmax+2);

	
	for (int i = 1 ; i < imax ; ++i)
		for( int j = 1 ; j <= jmax ; ++j) {
	
	u2_dx(i,j) = 0.25 * idx * (((grid_.u()(i,j) + grid_.u()(i+1,j)) * (grid_.u()(i,j) + grid_.u()(i+1,j))) - ((grid_.u()(i-1,j) + grid_.u()(i,j)) * (grid_.u()(i-1,j) + grid_.u()(i,j))) + 
			gamma_ * (std::fabs ((grid_.u()(i,j) + grid_.u()(i+1,j))) * (grid_.u()(i,j) - grid_.u()(i+1,j)) - std::fabs((grid_.u()(i-1,j) + grid_.u()(i,j))) * (grid_.u()(i-1,j) - grid_.u()(i,j)))) ;

	
	uv_dy(i,j) = 0.25 * idy * ((grid_.v()(i,j) + grid_.v()(i+1,j)) * (grid_.u()(i,j) + grid_.u()(i,j+1)) - 
			(grid_.v()(i,j-1) + grid_.v()(i+1,j-1)) * (grid_.u()(i,j-1) + grid_.u()(i,j)) + 	gamma_ * (std::fabs((grid_.v()(i,j) + grid_.v()(i+1,j))) * (grid_.u()(i,j) - grid_.u()(i,j+1)) - 
			 std::fabs((grid_.v()(i,j-1) + grid_.v()(i+1,j-1))) * (grid_.u()(i,j-1) - grid_.u()(i,j)))) ; 

	

	u2_x2(i,j) = idx2 * (grid_.u()(i+1,j) - 2 * grid_.u()(i,j) + grid_.u()(i-1,j));

	//PROGRESS("***********calculating u **********");
	u2_y2(i,j) = idy2 * (grid_.u()(i,j+1) - 2 * grid_.u()(i,j) + grid_.u()(i,j-1));

	
	p_dx(i,j) = idx * (grid_.p()(i+1,j) - grid_.p()(i,j));
	
	

		} // End for loop u values

	
	
	for (int i = 1 ; i <= imax ; ++i)
		for( int j = 1 ; j < jmax ; ++j) {
	

	uv_dx(i,j) = 0.25 * idx * ((grid_.u()(i,j) + grid_.u()(i,j+1)) * (grid_.v()(i,j) + grid_.v()(i+1,j)) - (grid_.u()(i-1,j) + grid_.u()(i-1,j+1)) * (grid_.v()(i-1,j) + grid_.v()(i,j)) 	+ gamma_ * (std::fabs((grid_.u()(i,j) + grid_.u()(i,j+1))) * (grid_.v()(i,j) - grid_.v()(i+1,j)) - 		std::fabs((grid_.u()(i-1,j) + grid_.u()(i-1,j+1))) * (grid_.v()(i-1,j) + grid_.v()(i,j)))) ; 

	
	v2_dy(i,j) = 0.25 * idy * ((grid_.v()(i,j) + grid_.v()(i,j+1)) * (grid_.v()(i,j) + grid_.v()(i,j+1)) - 
			(grid_.v()(i,j-1) + grid_.v()(i,j)) * (grid_.v()(i,j-1) + grid_.v()(i,j)) + 
		gamma_ * ((grid_.v()(i,j) + grid_.v()(i,j+1)) * (grid_.v()(i,j) - grid_.v()(i,j+1)) - 
		(grid_.v()(i,j-1) + grid_.v()(i,j)) * (grid_.v()(i,j-1) - grid_.v()(i,j))));


	v2_x2(i,j) = idx2 * (grid_.v()(i+1,j) - 2 * grid_.v()(i,j) + grid_.v()(i-1,j)) ;

	v2_y2(i,j) = idy2 * (grid_.v()(i,j+1) - 2 * grid_.v()(i,j) + grid_.v()(i,j-1)) ;
	
	p_dy(i,j) = idy * (grid_.p()(i,j+1) - grid_.p()(i,j)) ;	



		} //end for loop v values

	
	for (int i = 1 ; i < imax ; ++i)
		for( int j = 1 ; j <= jmax ; ++j) {
	grid_.f()(i,j) = grid_.u()(i,j) + dt_ * ( iRe * ( u2_x2(i,j) + u2_y2(i,j) ) - u2_dx(i,j) - uv_dy(i,j) + gx_ ); 
	//PROGRESS("*********in f************"); 
	} // End for loop f values	


	for (int i = 1 ; i <= imax ; ++i)
		for( int j = 1 ; j < jmax ; ++j) {
	grid_.g()(i,j) = grid_.v()(i,j) + dt_ * ( iRe * ( v2_x2(i,j) + v2_y2(i,j) ) - v2_dy(i,j) - uv_dx(i,j) + gy_ );
	//PROGRESS("*********in g************"); 
	} // End for loop g values
//PROGRESS("*********end g************"); 
	
/*****************Boundary Values ****************/
	
	for ( int j = 1 ; j <= jmax ; ++j ){
		
	grid_.f()(0,j) = grid_.u()(0,j) ;
	
	grid_.f()(imax,j) = grid_.u()(imax,j);}

	for ( int i = 1 ; i <= imax ; ++i){
	grid_.g()(i,0) = grid_.v()(i,0);
	grid_.g()(i,jmax) = grid_.v()(i,jmax);}
	PROGRESS("*********************************************************");
	 PROGRESS("***********End F and G Calculations **********");
	 PROGRESS("*********************************************************");
/*************************************************/


}

void FluidSimulator::testSimulate(){computeFG();}
