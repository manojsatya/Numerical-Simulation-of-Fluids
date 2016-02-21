#include "FluidSimulator.hh"
#include "Debug.hh"




FluidSimulator::FluidSimulator( const FileReader & conf ):grid_(conf),solver_(conf),
	gx_{conf.getParameter<real>("GX")},
	gy_{conf.getParameter<real>("GY")},
	Re_{conf.getParameter<real>("Re")},
        nu{conf.getParameter<real>("nu")},
        turbulenceMode{conf.getParameter<int>("Turbulence")},
	gamma_{conf.getParameter<real>("gamma")},
	dt_{conf.getParameter<real>("dt")},
	timesteps_{conf.getParameter<size_t>("timesteps")},
	tau_{conf.getParameter<real>("safetyfactor")},
	nf_{conf.getParameter<size_t>("normalizationfrequency")},
	xlength_{conf.getParameter<real>("xlength")},
	ylength_{conf.getParameter<real>("ylength")},
	outptinter_{conf.getParameter<size_t>("outputinterval")}
	
	{
	CHECK_MSG(Re_ > 0.0 ," Reynolds number cannot be negative ");
	CHECK_MSG(gamma_ > 0.0 && gamma_ < 1.0 ," Gamma should lie between 0 and 1 ");
	name_ = conf.getStringParameterBoundary("name");
	
	
	if (conf.checkParameter("boundary_condition_N")){
	//CHECK_MSG(boundary_condition_N_ != "inflow" || !(conf.getrealParameterBoundary("boundary_velocity_N")) , "North Boundary cannot have two boundary conditions") ;
	CHECK_MSG(boundary_condition_N_ != "outflow" || !(conf.getrealParameterBoundary("boundary_velocity_N")) , "North Boundary cannot have two boundary conditions") ;}
	
	if (conf.checkParameter("boundary_condition_S")){
	//CHECK_MSG(boundary_condition_S_ != "inflow" || !(conf.getrealParameterBoundary("boundary_velocity_S")) , "North Boundary cannot have two boundary conditions") ;
	CHECK_MSG(boundary_condition_S_ != "outflow" || !(conf.getrealParameterBoundary("boundary_velocity_S")) , "North Boundary cannot have two boundary conditions") ;}

	if (conf.checkParameter("boundary_condition_E")){
	//CHECK_MSG(boundary_condition_E_ != "inflow" || !(conf.getrealParameterBoundary("boundary_velocity_E")) , "North Boundary cannot have two boundary conditions") ;
	CHECK_MSG(boundary_condition_E_ != "outflow" || !(conf.getrealParameterBoundary("boundary_velocity_E")) , "North Boundary cannot have two boundary conditions") ;}
	
	if (conf.checkParameter("boundary_condition_W")){
	//CHECK_MSG(boundary_condition_W_ != "inflow" || !(conf.getrealParameterBoundary("boundary_velocity_W")) , "North Boundary cannot have two boundary conditions") ;
	CHECK_MSG(boundary_condition_W_ != "outflow" || !(conf.getrealParameterBoundary("boundary_velocity_W")) , "North Boundary cannot have two boundary conditions") ;}
	
	boundary_condition_S_ = conf.getStringParameterBoundary("boundary_condition_S");
	boundary_condition_N_ = conf.getStringParameterBoundary("boundary_condition_N");
	boundary_condition_W_ = conf.getStringParameterBoundary("boundary_condition_W");
	boundary_condition_E_ = conf.getStringParameterBoundary("boundary_condition_E");		
	boundary_velocity_S_ = conf.getrealParameterBoundary("boundary_velocity_S");
	boundary_velocity_N_ = conf.getrealParameterBoundary("boundary_velocity_N");
	boundary_velocity_W_ = conf.getrealParameterBoundary("boundary_velocity_W");
	boundary_velocity_E_ = conf.getrealParameterBoundary("boundary_velocity_E");
	
	if ( boundary_condition_S_ == "periodic" ) 
	  boundary_condition_N_ == "periodic" ; 
	else if (boundary_condition_N_ == "periodic")
	  boundary_condition_S_ == "periodic";

	 if ( boundary_condition_E_ == "periodic" ) 
	  boundary_condition_W_ == "periodic" ; 
	else if (boundary_condition_W_ == "periodic")
	  boundary_condition_E_ == "periodic";

         //if(conf.checkParameter("Turbulence")) {turbulenceMode = 0;PROGRESS(" Turbulence Mode set to Laminar");}
         //if(conf.checkParameter("nu")){nu = 0;}

         c_nu = 0.09;
         c_eps = 0.07; // Turbulence parameters
         c_1 = 0.126;
         c_2 = 1.92;
}


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


    Array<real> u2_dx(imax+2,jmax+2);
    Array<real> uv_dy(imax+2,jmax+2);
    Array<real> u2_x2(imax+2,jmax+2);
    Array<real> u2_y2(imax+2,jmax+2);
	//Array p_dx(imax+2,jmax+2);

    Array<real> uv_dx(imax+2,jmax+2);
    Array<real> v2_dy(imax+2,jmax+2);
    Array<real> v2_x2(imax+2,jmax+2);
    Array<real> v2_y2(imax+2,jmax+2);
	//Array p_dy(imax+2,jmax+2);

	
	for (int i = 1 ; i < imax ; ++i)
		for( int j = 1 ; j <= jmax ; ++j) {

                    if(grid_.isFluid(i,j) && grid_.isFluid(i+1,j)){
	 
	//u2_dx(i,j) = 0.25 * idx * (((grid_.u()(i,j) + grid_.u()(i+1,j)) * (grid_.u()(i,j) + grid_.u()(i+1,j))) - ((grid_.u()(i-1,j) + grid_.u()(i,j)) * (grid_.u()(i-1,j) + grid_.u()(i,j))) + 
	//		gamma_ * (std::fabs ((grid_.u()(i,j) + grid_.u()(i+1,j))) * (grid_.u()(i,j) - grid_.u()(i+1,j)) - std::fabs((grid_.u()(i-1,j) + grid_.u()(i,j))) * (grid_.u()(i-1,j) - grid_.u()(i,j)))) ;

	//u2_dx(i,j) = 0.25 * idx * ((grid_.u()(i,j) + grid_.u()(i+1,j)) * (grid_.u()(i,j) + grid_.u()(i+1,j)) - (grid_.u()(i-1,j) + grid_.u()(i,j) * grid_.u()(i-1,j) + grid_.u()(i,j))
	//	      + gamma_* (std::fabs(grid_.u()(i,j) + grid_.u()(i+1,j)) * (grid_.u()(i,j)-grid_.u()(i,j+1)) - std::fabs(grid_.u()(i-1,j) + grid_.u()(i,j)) * (grid_.u()(i-1,j) - grid_.u()(i,j))));
		  
    u2_dx(i,j) = 0.25 * idx * (((grid_.u(i,j,CENTER) + grid_.u(i,j,EAST))	*(grid_.u(i,j,CENTER)+ grid_.u(i,j,EAST))-(grid_.u(i,j,WEST)+ grid_.u(i,j,CENTER))*(grid_.u(i,j,WEST)+ grid_.u(i,j,CENTER)))
            + gamma_*(std::fabs(grid_.u(i,j,CENTER)+ grid_.u(i,j,EAST))*(grid_.u(i,j,CENTER)- grid_.u(i,j,EAST))- std::fabs(grid_.u(i,j,WEST)+ grid_.u(i,j,CENTER))*(grid_.u(i,j,WEST)- grid_.u(i,j,CENTER))));
		  
		//uv_dy(i,j) = 0.25 * idy * ((grid_.v()(i,j) + grid_.v()(i+1,j)) * (grid_.u()(i,j) + grid_.u()(i,j+1)) - 
	//		(grid_.v()(i,j-1) + grid_.v()(i+1,j-1)) * (grid_.u()(i,j-1) + grid_.u()(i,j)) + 	gamma_ * (std::fabs((grid_.v()(i,j) + grid_.v()(i+1,j))) * (grid_.u()(i,j) - grid_.u()(i,j+1)) - 
	//		 std::fabs((grid_.v()(i,j-1) + grid_.v()(i+1,j-1))) * (grid_.u()(i,j-1) - grid_.u()(i,j)))) ; 
       
    uv_dy(i,j) = 0.25 * idy * ( (grid_.v(i,j,CENTER) + grid_.v(i,j,EAST)) * (grid_.u(i,j,CENTER) + grid_.u(i,j,NORTH)) - (grid_.v(i,j,SOUTH) + grid_.v(i+1,j,SOUTH)) * (grid_.u(i,j,SOUTH) + grid_.u(i,j,CENTER))
              + gamma_ * ( std::fabs(grid_.v(i,j,CENTER) + grid_.v(i,j,EAST)) * (grid_.u(i,j,CENTER) - grid_.u(i,j,NORTH)) - std::fabs(grid_.v(i,j,SOUTH) + grid_.v(i+1,j,SOUTH)) * (grid_.u(i,j,SOUTH) - grid_.u(i,j,CENTER))));
    
    u2_x2(i,j) = idx2 * (grid_.u(i,j,EAST) - 2 * grid_.u(i,j,CENTER) + grid_.u(i,j,WEST));

	//PROGRESS("***********calculating u **********");
    u2_y2(i,j) = idy2 * (grid_.u(i,j,NORTH) - 2 * grid_.u(i,j,CENTER) + grid_.u(i,j,SOUTH));

	
	//p_dx(i,j) = idx * (grid_.p()(i+1,j) - grid_.p()(i,j));
	
                        }

		} // End for loop u values

	
	
	for (int i = 1 ; i <= imax ; ++i)
		for( int j = 1 ; j < jmax ; ++j) {
	
         if(grid_.isFluid(i,j) && grid_.isFluid(i,j+1)) {
	//uv_dx(i,j) = 0.25 * idx * ((grid_.u()(i,j) + grid_.u()(i,j+1)) * (grid_.v()(i,j) + grid_.v()(i+1,j)) - (grid_.u()(i-1,j) + grid_.u()(i-1,j+1)) * (grid_.v()(i-1,j) + grid_.v()(i,j)) 	+ gamma_ * (std::fabs((grid_.u()(i,j) + grid_.u()(i,j+1))) * (grid_.v()(i,j) - grid_.v()(i+1,j)) - 		std::fabs((grid_.u()(i-1,j) + grid_.u()(i-1,j+1))) * (grid_.v()(i-1,j) + grid_.v()(i,j)))) ; 

		  
		  
		  uv_dx(i,j) = 0.25 * idx * ((grid_.u(i,j,CENTER) + grid_.u(i,j,NORTH)) * (grid_.v(i,j,CENTER) + grid_.v(i,j,EAST)) - (grid_.u(i,j,WEST) + grid_.u(i,j+1,WEST)) * (grid_.v(i,j,WEST) + grid_.v(i,j,CENTER))
             + gamma_ * (std::fabs(grid_.u(i,j,CENTER) + grid_.u(i,j,NORTH)) * (grid_.v(i,j,CENTER) - grid_.v(i,j,EAST)) - std::fabs(grid_.u(i,j,WEST) + grid_.u(i,j+1,WEST)) * (grid_.v(i,j,WEST)-grid_.v(i,j,CENTER))));
		//PROGRESS("uv_dx");  
	//v2_dy(i,j) = 0.25 * idy * ((grid_.v()(i,j) + grid_.v()(i,j+1)) * (grid_.v()(i,j) + grid_.v()(i,j+1)) - 
	//		(grid_.v()(i,j-1) + grid_.v()(i,j)) * (grid_.v()(i,j-1) + grid_.v()(i,j)) + 
	//	gamma_ * ((grid_.v()(i,j) + grid_.v()(i,j+1)) * (grid_.v()(i,j) - grid_.v()(i,j+1)) - 
	//	(grid_.v()(i,j-1) + grid_.v()(i,j)) * (grid_.v()(i,j-1) - grid_.v()(i,j))));
	
    v2_dy(i,j) = 0.25 * idy * ((grid_.v(i,j,CENTER) + grid_.v(i,j,NORTH)) * (grid_.v(i,j,CENTER) + grid_.v(i,j,NORTH)) - (grid_.v(i,j,CENTER) + grid_.v(i,j,SOUTH)) * (grid_.v(i,j,CENTER) + grid_.v(i,j,SOUTH))
            + gamma_ * (std::fabs(grid_.v(i,j,CENTER) + grid_.v(i,j,NORTH)) * (grid_.v(i,j,CENTER) - grid_.v(i,j,NORTH)) - std::fabs(grid_.v(i,j,CENTER) + grid_.v(i,j,SOUTH)) * (grid_.v(i,j,SOUTH) - grid_.v(i,j,CENTER))));
//PROGRESS("v2_dy"); 
    v2_x2(i,j) = idx2 * (grid_.v(i,j,EAST) - 2 * grid_.v(i,j,CENTER) + grid_.v(i,j,WEST)) ;
//PROGRESS("v2_x2"); 
    v2_y2(i,j) = idy2 * (grid_.v(i,j,NORTH) - 2 * grid_.v(i,j,CENTER) + grid_.v(i,j,SOUTH)) ;
//	PROGRESS("v2_y2"); 
        //p_dy(i,j) = idy * (grid_.p()(i,j+1) - grid_.p()(i,j)) ;
                        }
		} //end for loop v values

	//PROGRESS("*********in f************"); 
	for (int i = 1 ; i < imax ; ++i)
        for( int j = 1 ; j <= jmax ; ++j)
            if(grid_.isFluid(i,j) && grid_.isFluid(i+1,j))
    grid_.f()(i,j) = grid_.u(i,j,CENTER) + dt_ * ( iRe * ( u2_x2(i,j) + u2_y2(i,j) ) - u2_dx(i,j) - uv_dy(i,j) + gx_ );
	
     // End for loop f values


	for (int i = 1 ; i <= imax ; ++i)
        for( int j = 1 ; j < jmax ; ++j)
            if(grid_.isFluid(i,j) && grid_.isFluid(i,j+1))
    grid_.g()(i,j) = grid_.v(i,j,CENTER) + dt_ * ( iRe * ( v2_x2(i,j) + v2_y2(i,j) ) - uv_dx(i,j) - v2_dy(i,j)  + gy_ );
	//PROGRESS("*********in g************"); 
     // End for loop g values
//PROGRESS("*********end g************"); 
	
/*****************Boundary Values ****************/
	
	for ( int j = 1 ; j <= jmax ; ++j ){
		
	grid_.f()(0,j) = grid_.u()(0,j) ;
	
	grid_.f()(imax,j) = grid_.u()(imax,j);}

	for ( int i = 1 ; i <= imax ; ++i){
	grid_.g()(i,0) = grid_.v()(i,0);
	grid_.g()(i,jmax) = grid_.v()(i,jmax);}
	 //PROGRESS("***********End f and g **********");
/*************************************************/


}

void FluidSimulator::computeFG_KE(){


        real idx = 1/real(grid_.dx());
        real idx2 = idx * idx ;
        real idy = 1/real(grid_.dy());
        real idy2 = idy * idy ;
        //real iRe = 1.0 / Re_ ;

        int imax = grid_.p().getSize(0)-2 ;
        int jmax = grid_.p().getSize(1)-2 ;

    Array<real> first_f(imax+2,jmax+2);
    Array<real> second_f(imax+2,jmax+2);
    Array<real> third_f(imax+2,jmax+2);

    Array<real> u2_dx(imax+2,jmax+2);
    Array<real> uv_dy(imax+2,jmax+2);
    //Array<real> u2_x2(imax+2,jmax+2);
    //Array<real> u2_y2(imax+2,jmax+2);
        //Array p_dx(imax+2,jmax+2);

    Array<real> first_g(imax+2,jmax+2);
    Array<real> second_g(imax+2,jmax+2);
    Array<real> third_g(imax+2,jmax+2);
    Array<real> uv_dx(imax+2,jmax+2);
    Array<real> v2_dy(imax+2,jmax+2);
    //Array<real> v2_x2(imax+2,jmax+2);
    //Array<real> v2_y2(imax+2,jmax+2);
        //Array p_dy(imax+2,jmax+2);


        for (int i = 1 ; i < imax ; ++i)
                for( int j = 1 ; j <= jmax ; ++j) {

                    if(grid_.isFluid(i,j) && grid_.isFluid(i+1,j)){

      first_f(i,j) = idx2* (nu_t_str(i+1,j)*(grid_.u(i,j,EAST) - grid_.u(i,j,CENTER)) - nu_t_str(i,j) *(grid_.u(i,j,CENTER) - grid_.u(i,j,WEST)));
//std::cout << " Inside Compute FG_KE....first_f " << std::endl;
     real nu_t_str_NE = 0.25 * (nu_t_str(i,j) + nu_t_str(i+1,j) + nu_t_str(i+1,j+1)+ nu_t_str(i,j+1)) ;
     real nu_t_str_SE = 0.25 * (nu_t_str(i,j) + nu_t_str(i+1,j) + nu_t_str(i+1,j-1) + nu_t_str(i,j-1));


     second_f(i,j) = idy * (nu_t_str_NE * (idy*(grid_.u(i,j,NORTH) - grid_.u(i,j,CENTER)) + idx * (grid_.v(i,j,EAST) - grid_.v(i,j,CENTER))) -
                    nu_t_str_SE * (idy*(grid_.u(i,j,CENTER) - grid_.u(i,j,SOUTH)) + idx * (grid_.v(i,j-1,EAST) - grid_.v(i,j,SOUTH))));

     third_f(i,j) = ((2/3) * idx) * (grid_.k(i,j,EAST) - grid_.k(i,j,CENTER));

    u2_dx(i,j) = (0.25 * idx) * (((grid_.u(i,j,CENTER) + grid_.u(i,j,EAST))	*(grid_.u(i,j,CENTER)+ grid_.u(i,j,EAST))-(grid_.u(i,j,WEST)+ grid_.u(i,j,CENTER))*(grid_.u(i,j,WEST)+ grid_.u(i,j,CENTER)))
            + gamma_*(std::fabs(grid_.u(i,j,CENTER)+ grid_.u(i,j,EAST))*(grid_.u(i,j,CENTER)- grid_.u(i,j,EAST))- std::fabs(grid_.u(i,j,WEST)+ grid_.u(i,j,CENTER))*(grid_.u(i,j,WEST)- grid_.u(i,j,CENTER))));



    uv_dy(i,j) = (0.25 * idy) * ( (grid_.v(i,j,CENTER) + grid_.v(i,j,EAST)) * (grid_.u(i,j,CENTER) + grid_.u(i,j,NORTH)) - (grid_.v(i,j,SOUTH) + grid_.v(i+1,j,SOUTH)) * (grid_.u(i,j,SOUTH) + grid_.u(i,j,CENTER))
              + gamma_ * ( std::fabs(grid_.v(i,j,CENTER) + grid_.v(i,j,EAST)) * (grid_.u(i,j,CENTER) - grid_.u(i,j,NORTH)) - std::fabs(grid_.v(i,j,SOUTH) + grid_.v(i+1,j,SOUTH)) * (grid_.u(i,j,SOUTH) - grid_.u(i,j,CENTER))));

    //u2_x2(i,j) = idx2 * (grid_.u(i,j,EAST) - 2 * grid_.u(i,j,CENTER) + grid_.u(i,j,WEST));


    //u2_y2(i,j) = idy2 * (grid_.u(i,j,NORTH) - 2 * grid_.u(i,j,CENTER) + grid_.u(i,j,SOUTH));
                        }

                } // End for loop u values



        for (int i = 1 ; i <= imax ; ++i)
                for( int j = 1 ; j < jmax ; ++j) {

                    if(grid_.isFluid(i,j) && grid_.isFluid(i,j+1)){

     real nu_t_str_NE = 0.25 * (nu_t_str(i,j) + nu_t_str(i+1,j) + nu_t_str(i+1,j+1)+ nu_t_str(i,j+1)) ;
     real nu_t_str_NW = 0.25 * (nu_t_str(i,j) + nu_t_str(i,j+1) + nu_t_str(i-1,j+1) + nu_t_str(i-1,j)) ;

     first_g(i,j) = idx * (nu_t_str_NE * (idx * (grid_.v(i,j,EAST) - grid_.v(i,j,CENTER)) + idy * (grid_.u(i,j,NORTH) - grid_.u(i,j,CENTER))) -
                             nu_t_str_NW * (idx * (grid_.v(i,j,CENTER) - grid_.v(i,j,WEST)) + idy * (grid_.u(i-1,j,NORTH) - grid_.u(i,j,WEST))));


     second_g(i,j) = idy2 * (nu_t_str(i,j+1) * (grid_.v(i,j,NORTH) - grid_.v(i,j,CENTER)) - nu_t_str(i,j)* (grid_.v(i,j,CENTER) - grid_.v(i,j,SOUTH)));


     third_g(i,j) = ((2/3) * idx) * (grid_.k(i,j,NORTH) - grid_.k(i,j,CENTER));

     uv_dx(i,j) = (0.25 * idx) * ((grid_.u(i,j,CENTER) + grid_.u(i,j,NORTH)) * (grid_.v(i,j,CENTER) + grid_.v(i,j,EAST)) - (grid_.u(i,j,WEST) + grid_.u(i,j+1,WEST)) * (grid_.v(i,j,WEST) + grid_.v(i,j,CENTER))
             + gamma_ * (std::fabs(grid_.u(i,j,CENTER) + grid_.u(i,j,NORTH)) * (grid_.v(i,j,CENTER) - grid_.v(i,j,EAST)) - std::fabs(grid_.u(i,j,WEST) + grid_.u(i,j+1,WEST)) * (grid_.v(i,j,WEST)-grid_.v(i,j,CENTER))));


    v2_dy(i,j) = (0.25 * idy) * ((grid_.v(i,j,CENTER) + grid_.v(i,j,NORTH)) * (grid_.v(i,j,CENTER) + grid_.v(i,j,NORTH)) - (grid_.v(i,j,CENTER) + grid_.v(i,j,SOUTH)) * (grid_.v(i,j,CENTER) + grid_.v(i,j,SOUTH))
            + gamma_ * (std::fabs(grid_.v(i,j,CENTER) + grid_.v(i,j,NORTH)) * (grid_.v(i,j,CENTER) - grid_.v(i,j,NORTH)) - std::fabs(grid_.v(i,j,CENTER) + grid_.v(i,j,SOUTH)) * (grid_.v(i,j,SOUTH) - grid_.v(i,j,CENTER))));

    //v2_x2(i,j) = idx2 * (grid_.v(i,j,EAST) - 2 * grid_.v(i,j,CENTER) + grid_.v(i,j,WEST)) ;

    //v2_y2(i,j) = idy2 * (grid_.v(i,j,NORTH) - 2 * grid_.v(i,j,CENTER) + grid_.v(i,j,SOUTH)) ;

                    }
                } //end for loop v values

        //PROGRESS("*********in f************");
        for (int i = 1 ; i < imax ; ++i)
        for( int j = 1 ; j <= jmax ; ++j)
            if(grid_.isFluid(i,j) && grid_.isFluid(i+1,j))
                grid_.f()(i,j) = grid_.u(i,j,CENTER) + dt_ * ((2 * first_f(i,j)) + second_f(i,j) - third_f(i,j) - u2_dx(i,j) - uv_dy(i,j) + gx_ );
    //grid_.f()(i,j) = grid_.u(i,j,CENTER) + dt_ * ( iRe * ( u2_x2(i,j) + u2_y2(i,j) ) - u2_dx(i,j) - uv_dy(i,j) + gx_ );

     // End for loop f values


        for (int i = 1 ; i <= imax ; ++i)
        for( int j = 1 ; j < jmax ; ++j)
            if(grid_.isFluid(i,j) && grid_.isFluid(i,j+1))
        grid_.g()(i,j) = grid_.v(i,j,CENTER) + dt_ * (first_g(i,j) + (2 * second_g(i,j)) - third_g(i,j) - uv_dx(i,j) - v2_dy(i,j)  + gy_) ;
     //grid_.g()(i,j) = grid_.v(i,j,CENTER) + dt_ * ( iRe * ( v2_x2(i,j) + v2_y2(i,j) ) - uv_dx(i,j) - v2_dy(i,j)  + gy_ );


        for ( int j = 1 ; j <= jmax ; ++j ){

        grid_.f()(0,j) = grid_.u()(0,j) ;

        grid_.f()(imax,j) = grid_.u()(imax,j);}

        for ( int i = 1 ; i <= imax ; ++i){
        grid_.g()(i,0) = grid_.v()(i,0);
        grid_.g()(i,jmax) = grid_.v()(i,jmax);}
         //PROGRESS("***********End f and g **********");

}


void FluidSimulator::composeRHS(){

	 real idx = 1/real(grid_.dx()); 
	 real idy = 1/real(grid_.dy());
	 real idt = 1.0/dt_ ;
	 int imax = grid_.p().getSize(0)-2 ;
	 int jmax = grid_.p().getSize(1)-2 ;

	for (int i = 1 ; i <= imax ; ++i)
        for( int j = 1 ; j <= jmax ; ++j)
            if(grid_.isFluid(i,j))
            grid_.rhs()(i,j) = idt * ( idx * (grid_.f(i,j,CENTER) - grid_.f(i,j,WEST)) + idy * (grid_.g(i,j,CENTER) - grid_.g(i,j,SOUTH))) ;

        // End for loop

} // End compose RHS

real FluidSimulator::f_nu(int i,int j){
    //int jmax = grid_.p().getSize(1)-2 ;

    real delta = grid_.dy();
    //real delta = 0.65;
       //std::cout << "delta =" << delta << std::endl;

    real R_delta,R_t,f_nu;

    R_delta = (std::sqrt(grid_.k(i,j,CENTER)) * delta )/nu ;

    if(grid_.e(i,j,CENTER)!=0){
        R_t = (grid_.k(i,j,CENTER) * grid_.k(i,j,CENTER))/(nu * grid_.e(i,j,CENTER));
    //std::cout << "grid_.e(i,j,CENTER) = "<<grid_.k(i,j,CENTER) << std::endl;
    }
    else
        R_t = 0;

    if(R_t!=0)
        f_nu = (1 - exp(-0.0165 * R_delta)) * (1 - exp(-0.0165 * R_delta)) * (1 + (20.5/R_t));
    else
        f_nu = 0;

    return f_nu;
}

 real FluidSimulator::f_1(int i,int j){

int jmax = grid_.p().getSize(1)-2 ;
      real delta = grid_.dy() * jmax;

     //real delta = 0.65;

     real R_delta,R_t,f_nu,f_1;

     R_delta = (std::sqrt(grid_.k(i,j,CENTER)) * delta )/nu ;


     //if(grid_.e(i,j,CENTER)!=0)
         R_t = (grid_.k(i,j,CENTER) * grid_.k(i,j,CENTER))/(nu * grid_.e(i,j,CENTER));
     //else
       //  R_t = 0;

     //if(R_t!=0)
         f_nu = (1 - exp(-0.0165 * R_delta)) * (1 - exp(-0.0165 * R_delta)) * (1 + (20.5/R_t));
     //else f_nu = 0;

     //if(f_nu!=0)
         f_1 = 1 + pow ((0.05/f_nu),3);
     //else
       //  f_1 = 0;

     return f_1;
 }

 real FluidSimulator::f_2(int i,int j){

     real R_t,f_2;

     //if(grid_.e(i,j,CENTER)!=0)
         R_t = (grid_.k(i,j,CENTER) * grid_.k(i,j,CENTER))/(nu * grid_.e(i,j,CENTER));
     //else
      //   R_t = 0;

     f_2 = 1 - exp(-(R_t*R_t));

     return f_2;
 }

 void FluidSimulator::compute_nu_t(){

     //Array<real> nu_t = grid_.nut();
     int imax = grid_.p().getSize(0)-2 ;
     int jmax = grid_.p().getSize(1)-2 ;

     for (int i = 1 ; i <= imax ; ++i)
     for( int j = 1 ; j <= jmax ; ++j){
         if(grid_.e(i,j,CENTER)!=0){
            grid_.nut()(i,j) = (c_nu * f_nu(i,j) * grid_.k(i,j,CENTER) * grid_.k(i,j,CENTER)) /grid_.e(i,j,CENTER);
            //std::cout << grid_.nut()(i,j) << std::endl;
         }
     else
           grid_.nut()(i,j) = 0;}
 }


 real FluidSimulator::nu_t(int i,int j){
//std::cout << " Inside Compute FG_KE....first_f " << std::endl;
    //Array<real> nu_t = grid_.nut();

    return grid_.nut()(i,j);

 }

 real FluidSimulator::nu_t_str(int i,int j){

     //Array<real> nu_t = grid_.nut();
     return nu + grid_.nut()(i,j);

 }

void FluidSimulator::computeKE(){

    real idx = 1/real(grid_.dx());
    real idx2 = idx * idx ;
    real idy = 1/real(grid_.dy());
    real idy2 = idy * idy ;
    //real iRe = 1.0 / Re_ ;

    int imax = grid_.p().getSize(0)-2 ;
    int jmax = grid_.p().getSize(1)-2 ;

    Array<real> k_new(imax+2,jmax+2);
    Array<real> e_new(imax+2,jmax+2);
    Array<real> nu_t(imax+2,jmax+2);


    for (int i = 1 ; i <=imax ; ++i){ // South

        k_new(i,0) = grid_.k()(i,0);
        e_new(i,0) = grid_.e()(i,0);

    }

    for (int i = 1 ; i <=imax ; ++i){ // North
        k_new(i,jmax + 1) = grid_.k()(i,jmax + 1);
        e_new(i,jmax + 1) = grid_.e()(i,jmax + 1);

    }

    for (int j = 1 ; j <=jmax ; ++j){ //West

        k_new(0,j)= grid_.k()(0,j);
        e_new(0,j)= grid_.e()(0,j);

    }

    for (int j = 1 ; j <=jmax ; ++j){ // East

        k_new(imax,j) = grid_.k()(imax,j);
        e_new(imax,j) = grid_.e()(imax,j);
    }



    for (int i = 1 ; i <= imax ; ++i)
            for( int j = 1 ; j <= jmax ; ++j) {

        //int i = 4; int j = 4;
                if(grid_.isFluid(i,j)){

                    real nu_t_east = (grid_.nut()(i,j) + grid_.nut()(i+1,j)) * 0.5; //std::cout << "nu_t_east =" << nu_t_east << std::endl;
                    real nu_t_west = (grid_.nut()(i,j) + grid_.nut()(i-1,j)) * 0.5; //std::cout << "nu_t_west =" << nu_t_west << std::endl;
                    real nu_t_north = (grid_.nut()(i,j) + grid_.nut()(i,j+1)) * 0.5;//std::cout << "nu_t_north =" << nu_t_north << std::endl;
                    real nu_t_south = (grid_.nut()(i,j) + grid_.nut()(i,j-1)) * 0.5;//std::cout << "nu_t_south =" << nu_t_south << std::endl;

                    real first_k = idx2 * (nu_t_east * (grid_.k(i,j,EAST) - grid_.k(i,j,CENTER)) -
                                    nu_t_west * (grid_.k(i,j,CENTER) - grid_.k(i,j,WEST)));

                    //std::cout << "first_k =" << first_k << std::endl;

                    real second_k = idy2 * (nu_t_north * (grid_.k(i,j,NORTH) - grid_.k(i,j,CENTER)) -
                                  nu_t_south * (grid_.k(i,j,CENTER) - grid_.k(i,j,SOUTH)));

                    //std::cout << "second_k =" << second_k << std::endl;

                    real third_k = (idx * 0.5)* ((grid_.u(i,j,CENTER) * (grid_.k(i,j,CENTER)+grid_.k(i,j,EAST))) -
                                    (grid_.u(i,j,WEST) * (grid_.k(i,j,WEST)+grid_.k(i,j,CENTER))) + gamma_ *
                                    ((std::fabs(grid_.u(i,j,CENTER)) * (grid_.k(i,j,CENTER) - grid_.k(i,j,EAST)))-
                                    (std::fabs(grid_.u(i,j,WEST)) * (grid_.k(i,j,WEST) - grid_.k(i,j,CENTER)))));

                    /*real third_k = (idx * 0.5) * (grid_.u(i,j,CENTER) * (grid_.k(i,j,CENTER)+grid_.k(i,j,EAST)) -
                                grid_.u(i,j,WEST) * (grid_.k(i,j,WEST)+grid_.k(i,j,CENTER)) + gamma_ *
                                (std::fabs(grid_.u(i,j,CENTER)) * (grid_.k(i,j,CENTER) - grid_.k(i,j,EAST)) -
                                std::fabs(grid_.u(i,j,WEST)) * (grid_.k(i,j,WEST) - grid_.k(i,j,CENTER))));*/

                    real four_k = (idy * 0.5) * (grid_.v(i,j,CENTER) * (grid_.k(i,j,CENTER)+grid_.k(i,j,NORTH)) -
                                grid_.v(i,j,SOUTH) * (grid_.k(i,j,SOUTH)+grid_.k(i,j,CENTER)) + gamma_ *
                                (std::fabs(grid_.v(i,j,CENTER)) * (grid_.k(i,j,CENTER) - grid_.k(i,j,NORTH)) -
                                std::fabs(grid_.v(i,j,SOUTH)) * (grid_.k(i,j,SOUTH) - grid_.k(i,j,CENTER))));

                    real du_dx = idx * (grid_.u(i,j,CENTER) - grid_.u(i,j,WEST));

                    real du_dy = (0.25 * idy) * (grid_.u(i,j,NORTH) + grid_.u(i-1,j,NORTH) - grid_.u(i,j,SOUTH) - grid_.u(i-1,j,SOUTH));

                    real dv_dy = idy * (grid_.v(i,j,CENTER) - grid_.u(i,j,SOUTH));
                    // Problem in next line
                    //std::cout << "i = " << i<< std::endl;
                    //std::cout << "j = " << j<< std::endl;
                    real dv_dx = (0.25 * idx) * (grid_.v(i,j,EAST) + grid_.v(i+1,j,SOUTH) - grid_.v(i,j,WEST) - grid_.v(i,j-1,WEST));


                    real gradugradut = (4 * du_dx * du_dx) + 2 * std::pow((du_dy + dv_dx),2) + (4 * dv_dy * dv_dy);

                    real five_k = 0.5 * nu_t(i,j) * pow(std::fabs(gradugradut),2);

                   //grid_.k()(i,j) = grid_.k(i,j,CENTER) + dt_* (first_k + second_k - third_k - four_k + five_k - grid_.e(i,j,CENTER));

                    k_new(i,j) = grid_.k(i,j,CENTER) + dt_* (first_k + second_k - third_k - four_k + five_k - grid_.e(i,j,CENTER));
                    if(k_new(i,j)<0){k_new(i,j) = 0.00001;}
                    //if(grid_.k()(i,j)<0){grid_.k()(i,j) = 0.00001;}
                    //std::cout << k_new(i,j) << std::endl;
                    //std::cout << "f_nu =" << f_nu(i,j) << std::endl;
                    real nu_t_f_nu_east = (f_nu(i,j)* grid_.nut()(i,j) + f_nu(i+1,j) * grid_.nut()(i+1,j)) * 0.5; //std::cout << "nu_t_fnu_east =" << nu_t_f_nu_east << std::endl;
                    real nu_t_f_nu_west = (f_nu(i,j)* grid_.nut()(i,j) + f_nu(i-1,j) * grid_.nut()(i-1,j)) * 0.5; //std::cout << "nu_t_fnu_east =" << nu_t_f_nu_west << std::endl;
                    real nu_t_f_nu_north = (f_nu(i,j)* grid_.nut()(i,j) + f_nu(i,j+1) *  grid_.nut()(i,j+1)) * 0.5; //std::cout << "nu_t_fnu_east =" << nu_t_f_nu_north << std::endl;
                    real nu_t_f_nu_south = (f_nu(i,j)* grid_.nut()(i,j) + f_nu(i,j-1) * grid_.nut()(i,j-1)) * 0.5;//std::cout << "nu_t_fnu_east =" << nu_t_f_nu_south << std::endl;

                    real first_e = (c_eps/c_nu) * (idx2 * (nu_t_f_nu_east * (grid_.e(i,j,EAST) - grid_.e(i,j,CENTER)) -
                                    nu_t_f_nu_west * (grid_.e(i,j,CENTER) - grid_.e(i,j,WEST))));

                    real second_e = (c_eps/c_nu) * (idy2 * (nu_t_f_nu_north * (grid_.e(i,j,NORTH) - grid_.e(i,j,CENTER)) -
                                  nu_t_f_nu_south * (grid_.e(i,j,CENTER) - grid_.e(i,j,SOUTH))));

                    real third_e = (idx * 0.5) * (grid_.u(i,j,CENTER) * (grid_.e(i,j,CENTER)+grid_.e(i,j,EAST)) -
                                (std::fabs(grid_.u(i,j,CENTER)) * (grid_.e(i,j,CENTER) - grid_.e(i,j,EAST)) -
                                std::fabs(grid_.u(i,j,WEST)) * (grid_.e(i,j,WEST) - grid_.e(i,j,CENTER))));

                    real four_e = (idy * 0.5) * (grid_.v(i,j,CENTER) * (grid_.e(i,j,CENTER)+grid_.e(i,j,NORTH)) -
                                                 grid_.u(i,j,WEST) * (grid_.e(i,j,WEST)+grid_.e(i,j,CENTER)) + gamma_ *
                                grid_.v(i,j,SOUTH) * (grid_.e(i,j,SOUTH)+grid_.e(i,j,CENTER)) + gamma_ *
                                (std::fabs(grid_.v(i,j,CENTER)) * (grid_.e(i,j,CENTER) - grid_.e(i,j,NORTH)) -
                                std::fabs(grid_.v(i,j,SOUTH)) * (grid_.e(i,j,SOUTH) - grid_.e(i,j,CENTER))));

                    real five_e = (c_1 * f_1(i,j) * grid_.k(i,j,CENTER)* 0.5) * gradugradut;

                    real six_e = (c_2 * f_2(i,j) * grid_.e(i,j,CENTER) * grid_.e(i,j,CENTER)) / grid_.k(i,j,CENTER);

                    //grid_.e()(i,j) = grid_.e(i,j,CENTER) + dt_ * (first_e + second_e - third_e - four_e + five_e - six_e);

                    e_new(i,j) = grid_.e(i,j,CENTER) + dt_ * (first_e + second_e - third_e - four_e + five_e - six_e);
                    //std::cout << e_new(i,j) << std::endl;
                }

                else
                {e_new(i,j) = 0.0; k_new(i,j) = 0.0;}

            }
grid_.k() = k_new;
grid_.e() = e_new;
//std::cout << " Inside Compute KE " << std::endl;
}
void FluidSimulator::updateVelocities(){

	 real idx = 1/real(grid_.dx()); 
	 real idy = 1/real(grid_.dy());
	//real idt = 1/dt_ ;
	int imax = grid_.p().getSize(0)-2 ;
	int jmax = grid_.p().getSize(1)-2 ;

	for (int i = 1 ; i < imax ; ++i)
        for( int j = 1 ; j <= jmax ; ++j)
            if(grid_.isFluid(i,j) && grid_.isFluid(i+1,j))
		grid_.u()(i,j) = grid_.f()(i,j) - dt_ * idx * ( grid_.p()(i+1,j) - grid_.p()(i,j) ) ; 
     // End for loop for u velocity


	for (int i = 1 ; i <=imax ; ++i)
        for( int j = 1 ; j < jmax ; ++j)
            if(grid_.isFluid(i,j) && grid_.isFluid(i,j+1))
		grid_.v()(i,j) = grid_.g()(i,j) - dt_ * idy * ( grid_.p()(i,j+1) - grid_.p()(i,j) ) ; 
     // End for loop for u velocity
} // End update Velocities 


void FluidSimulator::determineNextDT(){

	 real idx = 1.0/real(grid_.dx()); 
	 real dx_ = grid_.dx(); 
	 real idx2 = idx * idx ;
	 real idy = 1.0/real(grid_.dy()); 
	 real dy_ = grid_.dy();
	 real idy2 = idy * idy ; 	
	
	if (tau_ > 0.0){
		real abs_umax = grid_.u().get_abs_maxelement();
		//std::cout << "u max : " << abs_umax << std::endl;
		real abs_vmax = grid_.v().get_abs_maxelement();
		//std::cout << "v max : " << abs_vmax << std::endl;
		real second =  dx_ / abs_umax ;
		real third = dy_ / abs_vmax;
		real first = 0.5 * Re_/ (idx2 + idy2 );
		dt_ = std::min (second,third);
		dt_ = tau_ * std::min(first,dt_);		 
		} // end elseif	
	else{PROGRESS(" Time step length kept constant");}
	std::cout << "New time step taken :" << dt_ << std::endl;
} // end determineNextDT 

void FluidSimulator::initializeU(){
  
	int imax = grid_.p().getSize(0)-2 ;
	int jmax = grid_.p().getSize(1)-2 ;
	
	for (int i = 1 ; i < imax ; ++i)
        for( int j = 1 ; j <= (jmax * 0.5) ; ++j){
	  
	    grid_.u()(i,j) = 0.0 ;	  
	}
  
}

void FluidSimulator::boundaryCorrection(){
  int jmax = grid_.p().getSize(1)-2 ;
  for (int j = 1 ; j <(0.5* jmax) ; ++j){
    grid_.u()(0,j) = 0.0 ;  } 
}

/*void FluidSimulator::boundaryForKE(){
    int jmax = grid_.p().getSize(1)-2 ;
    if ( boundary_condition_W_ == "inflow") { // Inlet Turbulent properties
        real b;
        b = grid_.dy() * jmax ;
        for (int j = 1 ; j <=jmax ; ++j){

                    grid_.k()(0,j) = boundary_velocity_W_ * boundary_velocity_W_ * 0.003;
                    grid_.e()(0,j) =  (c_nu * pow((grid_.k()(0,j)),1.5) )/(0.03 * b);
            }
    }

}*/

void FluidSimulator::refreshBoundaries(){

	int imax = grid_.p().getSize(0)-2 ;
	int jmax = grid_.p().getSize(1)-2 ;

	if ( boundary_condition_S_ == "NOSLIP") {
		for (int i = 1 ; i <=imax ; ++i){
			grid_.u()(i,0) = 2 * boundary_velocity_S_ - grid_.u()(i,1);
                        grid_.v()(i,1) = 0.0 ;
                        //grid_.k()(i,1) = 0.0 ;
                        grid_.e()(i,0) = grid_.e()(i,1);
		}
		//PROGRESS("NOSLIP set to South boundary");
	}

	if ( boundary_condition_N_ == "NOSLIP") {
		for (int i = 1 ; i <=imax ; ++i){
			grid_.u()(i,jmax + 1) = 2 * boundary_velocity_N_ - grid_.u()(i,jmax);
			grid_.v()(i,jmax) = 0.0 ;
                       // grid_.k()(i,jmax) = 0.0 ;
                        grid_.e()(i,jmax +1 ) = grid_.e()(i,jmax);
		}
			//PROGRESS("NOSLIP set to North boundary");
	}

	if ( boundary_condition_W_ == "NOSLIP") {
		for (int j = 1 ; j <=jmax ; ++j){
			grid_.v()(0,j) = 2 * boundary_velocity_W_ - grid_.v()(1,j);
			grid_.u()(0,j) = 0.0 ;
                        //grid_.k()(0,j) = 0.0 ;
                        grid_.e()(0,j) = grid_.e()(1,j);
		}
		//PROGRESS("NOSLIP set to West boundary");
	}

	if ( boundary_condition_E_ == "NOSLIP") {
		for (int j = 1 ; j <=jmax ; ++j){
			grid_.v()(imax + 1,j) = 2 * boundary_velocity_E_ - grid_.v()(imax,j);
			grid_.u()(imax,j) = 0.0 ;
                        //grid_.k()(imax,j) = 0.0 ;
                        grid_.e()(imax+1,j) = grid_.e()(imax,j);
		}
		//PROGRESS("NOSLIP set to East boundary");
	}

	if ( boundary_condition_S_ == "inflow") {
		for (int i = 1 ; i <=imax ; ++i){
			grid_.u()(i,0) = - grid_.u()(i,1);
			grid_.v()(i,0) = boundary_velocity_S_ ;
		}
		//PROGRESS("Inflow set to South boundary");
	}

	if ( boundary_condition_W_ == "inflow") {
            real b;

            b = grid_.dy()* jmax;
		for (int j = 1 ; j <=jmax ; ++j){
			grid_.u()(0,j) = boundary_velocity_W_ ;
                        grid_.v()(0,j) =  - grid_.v()(1,j);
                        grid_.k()(0,j) = std::pow(boundary_velocity_W_,2) * 0.003;
                        grid_.e()(0,j) = (c_nu / (0.03 * b)) * pow((grid_.k()(0,j)),1.5);
		}
		//PROGRESS("Inflow set to West boundary");
	}
	
	if ( boundary_condition_N_ == "inflow") {
		for (int i = 1 ; i <=imax ; ++i){
			grid_.u()(i,jmax + 1) = - grid_.u()(i,jmax);
			grid_.v()(i,jmax) = boundary_velocity_N_ ;
		}
		//PROGRESS("Inflow set to North boundary");
	}

	if ( boundary_condition_E_ == "inflow") {
		for (int j = 1 ; j <=jmax ; ++j){
			grid_.u()(imax,j) = boundary_velocity_E_ ;
			grid_.v()(imax + 1,j) =  - grid_.v()(imax,j);			
		}
		//PROGRESS("Inflow set to East boundary");	  
	}

	if ( boundary_condition_E_ == "outflow") {
            //real b;
            //b = grid_.dy() * jmax ;
		for (int j = 1 ; j <=jmax ; ++j){
			grid_.u()(imax,j) = grid_.u()(imax - 1,j);
                        grid_.v()(imax + 1,j) = grid_.v()(imax,j);
                        grid_.e()(imax +1,j) = grid_.e()(imax,j);
                        grid_.k()(imax +1,j) = grid_.k()(imax,j);
                        //grid_.k()(imax,j) = boundary_velocity_W_ * boundary_velocity_W_ * 0.003;
                        //grid_.e()(imax +1,j) =  (c_nu * pow((grid_.k()(0,j)),1.5) )/(0.03 * b);
		}
		//PROGRESS("Outflow set to East boundary");
	}
	
	/*if ( flowField_ == "periodic") {
		for (int j = 1 ; j <=jmax ; ++j){
			grid_.u()(0,j) = grid_.u()(imax-1,j);
			grid_.u()(imax,j) = grid_.u()(1,j);
			grid_.v()(0,j) = grid_.v()(imax-1,j);
			grid_.v()(1,j) = grid_.v()(imax,j);
			grid_.v()(imax+1,j) = grid_.v()(2,j);
			grid_.p()(1,j) = grid_.p()(imax,j);	  
		}
		PROGRESS("Flow Field is periodic");
	}*/
	
	if ( boundary_condition_E_ == "periodic") {
		for (int j = 1 ; j <=jmax ; ++j){
			grid_.u()(0,j) = grid_.u()(imax-1,j);
			grid_.u()(imax,j) = grid_.u()(1,j);
			grid_.v()(0,j) = grid_.v()(imax-1,j);
			grid_.v()(1,j) = grid_.v()(imax,j);
			grid_.v()(imax+1,j) = grid_.v()(2,j);
			grid_.p()(1,j) = grid_.p()(imax,j);	  
		}
		//PROGRESS("Periodic Boundary at East face");
	}
	
	if ( boundary_condition_W_ == "periodic") {
		for (int j = 1 ; j <=jmax ; ++j){
			grid_.u()(0,j) = grid_.u()(imax-1,j);
			grid_.u()(imax,j) = grid_.u()(1,j);
			grid_.v()(0,j) = grid_.v()(imax-1,j);
			grid_.v()(1,j) = grid_.v()(imax,j);
			grid_.v()(imax+1,j) = grid_.v()(2,j);
			grid_.p()(1,j) = grid_.p()(imax,j);	  
		}
		//PROGRESS("Periodic Boundary at West face");
	}
	
	if ( boundary_condition_S_ == "SLIP") {
		for (int i = 1 ; i <=imax ; ++i){
			grid_.u()(i,0) = grid_.u()(i,1);
			grid_.v()(i,0) = 0.0 ;
		}
		//PROGRESS("SLIP set to South boundary");
	}

	if ( boundary_condition_N_ == "SLIP") {
		for (int i = 1 ; i <=imax ; ++i){
			grid_.u()(i,jmax + 1) = grid_.u()(i,jmax);
			grid_.v()(i,jmax) = 0.0 ;
		}
			//PROGRESS("SLIP set to North boundary");
	}
	
      /*  if ( boundary_condition_W_ == "inflow") { // Inlet Turbulent properties
            real b;
            b = grid_.dy() * jmax ;
            for (int j = 1 ; j <=jmax ; ++j){

                        grid_.k()(0,j) = boundary_velocity_W_ * boundary_velocity_W_ * 0.003;
                        grid_.e()(0,j) = (c_nu / (0.03 * b)) * pow((grid_.k()(0,j)),1.5);
                        //grid_.e()(0,j) =  (c_nu * pow((grid_.k()(0,j)),1.5) )/(0.03 * b);
                }
        }*/
	
} // End boundary loop 

void FluidSimulator::normalize(){
  
  
	int imax = grid_.p().getSize(0)-2 ;
	int jmax = grid_.p().getSize(1)-2 ;
	real Pavg , sum = 0.0 ;
	
	for (int i = 1 ; i <=imax ; i++)
        for( int j = 1 ; j <=jmax ; j++)
            if(grid_.isFluid(i,j))
		  sum += grid_.p()(i,j);		  		  

		  //std::cout << "summation in pressure : " << sum << std::endl ;
          //Pavg = sum / grid_.p().getSize();
            Pavg = sum / grid_.getNumFluid();
  
	for (int i = 1 ; i <=imax ; i++)
		for( int j = 1 ; j <=jmax ; j++) 
		    if(grid_.isFluid(i,j))
		  grid_.p()(i,j) = grid_.p()(i,j) - Pavg ;  
  
}

void FluidSimulator::printCourantNumber(){
  
    real C = 0.0 ;
    real idx = 1.0/real(grid_.dx()); 
    real idy = 1.0/real(grid_.dy()); 
    real abs_umax = grid_.u().get_abs_maxelement();
    real abs_vmax = grid_.v().get_abs_maxelement();
    
    C = dt_ * (abs_umax * idx + abs_vmax * idy );
    
    //CHECK_MSG( C < 1.0 , "Courant Number should be less than 1 ")
    
    std::cout << "Courant Number :" << C << std::endl;  
  
}

real FluidSimulator::getU(real x , real y){
  
	int imax = grid_.p().getSize(0)-2 ;
	int jmax = grid_.p().getSize(1)-2 ;
	int i = 0 ; int j = 0 ;
	real dx = grid_.dx();
	real dy = grid_.dy();
	real U;
	
	//std::cout << "length calculated "<<imax*dx << "Actual Length : " << xlength_ << std::endl; 
	CHECK_MSG((imax*dx <= xlength_) && (jmax*dy <= ylength_),"Given Coordinates are outside the domain ");
  
	do {
	  i = i + 1 ;
	  if ( dx * i >= x) {break ;} 
	  } while(i <= imax);
	do {
	  j = j + 1 ;
	  if ( dy * j >= y) {break ;} 
	   } while(j <= jmax);
  
	U =  0.5 * (grid_.u()(i+1,j) + grid_.u()(i-1,j));
  
	std::cout << "Value of U at said co-ordinates : " << U << std::endl;

	return U ;
}

real FluidSimulator::getV(real x , real y){
  
	int imax = grid_.p().getSize(0)-2 ;
	int jmax = grid_.p().getSize(1)-2 ;
	int i = 0 ; int j = 0 ;
	real dx = grid_.dx();
	real dy = grid_.dy();
	real V;
  
	CHECK_MSG((imax*dx <= xlength_) && (jmax*dy <= ylength_),"Given Coordinates are outside the domain ");
	
	do {
	i = i + 1 ;
	if ( dx * i >= x) {break ;} 
	} while(i <= imax);
	do {
	j = j + 1 ;
	if ( dy * j >= y) {break ;} 
	} while(j <= jmax);
  
	V =  0.25 * (grid_.v()(i,j+1) + grid_.v()(i,j-1));
  
	std::cout << "Value of V at said co-ordinates : " << V << std::endl;

    return V ;
}

real FluidSimulator::getP(real x , real y){
  
	int imax = grid_.p().getSize(0)-2 ;
	int jmax = grid_.p().getSize(1)-2 ;
	int i = 0 ; int j = 0 ;
	real dx = grid_.dx();
	real dy = grid_.dy();
	real P;
	
	CHECK_MSG((imax*dx <= xlength_) && (jmax*dy <= ylength_),"Given Coordinates are outside the domain ");
  
	do {
	i = i + 1 ;
	if ( dx * i >= x) {break ;} 
	} while(i <= imax);
	do {
	j = j + 1 ;
	if ( dy * j >= y) {break ;} 
	} while(j <= jmax);
  
	P =  0.25 * (grid_.p()(i+1,j) + grid_.p()(i-1,j) + grid_.p()(i,j+1) + grid_.p()(i,j-1));
  
	std::cout << "Value of P at said co-ordinates : " << P << std::endl;

	return P ;
}


void FluidSimulator::simulate(real duration){
	
	VTKWriter vtkWriter ( grid_, name_ , true, true );

	// Algorithm starts

	real t = 0.0; // 1:
	size_t n = 0; // 1:

	
	refreshBoundaries(); // Boundary Conditions on u and v 
	// line 3:

	while(t <= duration){ 

	// Select delta t line 4:
	determineNextDT();
	
	//Set boundary values for u and v (line 5:)
	refreshBoundaries(); 	
	
	// Compute F and G (line 6:)
	computeFG();
	grid_.f().print();
	
	//Compute right hand side of pressure equation (line 7:)
	composeRHS();
	
	//Solve pressure poisson equation (line 8 - 14 )
	solver_.solve(grid_);     
	//PROGRESS("SOR over");
	//Compute u and v in next time step(line 15:)
	updateVelocities();
	
	// Normalizing pressure after some iterations mentioned in par file
	if(n%nf_ == 0)
	  normalize();
	
		
	//grid_.u().print();

	// Writing output after some interval mentioned in par file
	if (n%outptinter_ == 0)
	vtkWriter.write();

	// Increment time(line 16:)
	t = t + dt_ ;

	std::cout << " Simulation completed for : " << t  << "seconds" <<  std::endl;
	//Increment step (line 17:)
	n++ ; 	

	} // end while loop

} // End of duration loop

void FluidSimulator::simulateTimeStepCount( unsigned int nrOfTimeSteps ){

        VTKWriter vtkWriter ( grid_, name_ , true, true,true,true );

	// Algorithm starts

	//real t = 0.0; // 1:
	size_t n = 0; // 1:

    //grid_.p().fill(fread.getParameter<real>("P_INIT"));
    //grid_.u().fill(fread.getParameter<real>("U_INIT"));
    //grid_.v().fill(fread.getParameter<real>("V_INIT"));
    grid_.k().fill(0.001);
    grid_.e().fill(0.001);

        refreshBoundaries(); // Boundary Conditions on u and v
	// line 3:
	if (name_ == "backstep"){
	boundaryCorrection();
	//Initialize lower domain to zero
	initializeU();}

        //if(turbulenceMode == 1){
          //  boundaryForKE();
        //}

	while(n <= nrOfTimeSteps){ 
	
	// Select delta t line 4:
	determineNextDT();
	
        if(turbulenceMode == 1){
            //std::cout << "Compute KE" << std::endl;
            computeKE();
            compute_nu_t();
        }
	
	//Set boundary values for u and v (line 5:)
	//refreshBoundaries(); 	
	
	// Compute F and G (line 6:)
        if(turbulenceMode ==1){
            computeFG_KE();}
        else
        computeFG();
	
	
	//Compute right hand side of pressure equation (line 7:)
	composeRHS();
	
	//Solve pressure poisson equation (line 8 - 14 )
	solver_.solve(grid_); 
	
	//Compute u and v in next time step(line 15:)
	updateVelocities();
	
	refreshBoundaries(); 
	
	// Boundary Correction at the West face of domain
	if (name_ == "backstep"){
	boundaryCorrection();}
	  
	// Normalizing pressure after some iterations mentioned in par file
	if(n%nf_ == 0)
	  normalize();
	
	//grid_.u().print();
	printCourantNumber();

	// Writing output after some interval mentioned in par file
	if (n%outptinter_ == 0)
	vtkWriter.write();

	// Increment time(line 16:)
	//t = t + dt_ ;

	std::cout << "Time Step finished: " << n  <<  std::endl;
	
	std::cout << "       " << std::endl;
	//Increment step (line 17:)
	n++ ; 	

	} // end while loop

}
