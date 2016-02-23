#include "SORSolver.hh"


 SORSolver::SORSolver (real xlength_, real ylength_, int imax_, int jmax_, real omega_, real eps_, real itermax_,std::string name_): xlength(xlength_),ylength(ylength_),
 imax(imax_),jmax(jmax_),omega(omega_),eps(eps_),itermax(itermax_), name(name_)
 {
	 //manual constructor for SORSolver
 }

 SORSolver::SORSolver ( FileReader & configuration )
{

	xlength     = configuration.getParameter("xlength");
	ylength     = configuration.getParameter("ylength");
	imax        = configuration.getParameter("imax");
	jmax        = configuration.getParameter("jmax");
	omega       = configuration.getParameter("omg");
	eps         = configuration.getParameter("eps");
	itermax     = configuration.getParameter("itermax");
	name        = configuration.getStringParameter("name");
	check_freq  = configuration.getParameter("checkfrequency");

 }

bool SORSolver::check_existence(Array<real> *rhs,int xSize, int ySize )
 {
	real check=0;
	 for(int j=1; j<ySize-1 ; ++j)
		 for(int i=1; i<xSize-1 ; ++i)
		 {
           check +=(*rhs)(j,i);
		 }
	 if(check <= 1e-5)
		 return true;
		 else
	 return false;
 }
 void SORSolver::boundary_handling(Array<real> *p_, int xSize, int ySize )
 {
	 for (int i=1; i<xSize-1; ++i)
	  {
	 	 (*p_)(0,i)       = (*p_)(1,i);
	 	 (*p_)(ySize-1,i) = (*p_)(ySize-2,i);

	  }

	  for (int j=0; j<ySize-1; ++j)
	  {
	 	 (*p_)(j,0)       = (*p_)(j,1);
	 	 (*p_)(j,xSize-1) = (*p_)(j,xSize-2);
	  }
 }

 bool SORSolver::solve( StaggeredGrid & grid )
 {
 clock_t start_time;
 start_time = clock();

// PRG_LEVEL("SORSolver started");
 ASSERT_MSG(grid.p().getSize(0)==imax+2,"Dimension mismatch between SORSolver p_ and Staggered Grid");
 ASSERT_MSG(grid.p().getSize(1)==jmax+2,"Dimension mismatch between SORSolver p_ and Staggered Grid");

 real dx = xlength/imax;
 real dy = ylength/jmax;

 //pointer to pressure and RHS for convenience
 Array<real> * p_   = &grid.p();
 Array<real> * rhs  = &grid.rhs();
 int xSize = p_->getSize(0);
 int ySize = p_->getSize(1);

 //initializing norm, residual and iterations
real res_norm = 100.;
real res_norm_check = res_norm;// for checking
real residual = 0;
size_t iteration = 0;

CHECK_MSG(check_existence(rhs,xSize, ySize),"SOR Solver: Solution does not exist for the given RHS");

// Initial Boundary Handling
 boundary_handling(p_, xSize, ySize);


// constants used in the loop
const real dx_sq_inv    = 1.0/(dx*dx);
const real dy_sq_inv    = 1.0/(dy*dy);
const real C_factor_inv = ( (dx*dx*dy*dy) / (2.0*(dx*dx+dy*dy)) );
const real C_factor= 1.0/C_factor_inv;
const real res_factor   = 1.0/grid.getNumFluid(); //divided by number of fluid cells
real neighbours;
std::stringstream msg1, msg2;

while(res_norm_check>eps && iteration<=itermax )
{
res_norm = 0;
 for(int j=1; j<ySize-1 ; ++j)
	 for(int i=1; i<xSize-1 ; ++i)
	 {
		 if(grid.isFluid(j,i))
		 {
		 neighbours =  ( grid.p(j,i,WEST) + grid.p(j,i,EAST) )*dx_sq_inv + ( grid.p(j,i,SOUTH) + grid.p(j,i,NORTH))*dy_sq_inv;
		 residual = C_factor*p_->getData(j,i) + rhs->getData(j,i) - neighbours ;
		 p_->getData(j,i) = (1-omega)*p_->getData(j,i) +  omega*C_factor_inv * ( neighbours - rhs->getData(j,i));
		 res_norm += residual*residual;
		 }

	 }



boundary_handling(p_,xSize,ySize);

// This sqrt is  only the expensive part in residual calculation; rest is more like free; here check_freq = 1 is also OK
if(iteration%check_freq)
res_norm_check = sqrt(res_norm*res_factor);

iteration++;
//msg1<<"SOR SOLVER: Iteration "<<iteration<<"/"<<itermax<<"\nCurrent residual:"<< res_norm<<std::endl;
//INFO(msg1.str());
}
msg2<<"SORSolver Finished Successfully.\n\t  Total Time Taken: "<<(clock()-start_time)*1.0e-6<<", seconds, Total iterations: "<<iteration<<"\n\t residual:"<< res_norm;
PRG_LEVEL(msg2.str());
if(iteration < itermax)
	return true;
else
	return false;

 }
