#ifndef SOR_SOLVER_HH
#define SOR_SOLVER_HH

#include "StaggeredGrid.hh"
#include "Array.hh"


class SORSolver
{
public:
   // Constructor to manually create SORSolver
   SORSolver (size_t itermax , real omg, real eps );

   // Constructor to create a SORSolver from a parsed configuration file
   SORSolver ( const FileReader & configuration );


   // solve the pressure equation on the staggered grid
   bool solve( StaggeredGrid & grid );

private:
   // TODO add solver parameters here as member
	size_t itermax_;
	real omg_;
	real eps_;
};






#endif //SOR_SOLVER_HH




