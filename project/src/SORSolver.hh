#ifndef SOR_SOLVER_HH
#define SOR_SOLVER_HH

#include "StaggeredGrid.hh"
#include "FileReader.hh"
#include <string>
#include "Debug.hh"
#include <cmath>
#include <sstream>
#include <time.h>

class SORSolver
{
public:
   // Constructor to manually create SORSolver
   SORSolver (real xlength_, real ylength_, int imax_, int jmax_, real omega_, real eps_, real itermax_, std::string name_);

   // Constructor to create a SORSolver from a parsed configuration file
   SORSolver (FileReader & configuration );

   // solve the pressure equation on the staggered grid
   bool solve( StaggeredGrid & grid );

   bool check_existence(Array<real> *rhs, int xSize, int ySize);

   real max_iter()
   {
	   return itermax;
   }

private:
   real xlength;
   real ylength;
   int  imax;
   int  jmax;
   real omega;
   real eps;
   real itermax;
   std::string name;
   int check_freq;
   void boundary_handling(Array<real> *data, int xSize, int ySize );
};



#endif //SOR_SOLVER_HH




