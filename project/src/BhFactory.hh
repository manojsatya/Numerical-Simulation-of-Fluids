#ifndef BH_Factory
#define BH_Factory

#include "StaggeredGrid.hh"
#include "FileReader.hh"
#include "Types.hh"
#include <string>
#include "Debug.hh"
#include <sstream>
#include <time.h>
#include <algorithm>
#include <vector>
//This class handles the boundaries

class BhFactory
{
public:
   // Constructor to create a BhFactory from a parsed configuration file
   BhFactory (FileReader & configuration );
   //checks whether boundary conditions make sense
   void Check_for_sense();
   // solve the pressure equation on the staggered grid;
   void refreshBoundaries( StaggeredGrid & grid );
   //for string comparison
   bool isequals(std::string& a, const std::string& b);
   int internal_time ; //for periodic boundary
   int turb_mode;
   real a,b; //dimension of the canal

private:

   real  boundary_velocity_N,boundary_velocity_E,boundary_velocity_W,boundary_velocity_S;
   std::string  boundary_condition_N, boundary_condition_E, boundary_condition_W,boundary_condition_S;

};

#endif //SOR_SOLVER_HH




