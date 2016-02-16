#include "StaggeredGrid.hh"





//Creating Staggered Grid for pressure(p_) and rhs for "(imax or jmax) + 2" elements because it includes extra layers on either side of the grid.
// Both pressure and rhs equations are 2D .Hence p_(i,j)
StaggeredGrid::StaggeredGrid( const FileReader & configuration):
p_(configuration.getParameter<int>("imax")+2,configuration.getParameter<int>("jmax")+2),
rhs_(configuration.getParameter<int>("imax")+2,configuration.getParameter<int>("jmax")+2)
{
	dx_ = configuration.getParameter<real>("xlength") / real (configuration.getParameter<int>("imax")) ; // data type conversion
	dy_ = configuration.getParameter<real>("ylength") / real (configuration.getParameter<int>("jmax")) ; // data type conversion

PROGRESS("******************Staggered Grid Created **************");
}

StaggeredGrid::StaggeredGrid( int xSize, int ySize, real dx, real dy ) :
p_(xSize+2,ySize+2),rhs_(xSize+2,ySize+2),dx_{dx},dy_{dy}{
}
