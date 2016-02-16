#include "StaggeredGrid.hh"





//Creating Staggered Grid for pressure(p_) and rhs for "(imax or jmax) + 2" elements because it includes extra layers on either side of the grid.
// Both pressure and rhs equations are 2D .Hence p_(i,j)
StaggeredGrid::StaggeredGrid( const FileReader & configuration):
p_(configuration.getParameter<int>("imax")+2,configuration.getParameter<int>("jmax")+2),
rhs_(configuration.getParameter<int>("imax")+2,configuration.getParameter<int>("jmax")+2),
u_(configuration.getParameter<int>("imax")+1,configuration.getParameter<int>("jmax")+2),
f_(configuration.getParameter<int>("imax")+1,configuration.getParameter<int>("jmax")+2),
v_(configuration.getParameter<int>("imax")+2,configuration.getParameter<int>("jmax")+1),
g_(configuration.getParameter<int>("imax")+2,configuration.getParameter<int>("jmax")+1)
{
	imax_ = configuration.getParameter<int>("imax");
	CHECK_MSG(imax_ > 20, "Please keep imax more than 20 for better results");
	jmax_ = configuration.getParameter<int>("jmax");
	CHECK_MSG(jmax_ > 20, "Please keep jmax more than 20 for better results");
	dx_ = configuration.getParameter<real>("xlength") / real (configuration.getParameter<int>("imax")) ; // data type conversion
	dy_ = configuration.getParameter<real>("ylength") / real (configuration.getParameter<int>("jmax")) ; // data type conversion
	CHECK_MSG(!dx_ == 0.0, " Width cannot be zero . use <xlength> as parameter");
	CHECK_MSG(!dy_ == 0.0, " Width cannot be zero . use <ylength> as parameter");	
	
		
	p_.fill(configuration.getParameter<real>("P_INIT"));
	rhs_.fill(0.0);
	u_.fill(configuration.getParameter<real>("U_INIT"));
	f_.fill(0.0);
	v_.fill(configuration.getParameter<real>("V_INIT"));
	g_.fill(0.0);
	

PROGRESS("******************Staggered Grid Created **************");
}

StaggeredGrid::StaggeredGrid( int xSize, int ySize, real dx, real dy ) :
p_(xSize,ySize),rhs_(xSize,ySize),u_(xSize-1,ySize),f_(xSize,ySize-1),v_(xSize-1,ySize),g_(xSize,ySize-1),dx_{dx},dy_{dy}{
}
