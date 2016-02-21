#include "StaggeredGrid.hh"





//Creating Staggered Grid for pressure(p_) and rhs for "(imax or jmax) + 2" elements because it includes extra layers on either side of the grid.
// Both pressure and rhs equations are 2D .Hence p_(i,j)
StaggeredGrid::StaggeredGrid( const FileReader & configuration):
p_(configuration.getParameter<int>("imax")+2,configuration.getParameter<int>("jmax")+2),
rhs_(configuration.getParameter<int>("imax")+2,configuration.getParameter<int>("jmax")+2),
u_(configuration.getParameter<int>("imax")+1,configuration.getParameter<int>("jmax")+2),
f_(configuration.getParameter<int>("imax")+1,configuration.getParameter<int>("jmax")+2),
v_(configuration.getParameter<int>("imax")+2,configuration.getParameter<int>("jmax")+1),
g_(configuration.getParameter<int>("imax")+2,configuration.getParameter<int>("jmax")+1),
k_(configuration.getParameter<int>("imax")+2,configuration.getParameter<int>("jmax")+2),
e_(configuration.getParameter<int>("imax")+2,configuration.getParameter<int>("jmax")+2),
nut_(configuration.getParameter<int>("imax")+2,configuration.getParameter<int>("jmax")+2),
Fluid(configuration.getParameter<int>("imax")+2,configuration.getParameter<int>("jmax")+2),
name_(configuration.getStringParameterBoundary("name"))
{
	imax_ = configuration.getParameter<int>("imax");
        //CHECK_MSG(imax_ >= 20, "Please keep imax more than 20 for better results");
	jmax_ = configuration.getParameter<int>("jmax");
        //CHECK_MSG(jmax_ >= 20, "Please keep jmax more than 20 for better results");
        if(configuration.checkParameter("Turbulence")){
        turbMode = configuration.getParameter<int>("Turbulence");}
	dx_ = configuration.getParameter<real>("xlength") / real (configuration.getParameter<int>("imax")) ; // data type conversion
	dy_ = configuration.getParameter<real>("ylength") / real (configuration.getParameter<int>("jmax")) ; // data type conversion
	CHECK_MSG(!dx_ == 0.0, " Width cannot be zero . use <xlength> as parameter");
	CHECK_MSG(!dy_ == 0.0, " Width cannot be zero . use <ylength> as parameter");	
	
	p_.fill(configuration.getParameter<real>("P_INIT"));
	u_.fill(configuration.getParameter<real>("U_INIT"));
	v_.fill(configuration.getParameter<real>("V_INIT"));
	rhs_.fill(0.0);	
	f_.fill(0.0);	
	g_.fill(0.0);
        //if(turbMode == 1){
        k_.fill(configuration.getParameter<real>("k_INIT"));
        e_.fill(configuration.getParameter<real>("e_INIT"));//}
        /*else{
        k_.fill(0.0);
        e_.fill(0.0);}*/

	Fluid.fill(10);
	
	if (name_ == "backstep"){
	  createRectangle(configuration.getParameter<real>("RectangleX1"),configuration.getParameter<real>("RectangleY1"),
                  configuration.getParameter<real>("RectangleX2"),configuration.getParameter<real>("RectangleY2"));}
        if (name_ == "ball"){
	  createCircle(configuration.getParameter<real>("CircleX"),configuration.getParameter<real>("CircleY"),
                configuration.getParameter<real>("CircleR"));}

//createRectangle(configuration.getParameter<real>("RectangleX1"),configuration.getParameter<real>("RectangleY1"),
                  //configuration.getParameter<real>("RectangleX2"),configuration.getParameter<real>("RectangleY2"));
	//createCircle(configuration.getParameter<real>("CircleX"),configuration.getParameter<real>("CircleY"),
                // configuration.getParameter<real>("CircleR"));
		
PROGRESS("******************Staggered Grid Created **************");
}

StaggeredGrid::StaggeredGrid( int xSize, int ySize, real dx, real dy ) :
p_(xSize+2,ySize+2),
rhs_(xSize+2,ySize+2),
u_(xSize+1,ySize+2),
f_(xSize+1,ySize+2),
v_(xSize+2,ySize+1),
g_(xSize+2,ySize+1),
k_(xSize+2,ySize+2),
e_(xSize+2,ySize+2),
nut_(xSize+2,ySize+2),
Fluid(xSize+2,ySize+2),dx_{dx},dy_{dy}{
Fluid.fill(100);
   p_.fill(1.0);
   u_.fill(1.0);
   v_.fill(1.0);

}


void StaggeredGrid::createRectangle(real x1,real y1,real x2,real y2){

    CHECK_MSG(x1 < (imax_ * dx_) && x2 < (imax_ * dx_), "X coordinates are outside the domain");
    CHECK_MSG(y1 < (jmax_ * dy_) && y2 < (jmax_ * dy_), "Y coordinates are outside the domain");
    //std::cout << "getSize(0) : " << Fluid.getSize(0) << std::endl;
 // std::cout << "getSize(1) : " << Fluid.getSize(1) << std::endl;
    //std::cout << "x2 : " << x2 << std::endl; 
  //std::cout << "y2 : " << y2 << std::endl; 

    for(int i = 1 ; i < Fluid.getSize(0) ; i++)
        for(int j = 1 ; j < Fluid.getSize(1) ; j++)
	  //if(dx_*(i) >= x1 && dx_*(i-0.5) <= x2 && dy_*(j-0.5) >= y1 && dy_*(j-0.5) <= y2)
            if(dx_*(i-0.5) >= x1 && dx_*(i-0.5) <= x2 && dy_*(j-0.5) >= y1 && dy_*(j-0.5) <= y2)
                setCellToObstacle(i,j);
}

void StaggeredGrid::createCircle(real x, real y, real r)
{
   for (int i=1; i<Fluid.getSize(0)-1; i++)
       for (int j=1; j<Fluid.getSize(1)-1; j++)
           if( ((x-(i-0.5)*dx_)*(x-(i-0.5)*dx_)) + ((y-j*dy_)*(y-j*dy_)) <= r*r)
              setCellToObstacle(i,j);
}

void StaggeredGrid::printfluid(){Fluid.print();}

/*StaggeredGrid ::StaggeredGrid(GrayScaleImage image):// real dx, real dy):
p_(image.width()+2,image.height()+2), 
rhs_(image.width()+2,image.height()+2),
u_(image.width()+1,image.height()+2),
f_(image.width()+1,image.height()+2),
v_(image.width()+2,image.height()+1),
g_(image.width()+2,image.height()+1),
k_(image.width()+2,image.height()+2),
e_(image.width()+1,image.height()+2),
Fluid(image.width()+2,image.height()+2),// dx_(dx), dy_(dy)
{
 PROGRESS("Reading Image File");
for (int i=1; i<Fluid.getSize(0); i++)
    for (int j=1; j<Fluid.getSize(1); j++)
        Fluid(i,j)= image.getElement(i,j); 
    
for (int i=1; i<Fluid.getSize(0); i++)
    for (int j=1; j<Fluid.getSize(1); j++)
        if(Fluid(i,j) == 0)
	  setCellToObstacle(i,j);
    //p_.fill(1.0);
   //u_.fill(1.0);
   //v_.fill(1.0);
}*/
