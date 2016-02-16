#include "Obstacle.hh"
#include "Debug.hh"


Obstacle::Obstacle( const FileReader & conf ):grid_(conf),
	x1_{conf.getParameter<real>("RectangleX1")},
	y1_{conf.getParameter<real>("RectangleY1")},
	x2_{conf.getParameter<real>("RectangleX2")},
	y2_{conf.getParameter<real>("RectangleY2")},
	x_{conf.getParameter<real>("CircleX")},
	y_{conf.getParameter<real>("CircleY")},
	r_{conf.getParameter<real>("CircleR")} {}
	
StaggeredGrid &Obstacle::grid() {return grid_;}

const StaggeredGrid &Obstacle::grid() const{return grid_;}
	
void Obstacle::Rectangle(){
  
  //std::cout << "x2_ : " << x2_ << std::endl; 
  //std::cout << "y2_ : " << y2_ << std::endl; 
   grid_.createRectangle(x1_ ,y1_ ,x2_ ,y2_ );
   
} 
 
void Obstacle::Circle(){
      
      grid_.createCircle(x_, y_, r_);     
      
    }	
	
	