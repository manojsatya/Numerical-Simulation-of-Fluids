#ifndef __OBSTACLE_H__
#define __OBSTACLE_H__


#include "FileReader.hh"
#include "StaggeredGrid.hh"
#include <string>
#include <iostream>


class Obstacle
{
  public:
    Obstacle( const FileReader & conf );
    
    StaggeredGrid & grid();
    const StaggeredGrid & grid() const;
    
    void Rectangle();
    void Circle();
    
    
    
  private:
    
  StaggeredGrid grid_; 
  
  real x1_ ; // Lower x co ord of rectangle
  real y1_ ; // Lower y co ord of rectangle
  real x2_ ; // upper right x co ord of rectangle
  real y2_ ; // upper right y co ord of rectangle
  real x_ ; // Center x co ord of circle
  real y_ ; // Center y co ord of circle
  real r_; // Radius of circle

    
    
    
    
};


#endif