#ifndef STAGGERED_GRID_HH
#define STAGGERED_GRID_HH


#include "Types.hh"
#include "Array.hh"
//#include "Array.cc"
#include "Debug.hh"
#include "FileReader.hh"
#include "GrayScaleImage.hh"

class Obstacle;


//*******************************************************************************************************************
/*! Class for storing all arrays required for the simulation
*
* For now it only contains an array for the pressure and another array for the
* right hand side of the pressure equation.
* In following assignments this will be extended and also contain
* arrays for x/y velocity and some arrays for intermediate values.
*
* Feel free to add member functions or variables to this class, but please don't
* delete or rename any existing functions, since future skeletons rely on these functions.
*
*/
//*******************************************************************************************************************
class StaggeredGrid
{
public:
   // Constructors to manually create staggered grid
   StaggeredGrid ( int xSize, int ySize, real dx, real dy ); // TODO implement!

   // Constructor to create a staggered grid from a parsed configuration file
   StaggeredGrid ( const FileReader & configuration );       // TODO implement!
   
   StaggeredGrid (GrayScaleImage image ); //, real dx, real dy);


   // Getters / Setters for member variables
    Array<real> & p()    { return p_;    }
    Array<real> & rhs()  { return rhs_;  }
    Array<real> & u()  { return u_;  }
    Array<real> & v()  { return v_;  } 
    Array<real> & f()  { return f_;  }
    Array<real> & g()  { return g_;  }
    Array<real> & k()  { return k_;  }
    Array<real> & e()  { return e_;  }
    Array<real> & nut()  { return nut_;  }

    //Array<unsigned char> & Fluid() {return Fluid;}

    const Array<real> & p()   const { return p_;   }
    const Array<real> & rhs() const { return rhs_; }
    const Array<real> & u() const { return u_; }
    const Array<real> & v() const { return v_; }
    const Array<real> & f() const { return f_; }
    const Array<real> & g() const { return g_; }
    const Array<real> & k() const { return k_; }
    const Array<real> & e() const { return e_; }
    const Array<real> & nut() const { return nut_; }

    //const Array<unsigned char> & Fluid() const {return Fluid;}

   real dx() const { return dx_; }
   real dy() const { return dy_; }

   //wrapped access
   inline real& u(const int x, const int y, Direction dir);
   inline real& v(const int x, const int y, Direction dir);
   inline real& p(const int x, const int y, Direction dir);
   inline real& f(const int x, const int y, Direction dir);
   inline real& g(const int x, const int y, Direction dir);
   inline real& k(const int x, const int y, Direction dir);
   inline real& e(const int x, const int y, Direction dir);
   inline real& nut(const int x, const int y, Direction dir);


   inline bool isFluid(const int x,const int y){return Fluid(x,y);}
   inline int getNumFluid();

   void printfluid();

   void setCellToObstacle(int x,int y){Fluid(x,y)=0;}
   
   
   void createRectangle(real x1,real y1,real x2,real y2);
   void createCircle(real x, real y,real r);


protected:
   Array<real> p_;   //< pressure field
   Array<real> rhs_; //< right hand side of the pressure equation
   Array<real> u_;   // Extended variables
   Array<real> f_;
   Array<real> v_;
   Array<real> g_;
   Array<real> k_;
   Array<real> e_;
   Array<real> nut_;
  Array<unsigned char> Fluid;
  
   int imax_;
   int jmax_;
   std::string name_ ;
   real dx_;   //< distance between two grid points in x direction
   real dy_;   //< distance between two grid points in y direction
};


//inline bool StaggeredGrid::isFluid(const int x,const int y){return Fluid(x,y);}

inline int StaggeredGrid::getNumFluid(){

    int NumberofFluidCells = 0;

    for(int i = 1 ; i < Fluid.getSize(0)-1;i++)
        for(int j = 1 ; j < Fluid.getSize(1)-1;j++)
            NumberofFluidCells += Fluid(i,j);
	//std::cout << "Number of Fluid cells : " << NumberofFluidCells << std::endl;
    return NumberofFluidCells;
}

inline real& StaggeredGrid::p(const int x, const int y, Direction dir){
  
  switch(dir){

    case NORTH: if (Fluid(x,y+1)) return p_(x,y + 1); else return p_(x , y);  
        //return Fluid(x,y+1) ? (p_(x,y + 1)): p_(x , y) ;

    case EAST: if (Fluid(x+1,y)) return p_(x+1,y ); else return p_(x , y);
        //return Fluid(x + 1,y) ? (p_(x + 1,y )): p_(x , y) ;

    case SOUTH: if (Fluid(x,y-1)) return p_(x,y - 1); else return p_(x , y);
        //return Fluid(x,y-1) ? (p_(x,y - 1)): p_(x , y) ;

    case WEST: if (Fluid(x-1,y)) return p_(x-1,y); else return p_(x , y);
        //return Fluid(x-1,y) ? (p_(x-1,y)): p_(x , y) ;

    default : return p_(x,y);
  }
}

inline real& StaggeredGrid::u(const int x, const int y, Direction dir){
  
	 
  switch(dir){
    case NORTH :
        //(!Fluid(x,y+1)) ? (u_(x,y+1) = - u_(x,y); return u_(x,y+1);) : return u_(x,y + 1);
      if (!Fluid(x,y+1))u_(x,y+1) = - u_(x,y); return u_(x,y+1);

    case EAST:
        //return (Fluid(x+1,y)) ? (u(x+1,y)) : (u_(x+1,y) = 0.0 );
      if(Fluid(x+1,y)) return u(x+1,y,CENTER); else u_(x+1,y) = 0.0 ;  return u_(x+1,y); 

    case SOUTH: 
        //return (!Fluid(x,y-1)) ? (u_(x,y-1) = - u_(x,y)) : u_(x,y - 1);
      if (!Fluid(x,y-1)) u_(x,y-1) = - u_(x,y); return u_(x,y-1);

    case WEST:
        //return (!Fluid(x-1,y)) ? (u_(x-1 ,y) = 0.0 ) : u_(x-1,y) ;
      if(!Fluid(x-1,y)) u_(x-1 ,y) = 0.0;  return u_(x-1 ,y); 

    default : if(Fluid(x+1,y)) return u_(x,y); else u_(x,y) = 0.0; return u_(x,y);
  }//return u_(x,y);

}

inline real& StaggeredGrid::v(const int x, const int y, Direction dir){
  
  switch(dir) {

    case NORTH:
        //return (!Fluid(x,y+1)) ? (v_(x ,y+1) = 0.0) : v(x ,y +1);
      if (Fluid(x,y+1))return v(x,y+1,CENTER); else v_(x,y+1) = 0.0 ; return v_(x,y+1);

    case EAST:
        //return (!Fluid(x+1,y)) ? (v_( x+1,y) = - v_(x,y)) : v_(x +1,y);
      if (!Fluid(x+1,y)) v_( x+1,y) = - v_(x,y); return v_( x+1,y); 

    case SOUTH:
        //return (!Fluid(x,y-1)) ? (v_(x,y-1) = 0.0) : v_(x,y-1) ;
      if(!Fluid(x,y-1))(v_(x,y-1) = 0.0); return v_(x,y-1);
      
    case WEST:
        //return (!Fluid(x-1,y)) ? (v_(x-1,y) = - v_(x,y)) : v_(x-1,y) ;
      if (Fluid(x-1,y))v_(x-1,y) = - v_(x,y);return v_(x-1,y); 

    default : if(!Fluid(x,y+1))v_(x,y)= 0.0;return v_(x,y); 
  }//    return v_(x,y);

}

inline real& StaggeredGrid::k(const int x, const int y, Direction dir){
    switch (dir) {
      case NORTH:
        if(Fluid(x,y+1)) return k_(x,y+1); else return k_(x,y+1)= 0.0;

      case SOUTH:
        if(Fluid(x,y-1)) return k_(x,y-1); else return k_(x,y-1)=0.0;

      case WEST:
        if (Fluid(x-1,y)) return k_(x-1,y); else return k_(x-1,y)=0.0;

      case EAST:
        if (Fluid(x+1,y)) return k_(x+1,y); else return k_(x+1,y)=0.0;

      default : return k_(x,y); }

}

inline real& StaggeredGrid::e(const int x, const int y, Direction dir){
    switch (dir) {
      case NORTH:
        if(Fluid(x,y+1)) return e_(x,y+1); else return e_(x,y);

      case SOUTH:
        if(Fluid(x,y-1)) return e_(x,y-1); else return e_(x,y);

      case WEST:
        if (Fluid(x-1,y)) return e_(x-1,y); else return e_(x,y);

      case EAST:
        if (Fluid(x+1,y)) return e_(x+1,y); else return e_(x,y);

      default : return e_(x,y); }

}

inline real& StaggeredGrid::f(const int x, const int y, Direction dir){
  switch (dir) {
    case CENTER: 
      if(Fluid(x+1,y)) return f_(x,y); else return u_(x,y);      
    

    case WEST:
        //return (!Fluid(x-1,y)) ? u_(x-1,y) : f_(x-1,y) ;
    if (Fluid(x-1,y)) return f_(x-1,y); else return u_(x-1,y); 

    default : return f_(x,y); }

}

inline real& StaggeredGrid::g(const int x, const int y, Direction dir){
  
  switch(dir) {
  
    case CENTER:
      if(Fluid(x,y+1)) return g_(x,y); else return v_(x,y);      
    

    case SOUTH:
        //return (!Fluid(x,y-1)) ? v_(x,y-1) : g_(x,y-1) ; 
      if (Fluid(x,y-1)) return g_(x,y-1);else return v_(x,y-1);

    default : return g_(x,y);}
}


#endif //STAGGERED_GRID_HH

