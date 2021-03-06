#ifndef STAGGERED_GRID_HH
#define STAGGERED_GRID_HH


#include "Types.hh"
#include "Array.hh"
#include "FileReader.hh"
#include "GrayScaleImage.hh"

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
   StaggeredGrid ( int xSize, int ySize, real dx, real dy );

   // Constructor to create a staggered grid from a parsed configuration file
   StaggeredGrid (FileReader & configuration );


   // Getters / Setters for member variables
   Array<real> & p()    { return p_;    }
   Array<real> & rhs()  { return rhs_;  }
   Array<real> & u()    { return u_;    }
   Array<real> & v()    { return v_;    }
   Array<real> & f()    { return f_;    }
   Array<real> & g()    { return g_;    }
   Array<real> & k_old()   { return k_old_;    }
   Array<real> & e_old()   { return e_old_;    }
   Array<real> & k_new()   { return k_new_;    }
   Array<real> & e_new()   { return e_new_;    }
   Array<real> & nu_t()    { return nu_t_;     }
   Array<real> & delta()   { return delta_;    }

   const Array<real> & p()   const { return p_;   }
   const Array<real> & rhs() const { return rhs_; }
   const Array<real> & u()   const { return u_;   }
   const Array<real> & v()   const{ return v_;    }
   const Array<real> & f()   const{ return f_;    }
   const Array<real> & g()   const{ return g_;    }
   const Array<real> & k_old()   const{ return k_old_;    }
   const Array<real> & e_old()   const{ return e_old_;    }
   const Array<real> & k_new()   const{ return k_new_;    }
   const Array<real> & e_new()   const{ return e_new_;    }
   const Array<real> & nu_t()    const{ return nu_t_;     }
   const Array<real> & delta()   const{ return delta_;    }

   real dx() const { return dx_; }
   real dy() const { return dy_; }
   int imax()const {return imax_;}
   int jmax()const {return jmax_;}
   inline bool isFluid(const int x, const int y);
   inline int getNumFluid();
   //wrapped accesses
   inline real u(const int x,const int y, DIRECTION dir);
   inline real v(const int x,const int y, DIRECTION dir);
   inline real p(const int x,const int y, DIRECTION dir);
   inline real k_old(const int y,const int x, DIRECTION dir);
   inline real e_old(const int y,const int x, DIRECTION dir);
   void setCellToObstacle(const int x, const int y);
   void createRectangle(real x1, real y1, real x2, real y2);
   void createCircle (real x, real y, real r);
   void load_image(const std::string & pngFilename);
   void save_image(const std::string & pngFilename);

protected:
   Array<real> p_;   //< pressure field
   Array<real> rhs_; //< right hand side of the pressure equation
   Array<real> u_;
   Array<real> v_;
   Array<real> f_;
   Array<real> g_;
   Array<real> k_old_;
   Array<real> e_old_;
   Array<real> k_new_;
   Array<real> e_new_;
   Array<real> nu_t_;
   Array<CELL_TYPE> flag_;
   Array<real> delta_;//to store closest wall distance

   real dx_;   //< distance between two grid points in x direction
   real dy_;   //< distance between two grid points in y direction
   int imax_;
   int jmax_;
};


bool StaggeredGrid::isFluid(const int y, const int x)
{
	if(flag_(y,x)==FLUID)
		return true;
	else
		return false;
}

int StaggeredGrid::getNumFluid()
{
	int ctr=0;
	for(int i=0; i<flag_.y_size; ++i)
	{
		for(int j=0; j<flag_.x_size;++j)
		{
			if(isFluid(i,j))
				++ctr;
		}
	}
	return ctr;
}

real StaggeredGrid::u(const int y,const int x, DIRECTION dir)
{
	switch (dir){
	case NORTH:
		if(isFluid(y+1,x) && isFluid(y+1,x+1))
			return u_(y+1,x);
		else if(!isFluid(y+1,x)&& !isFluid(y+1,x+1))
			return -u_(y,x);
		else
			return 0;
		break;

	case EAST:
		if(isFluid(y,x+1) && isFluid(y,x))
			return u_(y,x+1);
		else
			return 0;
		break;

	case WEST:
		if(isFluid(y,x-1) && isFluid(y,x))
			return u_(y,x-1);
		else
			return 0;
		break;

	case SOUTH:
		if(isFluid(y-1,x) && isFluid(y-1,x+1))
			return u_(y-1,x);
		else if(!isFluid(y-1,x) && !isFluid(y-1,x+1))
			return -u_(y,x);
		else
			return 0;
		break;

	case CENTER:
		return u_(y,x);
	default: ERROR("An unspecified direction in wrapped access");
	return 0;

	}

}
real StaggeredGrid::v(const int y,const int x, DIRECTION dir)
{
	switch (dir){
	case NORTH:
		if(isFluid(y,x) && isFluid(y+1,x))
			return v_(y+1,x);
		else
			return 0;
		break;

	case EAST:
		if(isFluid(y,x+1) && isFluid(y+1,x+1))
           return v_(y,x+1);
		else if(!isFluid(y,x+1) && !isFluid(y+1,x+1))
			return -v_(y,x);
		else
			return 0;
		break;

	case WEST:
		if(isFluid(y,x-1) && isFluid(y+1,x-1))
           return v_(y,x-1);
		else if(!isFluid(y,x-1) && !isFluid(y+1,x-1))
			return -v_(y,x);
		else
			return 0;
		break;

	case SOUTH:
		if(isFluid(y,x) && isFluid(y-1,x))
			return v_(y-1,x);
		else
			return 0;
		break;

	case CENTER:
		return v_(y,x);
	default:ERROR("An unspecified direction in wrapped access");
	return 0;
	}
}

real StaggeredGrid::p(const int y,const int x, DIRECTION dir)
{
	switch (dir){
	case NORTH:
		if(isFluid(y+1,x))
			return p_(y+1,x);
		else
			return p_(y,x);
		break;

	case EAST:
		if(isFluid(y,x+1))
			return p_(y,x+1);
		else
			return p_(y,x);
		break;

	case WEST:
		if(isFluid(y,x-1))
			return p_(y,x-1);
		else
			return p_(y,x);
		break;

	case SOUTH:
		if(isFluid(y-1,x))
			return p_(y-1,x);
		else
			return p_(y,x);
		break;

	case CENTER:
		return p_(y,x);

	default:ERROR("An unspecified direction in wrapped access");
	return 0;
	}
}
real StaggeredGrid::k_old(const int y,const int x, DIRECTION dir)
{
	switch (dir){
	case NORTH:
		if(isFluid(y+1,x))
			return k_old_(y+1,x);
		else
			return 0;//-k_old_(y,x);
		break;

	case EAST:
		if(isFluid(y,x+1))
			return k_old_(y,x+1);
		else
			return 0;//-k_old_(y,x);
		break;

	case WEST:
		if(isFluid(y,x-1))
			return k_old_(y,x-1);
		else
			return 0;//-k_old_(y,x);
		break;

	case SOUTH:
		if(isFluid(y-1,x))
			return k_old_(y-1,x);
		else
			return 0;//-k_old_(y,x);
		break;

	case CENTER:
		return k_old_(y,x);

	default:ERROR("An unspecified direction in wrapped access");
	return 0;
	}
}

real StaggeredGrid::e_old(const int y,const int x, DIRECTION dir)
{
	switch (dir){
	case NORTH:
		if(isFluid(y+1,x))
			return e_old_(y+1,x);
		else
			return e_old_(y,x);
		break;

	case EAST:
		if(isFluid(y,x+1))
			return e_old_(y,x+1);
		else
			return e_old_(y,x);
		break;

	case WEST:
		if(isFluid(y,x-1))
			return e_old_(y,x-1);
		else
			return e_old_(y,x);
		break;

	case SOUTH:
		if(isFluid(y-1,x))
			return e_old_(y-1,x);
		else
			return e_old_(y,x);
		break;

	case CENTER:
		return e_old_(y,x);

	default:ERROR("An unspecified direction in wrapped access");
	return 0;
	}
}




#endif //STAGGERED_GRID_HH

