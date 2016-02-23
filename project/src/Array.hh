#ifndef ARRAY_HH
#define ARRAY_HH


#include "Types.hh"
#include <memory>
#include "Debug.hh"
#include <sstream>
#include <vector>
#include <iostream>
#include <string>
#include <functional>
#include <math.h>
#include <limits>
#include <vector>
//*******************************************************************************************************************
/*!  Array class for 1,2 and 3 dimensions
*
*    - all elements should be stored in a contiguous chunk of memory ( no vector<vector> ! )
*/
//*******************************************************************************************************************
//New cange in Array class; switches bounds check 'on' only if debug mode is switched on; (performance reasons)

template<typename T> class Array
{

public:
   T *Data;
   int z_size;
   int y_size;
   int x_size;


   //===================================================================================================================
   //
   //  Constructors
   //
   //===================================================================================================================

   //default constructor
   Array():Data(nullptr),z_size(0),y_size(0),x_size(0)
   {}

   //1D constructor
   Array(int xSize):Data (nullptr),z_size(1),y_size(1),x_size(xSize)
   {
   	if (xSize>=1)
   	{
   	 Data = new T[x_size];
   	 fill(T());
     }
   else
    {
    std::stringstream msg;
    msg<<"Dimensions less than 1 not possible: Here "<<xSize<<" ,use default constructor if reqd.";
    ERROR(msg.str());
    }
   }
  //2D constructor
   Array( int ySize , int xSize ):Data (nullptr),z_size(1),y_size(ySize),x_size(xSize)
   {
   	 if (xSize>=1 && ySize>=1)
   	 {
   		 Data = new T[x_size*ySize];
   		 fill(T());

   	 }
   	 else
   	 {
   		 std::stringstream msg;
   		 msg<<"dimensions less than 1 not possible: Here ("<<xSize<<","<<ySize<<")";
   		 ERROR(msg.str());
   	 }

   }
//3D constructor
   Array( int zSize, int ySize, int xSize ):Data (nullptr),z_size(zSize),y_size(ySize),x_size(xSize)
   {
   	 if (xSize>=1 && (ySize>=1 && zSize>=1))
   	 {
   		 Data = new T[x_size*ySize*zSize];
   		 fill(T());
   	 }
   	 else
   	 {
   		 std::stringstream msg;
   		 msg<<"dimensions less than 1 not possible: Here ("<<xSize<<","<<ySize<<","<<zSize<<")";
   		 ERROR(msg.str());
   	 }

   }


  Array(const Array& s):Data(nullptr),z_size(s.z_size),y_size(s.y_size),x_size(s.x_size)
   {
   	*this = s;
   }

   Array& operator= (const Array& s)
   {
   	if(this != &s)
   	{
   	   cleanup();
   	   z_size = s.z_size;
           y_size = s.y_size;
           x_size = s.x_size;
           // performing a deep-copy
           Data = new T[x_size*y_size*z_size];
           for(int i=0; i<s.getSize();++i)
           {
           	Data[i] = s.Data[i];

           }
      	}
   	return *this;

   }


    ~Array()
    {
    	cleanup();
    }

   //===================================================================================================================
   //
   //  Inline Access Operators and Sizes
   //
   //===================================================================================================================

    void cleanup()
    {
    	if(Data != nullptr)
    	{
    	x_size = -1;
    	y_size = -1;
    	z_size = -1;
    	delete[] Data;
    	Data = nullptr;
    	}
    }

   void bounds_check(int i ,int j , int k ) const
   {
      std::stringstream msg;
      msg<<"BOUND ERROR: The given indexes: ("<<i<<","<<j<<","<<k<<") DOES NOT EXIST; Array size is:"<<z_size<<"x"<<y_size<<"x"<<x_size;
      ASSERT_MSG(i>=0,msg.str());
      ASSERT_MSG(j>=0,msg.str());
      ASSERT_MSG(k>=0,msg.str());
      ASSERT_MSG(i<z_size,msg.str());
      ASSERT_MSG(j<y_size,msg.str());
      ASSERT_MSG(k<x_size,msg.str());
   }

   // Operator() 1D
   inline T& operator ()(int i)
   {
      return(getData(i));
   }

   // Operator() 2D
   inline T& operator ()(int i,int j)
   {
      return(getData(i,j));

   }

   // Operator() 3D
   inline T& operator ()(int i, int j, int k)
   {
        return(getData(i,j,k));
   }

   //for const arrays

   // Operator() 1D
   inline const T& operator ()(int i) const
   {
      return(getData(i));
   }

   // Operator() 2D
   inline const T& operator ()(int i,int j) const
   {
      return(getData(i,j));

   }

   // Operator() 3D
   inline const T& operator ()(int i, int j, int k) const
   {
      return(getData(i,j,k));
   }


   inline T& getData(int i)
     {
        ASSERT_MSG(!isempty(),"TRYING TO ACCESS UNDEFINED ARRAY");

        #ifndef NDEBUG
        bounds_check(0,0,i);
        #endif
        return Data[i];
     }


     inline T& getData(int i,int j)
     {
        ASSERT_MSG(!isempty(),"TRYING TO ACCESS UNDEFINED ARRAY");

        #ifndef NDEBUG
        bounds_check(0,i,j);          //bounds check only if DEBUG mode is on
        #endif
        return Data[i*x_size + j];

     }

     inline T& getData(int i, int j, int k)
     {
          ASSERT_MSG(!isempty(),"TRYING TO ACCESS UNDEFINED ARRAY");

          #ifndef NDEBUG
          bounds_check(i,j,k);
          #endif
          return Data[i*x_size*y_size + j*x_size + k];
     }

     //for const arrays

     inline const T& getData(int i) const
     {
        ASSERT_MSG(!isempty(),"TRYING TO ACCESS UNDEFINED ARRAY");
        bounds_check(0,0,i);
        return Data[i];
     }

     inline const T& getData(int i,int j) const
     {
        ASSERT_MSG(!isempty(),"TRYING TO ACCESS UNDEFINED ARRAY");
        bounds_check(0,i,j);
        return Data[i*x_size+j];

     }

     inline const T& getData(int i, int j, int k) const
     {
     	 ASSERT_MSG(!isempty(),"TRYING TO ACCESS UNDEFINED ARRAY");
     	 bounds_check(i,j,k);
         return Data[i*x_size*y_size + j*x_size + k];
     }


   // initialize the whole array with a constant value
   void fill( T value )
   {
        for(int i=0; i< getSize();++i)
        {
   	     Data[i] = value;
        }
   }
   // only applicable for 2D; corner if needed to fill value in corner
   void fill_left_boundary( T value, bool corner=true )
   {
	   int corner_ = 0;
	   if(corner == true)
		   corner_= 0;
	   else
		   corner_ = 1;

        for(int i=corner_; i< getSize(1)- corner_;++i)
        {
   	     getData(i,0) = value;
        }
   }


   // only applicable for 2D; corner if needed to fill value in corner
   void fill_right_boundary( T value, bool corner=true )
   {
	   int corner_ = 0;
	   if(corner == true)
		   corner_= 0;
	   else
		   corner_ = 1;

        for(int i=corner_; i< getSize(1)- corner_;++i)
        {
   	     getData(i,getSize(0)-1) = value;
        }
   }

   // only applicable for 2D; corner if needed to fill value in corner
   void fill_bottom_boundary( T value, bool corner=true)
   {
	   int corner_ = 0;
	   if(corner == true)
		   corner_= 0;
	   else
		   corner_ = 1;

        for(int i=corner_; i< getSize(0)- corner_;++i)
        {
   	     getData(0,i) = value;
        }
   }

   // only applicable for 2D; corner if needed to fill value in corner
   void fill_top_boundary( T value, bool corner=true )
   {
	   int corner_ = 0;
	   if(corner == true)
		   corner_= 0;
	   else
		   corner_ = 1;

        for(int i=corner_; i< getSize(0)- corner_;++i)
        {
   	     getData(getSize(1)-1,i) = value;
        }
   }

   // only applicable for 2D; corner if needed to fill value in corner; only works if both array is of same dimension
     void fill_left_boundary( Array<T> value, bool corner=true )
     {
  	   int corner_ = 0;
  	   if(corner == true)
  		   corner_= 0;
  	   else
  		   corner_ = 1;
            std::cout<< "Size = " << getSize(1)-corner_ << std::endl;
          for(int i=corner_; i< getSize(1)- corner_;++i)
          {
     	     getData(i,0) = value(i,0);
          }
     }


     // only applicable for 2D; corner if needed to fill value in corner
     void fill_right_boundary( Array<T> value, bool corner=true )
     {
  	   int corner_ = 0;
  	   if(corner == true)
  		   corner_= 0;
  	   else
  		   corner_ = 1;

          for(int i=corner_; i< getSize(1)- corner_;++i)
          {
     	     getData(i,getSize(0)-1) = value(i,getSize(0)-1);
          }
     }

     // only applicable for 2D; corner if needed to fill value in corner
     void fill_bottom_boundary( Array<T> value, bool corner=true)
     {
  	   int corner_ = 0;
  	   if(corner == true)
  		   corner_= 0;
  	   else
  		   corner_ = 1;

          for(int i=corner_; i< getSize(0)- corner_;++i)
          {
     	     getData(0,i) = value(0,i);
          }
     }

     // only applicable for 2D; corner if needed to fill value in corner
     void fill_top_boundary( Array<T> value, bool corner=true )
     {
  	   int corner_ = 0;
  	   if(corner == true)
  		   corner_= 0;
  	   else
  		   corner_ = 1;

          for(int i=corner_; i< getSize(0)- corner_;++i)
          {
     	     getData(getSize(1)-1,i) = value(getSize(1)-1,i);
          }
     }


   void fill_inner(T value)
      {
   	   for(int j=1; j<getSize(1)-1; ++j)
   		   for(int i=1; i<getSize(0)-1; ++i)
   		   {
   			   getData(j,i) = value;
   		   }
      }

   // initialize the whole array with a constant function, that takes no arguments for eg:rand
   void fill(std::function<T()> func)
   {
        for(int i=0; i< getSize();++i)
        {
   	     Data[i] = func();
        }
   }

   // initialize the array with a function take 2 arguments(x,y), and evaluates the function on each of the grid point to fill in the array -2D
   void fill(std::function<T(int,int)> func)
   {
	   for(int j=0; j<getSize(1); ++j)
		   for(int i=0; i<getSize(0); ++i)
			  {
			   getData(j,i) = func(i,j);
			  }
   }

   // same as above but only in inner domain; ghost layer is excluded
   void fill_inner(std::function<T(int,int)> func)
   {
	   for(int j=1; j<getSize(1)-1; ++j)
		   for(int i=1; i<getSize(0)-1; ++i)
		   {
			   getData(j,i) = func(i,j);
		   }
   }

   //copy to ghost ; corners would be excluded
   //factor is the scaling (optional
   void copy_right_ghost(real scale =1, real shift = 0.0)
   {
	   int imax = getSize(0)-1;
	   for(int j=1; j<getSize(1)-1; ++j)
		   {
			   getData(j,imax) = scale*getData(j,imax-1) + shift;
		   }
   }

   void copy_left_ghost(real scale=1, real shift=0)
   {
	  for(int j=1; j<getSize(1)-1; ++j)
		   {
			   getData(j,0) = scale*getData(j,1) + shift;
		   }
   }

   void copy_top_ghost(real scale=1, real shift=0)
   {
	   int jmax = getSize(1)-1;
	   for(int i=1; i<getSize(0)-1; ++i)
		   {
			   getData(jmax,i) = scale*getData(jmax-1,i) + shift;
		   }
   }

   void copy_bottom_ghost(real scale=1, real shift =0)
   {
	   for(int i=1; i<getSize(0)-1; ++i)
		   {
			   getData(0,i) = scale*getData(1,i) + shift;
		   }
   }
   //wraps the array to periodic in x direction;
   //direction specifies copying direction:  +ve => EAST -> WEST  ;  -ve=> WEST->EAST
   // 0 direction =>nothing to do
   //vector::direction specify the sequence of direction
      void wrap_x(std::vector<int> direction)
  {
     int overlap = direction.size();
	 for(int over=0; over<overlap; ++over )
	  {
	   if(direction[over] > 0)
	   {
	       for(int j=0; j<getSize(1); ++j)
	       {

		        getData(j,over) = getData(j,(getSize(0)-1)-(overlap-over));
	       }
	   }
	   else if(direction[over]<0)
	   {
	 	   for(int j=0; j<getSize(1); ++j)
	 	   {
	 		    getData(j,(getSize(0)-1)-(overlap-over))=   getData(j,over) ;
	 	   }
	   }

	  }
   }
   //wraps the array to periodic in y direction
   //direction specifies copying direction:  +ve => NORTH -> SOUTH  ;  -ve => SOUTH->NORTH
   // 0 direction =>nothing to do
   void wrap_y( std::vector<int> direction)
   {
	 int overlap = direction.size();
	 for(int over=0; over<overlap; ++over)
	  {
	   if(direction[over]>0)
	   {
	    for(int i=0; i<getSize(0); ++i)
	    {
			   getData(over,i) = getData((getSize(1)-1)-(overlap-over),i);
	    }
	   }
	   else if(direction[over]<0)
	   {
	    for(int i=0; i<getSize(0); ++i)
	    {
			   getData((getSize(1)-1)-(overlap-over),i)= getData(over,i) ;
	    }
	   }

	  }

   }


// finds the maximum value of absolute values in an array, only for 2D
   real abs_max()
   {
	   real max = -1.0;
	   for(int i =0; i<getSize(); ++i)
		   {
		     max = std::max(fabs((Data[i])),max);
		   }
	   return max;
   }

   // finds the minimum value of absolute values in an array, only for 2D
      real abs_min() const
      {
   	   real min = std::numeric_limits<real>::max();
   	   for(int i =0; i<getSize(); ++i)
   		   {
   			   min = std::min(fabs(Data[i]),min);
   		   }
   	   return min;
      }

      // finds the maximum value of values in an array, only for 2D
         real max() const
         {
      	   real max = std::numeric_limits<real>::min();
      	   for(int i =0; i<getSize(); ++i)
      		   {
      		     max = std::max(Data[i],max);
      		   }
      	   return max;
         }

        // finds the minimum value of values in an array, only for 2D
           real min()
           {
            real min = std::numeric_limits<real>::max();
            for(int i =0; i<getSize(); ++i)
         	   {
         		   min = std::min(Data[i],min);
         	   }
            return min;
           }

      //finds sum of all values in array
      real sum() const
      {
    	  real sum = 0.0;
    	  for(int j=0;j<y_size;++j)
    		  for(int i=0; i<x_size;++i)
    		  {
    			  sum+=getData(j,i);
    		  }
    	  return sum;
      }


   // return total size of array
   int getSize() const
   {
   	if(!isempty())
   	{
   		return x_size*y_size*z_size;
   	}
   	else
   	{
   		return 0;
   	}
   }

   // return xSize for dimension==0, ySize for dimension==1 and zSize for dimension==2
   // other dimension values are not allowed
   int getSize(int dim) const
   {
   	std::stringstream msg;
   	if(!isempty())
   	{
   		if(dim == 0)
   			return x_size;
   		else if(dim ==1)
   			return y_size;
   		else if(dim == 2)
   			return z_size;
   		else
   		{
   			msg<<"The specified dimension "<<dim<<" does not exist, please pass in a value between 0-2";
   			ERROR(msg.str());
   			return -1;
   		}

   	}
   	else
   	{
   		return 0;
   	}
   }

//swaps arrays efficiently
   void swap(Array<T>& other)
   {
	   T *temp = this->Data;
	   this->Data = other.Data;
	   other.Data = temp;
   }


   //checks if Array is empty or not
   bool isempty() const
   {
     if(Data != nullptr)
   	   return false;
     else
       return true;
   }
   // Print the whole array ( for debugging purposes )

   // For 2D Arrays the positive x-coordinate goes to the right
   //                   positive y-coordinate goes upwards
   //      -> the line with highest y-value should be printed first
   void print() const
   {
	   if(isempty())
	   {
	   INFO("Nothing to be printed");
	   }
	   else
	   {
	   for(int i = 0; i<z_size; ++i)
	   {
		   std::cout<<std::endl;
		   for (int j=y_size-1; j>=0;--j)
		   {
			   for(int k=0; k<x_size;++k)
			   {
				  std::cout<< this->operator ()(i,j,k)<<" ";
			   }
                      std::cout<<std::endl;
		   }
	   }
	   }
   }


};


#endif //ARRAY_HH

