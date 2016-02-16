#ifndef ARRAY_HH
#define ARRAY_HH

#include <vector>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include "Debug.hh"
#include "Types.hh"



//*******************************************************************************************************************
/*!  Array class for 1,2 and 3 dimensions
*
*    - all elements should be stored in a contiguous chunk of memory ( no vector<vector> ! )
*/
//*******************************************************************************************************************
template<class T>
class Array
{
public:
   // Constructors for 1D,2D and 3D
   Array( int xSize );
   Array( int xSize, int ySize );
   Array( int xSize, int ySize, int zSize );


   // Depending on your implementation you might need the following:
   // ~Array();
   // Array(const Array& s);
   // Array& operator= (const Array& s);


   // Access Operators for 1D, 2D and 3D
   inline T & operator () ( int i );
   inline T & operator () ( int i ,int j );
   inline T & operator () ( int i, int j, int k );

   // for const Arrays the following access operators are required
    inline const T & operator () ( int i ) const;
    inline const T & operator () ( int i ,int j ) const;
    inline const T & operator () ( int i, int j, int k ) const;



   // initialize the whole array with a constant value
   void fill( T value );


   // return total size of the array
   int getSize() const;

   // return xSize for dimension==0, ySize for dimension==1 and zSize for dimension==2
   // other dimension values are not allowed
   int getSize(int dimension ) const;


   // Print the whole array ( for debugging purposes )
   void print();

      
   // To get the absolute maximum element
   real get_abs_maxelement();
   

private:
    std::vector<T> matrix;// Object Array
	int imax,jmax,kmax; 

};


//===================================================================================================================
//
//  Inline Access Operators and Sizes
//
//===================================================================================================================


// Operator() 1D
template<class T>
inline T& Array<T>::operator ()(int i)
{
   //TODO
	CHECK_MSG(jmax==1 && kmax==1, "For 1D Simulation");
	CHECK_MSG(i >= 0 && i < imax, "i is bounded between 0(included) and imax");	
   	return matrix[i];
}

// Operator() 2D
template<class T>
inline T& Array<T>::operator ()(int i,int j)
{
   //TODO
   	CHECK_MSG(kmax==1, "For 2D Simulation");
	CHECK_MSG(i < imax && i >= 0 && j < jmax && j >= 0 , "i and j are bounded between 0(included) and respective max values");	
   	return matrix[i + imax*j];
}

// Operator() 3D
template<class T>
inline T& Array<T>::operator ()(int i, int j, int k)
{
   //TODO
   	CHECK_MSG(i >= 0 && i < imax && j >= 0 && j < jmax && k >= 0 && k < kmax, "i ,j ,k are bounded between 0(included) and respective max values");	
   	return matrix[i + imax*j + imax* jmax * k];
   
}

// Operator() 1D
template<class T>
inline const T& Array<T>::operator ()(int i) const
{
   //TODO
	CHECK_MSG(jmax==1 && kmax==1, "For 1D Simulation");
	CHECK_MSG(i >= 0 && i < imax, "i is bounded between 0(included) and imax");	
   	return matrix[i];
}

// Operator() 2D
template<class T>
inline const T& Array<T>::operator ()(int i,int j) const
{
   //TODO
   	CHECK_MSG(kmax==1, "For 2D Simulation");
	CHECK_MSG(i < imax && i >= 0 && j < jmax && j >= 0 , "i and j are bounded between 0(included) and respective max values");	
   	return matrix[i + imax*j];
}

// Operator() 3D
template<class T>
inline const T& Array<T>::operator ()(int i, int j, int k) const
{
   //TODO
   	CHECK_MSG(i >= 0 && i < imax && j >= 0 && j < jmax && k >= 0 && k < kmax, "i ,j ,k are bounded between 0(included) and respective max values");	
   	return matrix[i + imax*j + imax* jmax * k];
   
}

//#include "Array.cc"

#endif //ARRAY_HH

