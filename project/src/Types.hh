#ifndef TYPES_HH
#define TYPES_HH


// This typedef makes it possible to switch between float and double accuracy
// please do not use "float" or "double" directly in your code, but use real instead
typedef double real;


// Enumeration of boundary conditions
typedef enum { NOSLIP, SLIP, OUTFLOW, PERIODIC } BCTYPE;

//Directions
typedef enum { NORTH, EAST, WEST, SOUTH, CENTER } DIRECTION;

typedef enum{SOLID, FLUID} CELL_TYPE;


#endif //TYPES_HH
