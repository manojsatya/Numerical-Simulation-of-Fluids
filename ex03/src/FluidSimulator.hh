#ifndef __FLUID_SIMULATOR_H__
#define __FLUID_SIMULATOR_H__


#include "FileReader.hh"
#include "StaggeredGrid.hh"
#include "SORSolver.hh"
#include <cmath>

class FluidSimulator
{
  public:
      FluidSimulator( const FileReader & conf );

      /// Simulates a given time-length
      void simulate             ( real duration              );
      void simulateTimeStepCount( unsigned int nrOfTimeSteps );
	
      void testSimulate();


      // Getter functions for the internally stored StaggeredGrid
            StaggeredGrid & grid();
      const StaggeredGrid & grid() const;

  private:
      void computeFG();
      StaggeredGrid grid_; 
      SORSolver solver_;
      real gx_;
      real gy_;
      real Re_;
      real gamma_;
      real dt_;
      size_t timesteps_;

};



#endif
