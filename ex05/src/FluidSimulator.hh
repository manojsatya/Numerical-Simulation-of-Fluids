#ifndef __FLUID_SIMULATOR_H__
#define __FLUID_SIMULATOR_H__

#include <string>
#include <iostream>
#include <cmath>
#include "FileReader.hh"
#include "StaggeredGrid.hh"
#include "SORSolver.hh"

#include "VTKWriter.hh"

class FluidSimulator
{
  public:
      FluidSimulator( const FileReader & conf );

      /// Simulates a given time-length
      void simulate             ( real duration              );
      void simulateTimeStepCount( unsigned int nrOfTimeSteps );      


      // Getter functions for the internally stored StaggeredGrid
            StaggeredGrid & grid();
      const StaggeredGrid & grid() const;
      
      real getU(real x , real y);
      real getV(real x , real y);
      real getP(real x , real y);

  private:
      void computeFG();
      void composeRHS();
      void updateVelocities();
      void determineNextDT();
      void refreshBoundaries(); 
      void normalize();
      void initializeU();
      void printCourantNumber();
      void boundaryCorrection();

      StaggeredGrid grid_; 
      SORSolver solver_;
      real gx_;
      real gy_;
      real Re_;
      real gamma_;
      real dt_;
      size_t timesteps_;
      real tau_ ;
      size_t nf_;
      real xlength_ ;
      real ylength_ ;
      size_t outptinter_ ;
      

      std::string boundary_condition_S_ ;
      std::string boundary_condition_N_ ;
      std::string boundary_condition_W_ ;
      std::string boundary_condition_E_ ;
      std::string name_;
      std::string flowField_;

      real boundary_velocity_S_ ;
      real boundary_velocity_N_ ;
      real boundary_velocity_W_ ;
      real boundary_velocity_E_ ;		

};



#endif
