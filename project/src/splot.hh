#ifndef SPLOT_HH
#define SPLOT_HH

#include "Debug.hh"
#include "Array.hh"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

void splot(Array<real> &X, Array<real> &Y, Array<real> &Z )
{
	PRG_LEVEL("Writing Plotting Data");
	ASSERT_MSG(X.getSize(1)==1,"SPLOT :ERROR in plotting X should be 1D");
	ASSERT_MSG(Y.getSize(1)==1,"SPLOT: ERROR in plotting Y should be 1D");
	ASSERT_MSG((X.getSize(0)==Z.getSize(0)&&(Y.getSize(0)==Z.getSize(1))),"SPLOT: X,Y,Z not compatible, please check the sizes");
	std::fstream file;
	file.open("solution.dat",std::fstream::out);
	for(int i=0; i<Z.getSize(1);++i)
	{
		for(int j=0; j<Z.getSize(0);++j)
		{
			file<<X(j)<<" "<<Y(i)<<" "<<Z(i,j)<<"\n";
		}
		file<<"\n";
	}
	INFO("Plotting Data Written; Please run the following command to visualize: 'gnuplot splot.plt'.A png file would be saved in your current directory");
	file.close();
}

void splot(real x_begin , real h_x, real y_begin,  real h_y, Array<real> &Z, std::string  filename="solution.dat")
{
	PRG_LEVEL("Writing Plotting Data");
	std::fstream file;
	std::stringstream out;
	out<<"plots/"<<filename<<".dat";
	file.open(out.str(),std::fstream::out);
	for(int i=0; i<Z.getSize(1);++i)
	{
		for(int j=0; j<Z.getSize(0);++j)
		{
			file<<x_begin+j*h_x<<" "<<y_begin+i*h_y<<" "<<Z(i,j)<<"\n";
		}
		file<<"\n";
	}
	file.close();
}



#endif //SPLOT
