#include "StaggeredGrid.hh"

StaggeredGrid::StaggeredGrid ( int jmax, int imax, real dx, real dy ) : p_(Array<real> (jmax+2,imax+2)),rhs_(Array<real> (jmax+2,imax+2)),
u_(Array<real> (jmax+2,imax+1)),v_(Array<real> (jmax+1,imax+2)),f_(Array<real>(jmax+2,imax+1)), g_(Array<real>(jmax+1,imax+2)),
k_old_(Array<real> (jmax+2,imax+2)),e_old_(Array<real> (jmax+2,imax+2)),k_new_(Array<real> (jmax+2,imax+2)),e_new_(Array<real> (jmax+2,imax+2)),
nu_t_(jmax+2,imax+2),flag_(jmax+2,imax+2),delta_(jmax+2,imax+2),dx_(dx),dy_(dy),imax_(imax),jmax_(jmax)
{
	flag_.fill(FLUID);
// manual constructor
//ghost layer considered
}

StaggeredGrid::StaggeredGrid (FileReader & configuration)
{
  real xSize = configuration.getParameter("xlength");
  real ySize = configuration.getParameter("ylength");
  imax_   = configuration.getParameter("imax");
  jmax_   = configuration.getParameter("jmax");

  dx_        = static_cast<real> (xSize) /imax_;
  dy_        = static_cast<real> (ySize)/jmax_;
  std::cout<<"dx"<<dx_;

  p_         = Array<real> (jmax_+2,imax_+2);
  rhs_       = Array<real> (jmax_+2,imax_+2);
  u_         = Array<real> (jmax_+2,imax_+1);
  v_		 = Array<real> (jmax_+1,imax_+2);
  f_ 		 = Array<real> (jmax_+2,imax_+1); //really only f_(jmax,imax+1) required but for convenience, so ghost layer on top and bottom
  g_		 = Array<real> (jmax_+1,imax_+2); //really only f_(jmax,imax+1) required but for convenience, so ghost layer on right and left
  k_old_     = Array<real> (jmax_+2,imax_+2);//for turbulence
  e_old_     = Array<real> (jmax_+2,imax_+2);//for turbulence
  k_new_     = Array<real> (jmax_+2,imax_+2);//for turbulence
  e_new_     = Array<real> (jmax_+2,imax_+2);//for turbulence
  nu_t_      = Array<real> (jmax_+2,imax_+2);//for turbulent viscosity

  flag_      = Array<CELL_TYPE> (jmax_+2,imax_+2);
  delta_     = Array<real> (jmax_+2,imax_+2);//for storing nearest wall distance

  flag_.fill(FLUID); //not required but still for readability
}

void StaggeredGrid::setCellToObstacle(const int y, const int x)
{
	flag_(y,x) = SOLID;
}

//sets rectangle obstacle
//@param x1,y1 : bottom-left coordinates
//@param x2,y2 : top-right coordinates
void StaggeredGrid::createRectangle(real x1, real y1, real x2, real y2)
{
  real x_origin = x1<x2 ? x1:x2;
  real y_origin = y1<y2 ? y1:y2;

  int x_start = static_cast<int>(x_origin/dx_);
  int y_start = static_cast<int>(y_origin/dy_);
  int length = static_cast<int>( std::abs(x2-x1) / dx_);
  int height = static_cast<int>( std::abs(y2-y1) / dy_);

  x_start += 1; //Exclude ghost layer
  y_start += 1;


  CHECK_MSG(x_start>0 && x_start<flag_.x_size,"The dimension of obstacle exceeds the domain size");
  CHECK_MSG(y_start>0 && y_start<flag_.x_size,"The dimension of obstacle exceeds the domain size");
  CHECK_MSG(x_start+length>0 && x_start+length<flag_.x_size,"The dimension of obstacle exceeds the domain size");
  CHECK_MSG(y_start+height>0 && y_start+height<flag_.x_size,"The dimension of obstacle exceeds the domain size");


  for(int y = y_start; y<y_start+height; ++y)
	  for(int x = x_start; x<x_start+length; ++x )
		  setCellToObstacle(y,x);

}

void StaggeredGrid::createCircle (real x, real y, real r)
{
	CHECK_MSG(r>0,"Radius of circle should be greater than 0");
	std::cout<<"\n"<<r<<"\n";
    int x_origin = static_cast<int> (x/dx_)+1; //+1 to exclude ghost layer
    int y_origin = static_cast<int> (y/dy_)+1;
    int radius_x   = static_cast<int> (r/dx_); //radius in x and y might be mapped differently, depending on discretisation
    int radius_y   = static_cast<int> (r/dy_);
    std::cout<<"\n"<<radius_y<<"\n";

    int x_start = x_origin-radius_x;
    int x_end   = x_origin+radius_x;
    int y_start = y_origin-radius_y;
    int y_end   = y_origin+radius_y;

    CHECK_MSG(x_start>0 && x_start<flag_.x_size,"The dimension of obstacle exceeds the domain size");
    CHECK_MSG(y_start>0 && y_start<flag_.x_size,"The dimension of obstacle exceeds the domain size");
	CHECK_MSG(x_end>0 && x_end<flag_.x_size,"The dimension of obstacle exceeds the domain size");
	CHECK_MSG(y_end>0 && y_end<flag_.x_size,"The dimension of obstacle exceeds the domain size");
	CHECK_MSG(radius_x > 0 && radius_y > 0,"Discretisation too coarse for the specified radius of circle");

	for(int y= y_start; y<=y_end ; ++y)
    	for(int x = x_start; x<=x_end ; ++x)
    	{
    		if( ( (x-x_origin)*(x-x_origin)/(radius_x*radius_x) ) + ( (y-y_origin)*(y-y_origin)/(radius_y*radius_y)) <= 1) //equation of ellipse
    		{
    			setCellToObstacle(y,x);
    		}
    	}

}

void StaggeredGrid::load_image(const std::string & pngFilename)
{
	GrayScaleImage img(pngFilename );

	img = img.getResizedImage(flag_.x_size,flag_.y_size); //resize the image to our domain size

	for(int y=0; y<img.height(); ++y)
        {
		for(int x=0; x<img.width(); ++x)
		{

			if(img(x,y) < 0.5)
			{
				setCellToObstacle(y,x);
			}
                }
		
         }
}

void StaggeredGrid::save_image(const std::string & pngFilename)
{

	GrayScaleImage img2(pngFilename );


	img2 = img2.getResizedImage(flag_.x_size,flag_.y_size);


	for(int y=0; y<flag_.y_size; ++y)
	{
		for(int x=0; x<flag_.x_size; ++x)
		{
			if(!isFluid(y,x))
		        img2.getElement(x,flag_.y_size-y-1) = 0; //flip and store
			else
				img2.getElement(x,flag_.y_size-y-1) = 255;//flip and store

		}


	}

	img2.save(pngFilename);
}
