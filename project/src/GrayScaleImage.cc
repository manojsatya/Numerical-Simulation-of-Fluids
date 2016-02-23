//======================================================================================================================
/*!
 *  \file   GrayScaleImage.cpp
 *  \author Martin Bauer <martin.bauer@fau.de>
 */
//======================================================================================================================


// Local includes
#include "GrayScaleImage.hh"


// Extern
#include "lodepng.h"

#include <limits>


GrayScaleImage::GrayScaleImage( const std::string & pngFilename )
{
   unsigned int tmpWidth, tmpHeight;
   unsigned int error = lodepng::decode( image_, tmpWidth, tmpHeight, pngFilename, LCT_GREY, 8 );
   size_[0] = tmpWidth;
   size_[1] = tmpHeight;

   CHECK_MSG( error == 0 , "Error while loading PNG file: " << lodepng_error_text(error) );
}

void GrayScaleImage::save( const std::string & pngFilename )
{
   unsigned int error = lodepng::encode( pngFilename, image_,
                                         static_cast<unsigned int>( size_[0] ),
                                         static_cast<unsigned int>( size_[1] ),
                                         LCT_GREY, 8 );

   CHECK_MSG( error ==0 , "Error while loading PNG file: " << lodepng_error_text(error) );
}


real GrayScaleImage::operator() ( int x, int y ) const
{
   static const real maxVal = real( std::numeric_limits<unsigned char>::max() );
   const int yFlip = size_[1] - y - 1;

   return real ( getElement(x, yFlip) ) / maxVal;
}

//BUG in this function cleared
GrayScaleImage GrayScaleImage::getResizedImage( int newWidth, int newHeight ) const
{
   if ( newWidth == size_[0]  && newHeight == size_[1] )
      return *this;


   GrayScaleImage resizedImage;
   resizedImage.size_[0] = newWidth;
   resizedImage.size_[1] = newHeight;

   resizedImage.image_.resize( newWidth * newHeight );

   real scaleX = real( size_[0] ) / real( newWidth );
   real scaleY = real( size_[1] ) / real( newHeight );

   for( int y = 0; y <  newHeight; ++y )
      for( int x = 0; x <  newWidth; ++x )
      {
         real oldX = x * scaleX;
         real oldY = y * scaleY;
         int oldXi = static_cast<int>( oldX );
         int oldYi = static_cast<int>( oldY );
         int wrap_oldXi_east,wrap_oldYi_north;

         real xDiff = oldX - oldXi;
         real yDiff = oldY - oldYi;


         //we need a wrapper to avoid out of bounds
         if(oldXi+1>=size_[0])
        	 wrap_oldXi_east = oldXi;
         else
        	 wrap_oldXi_east = oldXi+1;

         if(oldYi+1>=size_[0])
        	 wrap_oldYi_north = oldYi;
         else
        	 wrap_oldYi_north = oldYi+1;


         // bilinear interpolation
         //this interpolation was giving OUT OF BOUNDS due to  oldXi + 1 and oldYi + 1; avoided with wrapped access
         resizedImage.getElement( x, y ) =
                  static_cast<unsigned char> (
                  (1 - xDiff) * (1 - yDiff ) * getElement( oldXi    , oldYi    ) +
                       xDiff  * (1 - yDiff ) * getElement( wrap_oldXi_east, oldYi    ) +
                  (1 - xDiff) *      yDiff   * getElement( oldXi    , wrap_oldYi_north ) +
                       xDiff  *      yDiff   * getElement( wrap_oldXi_east, wrap_oldYi_north ) );
      }

   return resizedImage;
}










