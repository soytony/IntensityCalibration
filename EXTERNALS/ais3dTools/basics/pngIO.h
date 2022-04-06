/* \author Bastian Steder */

#ifndef AIS3TOOL_PNG_IO_H
#define AIS3TOOL_PNG_IO_H

#include <iostream>
#include <vector>

namespace Ais3dTools {

/**
 * \brief Capsulates libPNG functionalities
 */
class PngIO
{
  public:
    //-----PUBLIC STATIC FUNCTIONS-----
    //! Create a PNG from the input data - keys and values can be used to write meta data to the PNG
    static bool createPng(const unsigned char* imageData, int width, int height,
                          int bitDepth, int colorType, std::vector<unsigned char>& outputData,
                          const std::vector<std::string>* keys=NULL, const std::vector<std::string>* values=NULL);
    
    //! Create PNG from the 16Bit grayscale input data - useful to compress Openni depth images
    static bool create16BitGrayscalePng(const unsigned char* imageData, int width, int height, std::vector<unsigned char>& outputData,
                                        const std::vector<std::string>* keys=NULL, const std::vector<std::string>* values=NULL);
    
    //! Decompress a PNG - keys and values can be used to read meta data from the PNG
    static bool decompressPng(const unsigned char* data, std::vector<unsigned char>& outputData,
                              int& width, int& height, int& colorType, int& bitDepth,
                              std::vector<std::string>* keys=NULL, std::vector<std::string>* values=NULL);
    
    //! Same as above, but using a file from the disk
    static bool decompressPngFromDisk(const std::string& filename, std::vector<unsigned char>& outputData,
                                      int& width, int& height, int& colorType, int& bitDepth,
                                      std::vector<std::string>* keys=NULL, std::vector<std::string>* values=NULL);
    
  protected:
    //-----PROTECTED MEMBER VARIABLES-----
    //std::vector<unsigned char> outputData_;
};

} // namespace end

#endif
