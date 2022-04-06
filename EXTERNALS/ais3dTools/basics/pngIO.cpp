/* \author Bastian Steder */

#include "pngIO.h"

#include <png.h>
#include <fstream>
#include <string.h>
#include "basics/macros.h"
#include "basics/timeutil.h"

void pngWriteDataHandler(png_structp pngPtr, png_bytep data, png_size_t length)
{
  //std::vector<char>& output = *(std::vector<char>*)png_ptr->io_ptr;
  std::vector<char>& output = *(std::vector<char>*)png_get_io_ptr(pngPtr);
  output.insert(output.end(), data, data+length);
}
void pngFlushHandler(png_structp pngPtr)
{
  (void)pngPtr; // No warning for not being used
}

struct InputData {
  const unsigned char* data;
  size_t pos;
};

void pngReadDataHandler(png_structp pngPtr, png_bytep data, png_size_t length) {
  InputData& input = *(InputData*)png_get_io_ptr(pngPtr);
  memcpy(data, input.data+input.pos, length);
  input.pos += length;
}

bool Ais3dTools::PngIO::createPng(const unsigned char* imageData, int width, int height,
                                  int bitDepth, int colorType, std::vector<unsigned char>& outputData,
                                  const std::vector<std::string>* keys, const std::vector<std::string>* values)
{
  //AIS3DTOOLS_MEASURE_FUNCTION_TIME;
  
  png_structp png_ptr = NULL;
  png_infop info_ptr = NULL;
  png_bytepp row_pointers = NULL;
  
  png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  if (png_ptr == NULL) {
    return false;
  }
  
  outputData.clear();
  png_set_write_fn(png_ptr, &outputData, pngWriteDataHandler, pngFlushHandler);
  
  info_ptr = png_create_info_struct(png_ptr);
  if (info_ptr == NULL) {
    png_destroy_write_struct(&png_ptr, NULL);
    return false;
  }
  
  if (setjmp(png_jmpbuf(png_ptr))) {
    png_free(png_ptr, row_pointers);
    png_destroy_write_struct(&png_ptr, &info_ptr);
    return false;
  }
  
  png_set_IHDR(png_ptr, info_ptr,
    width, height,
    bitDepth,      // Bit depth
    colorType,     // Color encoding
    PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT  // Default values
  );
  
  if (keys!=NULL && values!=NULL) {
    if (keys->size()!=values->size()) {
      std::cerr << __PRETTY_FUNCTION__<<": Error, keys->size()!=values->size()\n";
    }
    else if (!keys->empty()) {
      std::vector<png_text> texts;
      texts.resize(keys->size());
      for (size_t textIdx=0; textIdx<texts.size(); ++textIdx) {
        texts[textIdx].compression = -1;
        texts[textIdx].key  = const_cast<char*>(keys->at(textIdx).c_str());
        texts[textIdx].text = const_cast<char*>(values->at(textIdx).c_str());
        texts[textIdx].text_length = values->at(textIdx).size();
      }
      png_set_text(png_ptr, info_ptr, &texts[0], keys->size());
    }
  }


//[> png_get_text also returns the number of text chunks in *num_text <]
//extern PNG_EXPORT(png_uint_32,png_get_text) PNGARG((png_structp png_ptr,
   //png_infop info_ptr, png_textp *text_ptr, int *num_text));
//#endif
  
  row_pointers = (png_bytepp)png_malloc(png_ptr, height*sizeof(png_bytep));
  size_t rowSize = png_get_rowbytes(png_ptr,info_ptr);
  for (int y=0; y<height; y++)
    row_pointers[y] = (png_byte*)&imageData[y*rowSize];
  
  png_set_rows(png_ptr, info_ptr, row_pointers);
  png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);
  
  // Free allocated memory
  png_free(png_ptr, row_pointers);
  png_destroy_write_struct(&png_ptr, &info_ptr);
  
  //std::cout << "Created PNG of size "<<outputData.size()<<"\n";
  return true;
}

bool Ais3dTools::PngIO::create16BitGrayscalePng(const unsigned char* imageData, int width, int height, std::vector<unsigned char>& outputData,
                                                const std::vector<std::string>* keys, const std::vector<std::string>* values)
{
  return createPng(imageData, width, height, 16, PNG_COLOR_TYPE_GRAY, outputData, keys, values);
}

bool Ais3dTools::PngIO::decompressPngFromDisk(const std::string& filename, std::vector<unsigned char>& outputData,
                                              int& width, int& height, int& colorType, int& bitDepth,
                                              std::vector<std::string>* keys, std::vector<std::string>* values)
{
  std::ifstream pngFile(filename.c_str());
  if (!pngFile)
    return false;
  
  pngFile.seekg(0,std::ios_base::end);
  int fileSize = pngFile.tellg();
  pngFile.seekg(0,std::ios_base::beg);
  
  std::vector<unsigned char> data;
  data.resize(fileSize);
  pngFile.read((char*)&data[0], fileSize);
  pngFile.close();
  bool ret = decompressPng(&data[0], outputData, width, height, colorType, bitDepth, keys, values);
  //if (keys!=NULL && values!=NULL) {
    //if (keys->empty())
      //std::cout << "No text information stored in PNG.\n";
    //for (size_t keyIdx=0; keyIdx<keys->size(); ++keyIdx)
      //std::cout << keys->at(keyIdx)<<" = "<<values->at(keyIdx)<<"\n";
  //}
  return ret;
}

bool Ais3dTools::PngIO::decompressPng(const unsigned char* data, std::vector<unsigned char>& outputData,
                                      int& width, int& height, int& colorType, int& bitDepth,
                                      std::vector<std::string>* keys, std::vector<std::string>* values)
{
  //AIS3DTOOLS_MEASURE_FUNCTION_TIME;

  if (keys!=NULL)
    keys->clear();
  if (values!=NULL)
    values->clear();
  
  png_structp png_ptr;
  png_infop info_ptr;
  //int number_of_passes;
  png_bytep * row_pointers;
  
  //unsigned char header[8];    // 8 is the maximum size that can be checked
  
  InputData source;
  source.data = data;
  source.pos = 0;
  
  if (png_sig_cmp((unsigned char*) source.data, 0, 8)) {
    std::cerr << __PRETTY_FUNCTION__<<": Not a valid PNG header.\n";
    return false;
  }
  
  png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  if (png_ptr==NULL) {
    std::cerr << __PRETTY_FUNCTION__<<": Could not create read struct.\n";
    return false;
  }
  
  info_ptr = png_create_info_struct(png_ptr);
  if (info_ptr==NULL) {
    std::cerr << __PRETTY_FUNCTION__<<": Could not create info struct.\n";
    png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
    return false;
  }
  
  png_set_read_fn(png_ptr, (png_voidp)&source, pngReadDataHandler);
  
  if (setjmp(png_jmpbuf(png_ptr))) {
    std::cerr << __PRETTY_FUNCTION__<<": Something went wrong while reading the PNG information.\n";
    png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
    return false;
  }
  png_read_info(png_ptr, info_ptr);
  
  width  = png_get_image_width(png_ptr, info_ptr);
  height = png_get_image_height(png_ptr, info_ptr);
  colorType = png_get_color_type(png_ptr, info_ptr);
  bitDepth  = png_get_bit_depth(png_ptr, info_ptr);
  //std::cout << PVARC(width)<<PVARC(height)<<PVARC(colorType)<<PVARN(bitDepth);
  
  if (keys!=NULL && values!=NULL) {
    png_textp texts;
    int textsSize;
    png_get_text(png_ptr, info_ptr, &texts, &textsSize);
    for (size_t textIdx=0; textIdx<(size_t)textsSize; ++textIdx) {
      keys->push_back(texts[textIdx].key);
      values->push_back(texts[textIdx].text);
    }
  }
  
  //number_of_passes = png_set_interlace_handling(png_ptr);
  png_read_update_info(png_ptr, info_ptr);
  
  if (setjmp(png_jmpbuf(png_ptr))) {
    std::cerr << __PRETTY_FUNCTION__<<": Something went wrong while decompressing the PNG.\n";
    png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
    return false;
  }
  
  row_pointers = (png_bytepp)png_malloc(png_ptr, height*sizeof(png_bytep));
  size_t rowSize = png_get_rowbytes(png_ptr,info_ptr);
  outputData.resize(height*rowSize);
  for (int y=0; y<height; y++)
    row_pointers[y] = (png_byte*)&outputData[y*rowSize];
  
  png_read_image(png_ptr, row_pointers);
  
  png_free(png_ptr, row_pointers);
  png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
  
  return true;
}
