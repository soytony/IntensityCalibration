/***************************************************************************
 *            filesysTools.cpp
 *
 *  Fr 02 Mär 2007 23:14:08 CET
 *  Copyright 2007 Rainer Kümmerle
 *  Email rk@raikue.net
 ****************************************************************************/
#include "filesys_tools.h"
#include "macros.h"

#include <sys/stat.h>
#include <ctime>
#include <sys/types.h>

#include <cstdio>
#include <iostream>

#ifdef UNIX
  #include <wordexp.h>
  #include <dirent.h>
#endif //UNIX

#ifdef WINDOWS
  #include <Windows.h>
  #include <WinBase.h>
  typedef unsigned int uint;
#endif //windows

namespace Ais3dTools {

std::string getFileExtension(const std::string& filename)
{
  std::string::size_type lastDot = filename.find_last_of('.');
  if (lastDot != std::string::npos) 
    return filename.substr(lastDot + 1);
  else
    return "";
}

std::string getPureFilename(const std::string& filename)
{
  std::string::size_type lastDot = filename.find_last_of('.');
  if (lastDot != std::string::npos) 
    return filename.substr(0, lastDot);
  else
    return filename;
}

std::string getBasename(const std::string& filename)
{
  char dirDivider = '/';
  #if WINDOWS
    dirDivider = '\\';
  #endif
  std::string::size_type lastSlash = filename.find_last_of(dirDivider);
  if (lastSlash != std::string::npos) 
    return filename.substr(lastSlash + 1);
  else
    return filename;
}

std::string getDirname(const std::string& filename)
{
  char dirDivider = '/';
  #if WINDOWS
    dirDivider = '\\';
  #endif
  std::string::size_type lastSlash = filename.find_last_of(dirDivider);
  if (lastSlash != std::string::npos) 
    return filename.substr(0, lastSlash);
  else
    return "";
}

std::string changeFileExtension(const std::string& filename, const std::string& newExt, bool stripDot)
{
  std::string::size_type lastDot = filename.find_last_of('.');
  if (lastDot != std::string::npos) {
    if (stripDot)
      return filename.substr(0, lastDot) + newExt;
    else
      return filename.substr(0, lastDot + 1) + newExt;
  } else
    return filename;
}

bool fileExists(const char* filename)
{
  struct stat statInfo;
  return (stat(filename, &statInfo) == 0);
}

bool isRegularFile(const char* filename)
{
  #ifdef WINDOWS
    std::cerr << __PRETTY_FUNCTION__ << ": not implemented on windows" << std::endl;
    return false;
  #else
    struct stat statInfo;
    return (stat(filename, &statInfo) == 0 && S_ISREG(statInfo.st_mode));
  #endif 
}

bool isDirectory(const char* filename)
{
  #ifdef WINDOWS
    std::cerr << __PRETTY_FUNCTION__ << ": not implemented on windows" << std::endl;
    return false;
  #else
    struct stat statInfo;
    return (stat(filename, &statInfo) == 0 && S_ISDIR(statInfo.st_mode));
  #endif
}

bool isSymbolicLink(const char* filename)
{
#ifdef WINDOWS
  return false;
#else
  struct stat statInfo;
  return (lstat(filename, &statInfo) == 0 && S_ISLNK(statInfo.st_mode));
#endif
}

std::string getCurrentDateAsFilename()
{
#ifdef WINDOWS
  std::cerr << __PRETTY_FUNCTION__ << ": not implemented on windows" << std::endl;
  return std::string();
#else
  time_t t = time(NULL);
  const size_t dateStrSize = 1024;
  char dateStr[dateStrSize];
  if (strftime(dateStr, dateStrSize, "%Y%m%d_%H%M%S", localtime(&t)) == 0)
    fprintf(stderr, "Error (%s: %s) Date: %s\n", __FILE__, __func__, dateStr);
  return std::string(dateStr);
#endif
}

time_t getLastModificationDate(const char* filename)
{
  struct stat statInfo;
  if (stat(filename, &statInfo) == 0) {
    return statInfo.st_mtime;
  } else {
    return 0;
  }
}

time_t getLastAccessDate(const char* filename)
{
  struct stat statInfo;
  if (stat(filename, &statInfo) == 0) {
    return statInfo.st_atime;
  } else {
    return 0;
  }
}

time_t getLastStatusChangeDate(const char* filename)
{
  struct stat statInfo;
  if (stat(filename, &statInfo) == 0) {
    return statInfo.st_ctime;
  } else {
    return 0;
  }
}

#ifdef UNIX
bool createDirectory(const char* dirName, bool pub)
{
  bool status = true;
  status &= (mkdir(dirName, 0) == 0);
  if (pub)
    status &= (0 == chmod(dirName, // set directory to rwxrwxrwx
          S_IRUSR | S_IWUSR | S_IXUSR |
          S_IRGRP | S_IWGRP | S_IXGRP |
          S_IROTH | S_IWOTH | S_IXOTH ));
  else
    status &= (0 == chmod(dirName, // set directory to rwxr-xr-x
          S_IRUSR | S_IWUSR | S_IXUSR |
          S_IRGRP | S_IXGRP |
          S_IROTH | S_IXOTH ));
  return status;
}
#endif

#ifdef WINDOWS
bool createDirectory(const char* dirName, bool pub)
{
  std::cerr << __PRETTY_FUNCTION__ << ": not implemented on windows" << std::endl;
  bool status = true;
  //status &= (mkdir(dirName) == 0);
  return status;
}
#endif

off_t getFileSize(const char* filename)
{
  struct stat statInfo;
  if (stat(filename, &statInfo) == 0) {
    return statInfo.st_size;
  } else {
    return -1;
  }
}

std::vector<std::string> getDirectoryElements(const char* dir, bool onlyFiles)
{
  std::vector<std::string> dirEntries;
#ifdef UNIX
  DIR* curDir = opendir(dir);
  if(curDir == NULL) {
    return dirEntries;
  }

  struct dirent* dirEnt = readdir(curDir);
  while(dirEnt != NULL) {
    std::string name = dirEnt->d_name;

    if (name != "." && name != ".." && (!onlyFiles || dirEnt->d_type == DT_REG))
      dirEntries.push_back(name);

//#ifdef WINDOWS
//    if (name != "." && name != ".." && (!onlyFiles))
//      dirEntries.push_back(name);
//#endif //windows
    dirEnt = readdir(curDir);
  }
  closedir(curDir);
#endif //linux
#ifdef WINDOWS
  std::cerr << __PRETTY_FUNCTION__ << ": not implemented on windows" << std::endl;   
#endif
  return dirEntries;
}

std::vector<std::string> getFilesByPattern(const char* pattern)
{
  #ifdef UNIX
  std::vector<std::string> result;
  wordexp_t p;
  wordexp(pattern, &p, 0);
  result.reserve(p.we_wordc);
  for (size_t i = 0; i < p.we_wordc; ++i)
    result.push_back(p.we_wordv[i]);
  wordfree(&p);
  return result;
  #endif //UNIX

  
  #ifdef WINDOWS
  std::cerr << "WARNING: " << __PRETTY_FUNCTION__ << " not implemented on windows target" << std::endl;
  return std::vector<std::string>();
  #endif// WINDOWS
}

std::string expandFilename(const std::string& filename) {
  std::vector<std::string> filenames = getFilesByPattern(filename.c_str());
  if (filenames.empty())
    return "";
  else
    return filenames.back();
}

bool copyFile(const char* source, const char* dest)
{
  FILE* in = fopen(source, "r");
  FILE* out = fopen(dest, "w");
  if (in == NULL || out == NULL) {
    if (in)
      fclose(in);
    if (out)
      fclose(out);
    return false;
  }

  // copy over the data
  int ch = EOF;
  while((ch = fgetc(in)) != EOF) {
    int wrote = fputc(ch, out);
    if (wrote == EOF) {
      fclose(in);
      fclose(out);
      return false;
    }
  }

  // clean up
  fclose(in);
  fclose(out);
  return true;
}

}
