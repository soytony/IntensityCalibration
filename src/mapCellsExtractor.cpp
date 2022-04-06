// This software is licenced under the Creative Commons
// (Attribution-NonCommercial-ShareAlike):
// http://creativecommons.org/licenses/by-nc-sa/4.0

// Author: Bastian Steder <steder@informatik.uni-freiburg.de>


//#include <qapplication.h>
#include <fstream>
#include "remissionCalibrationHelper.h"


#define BACKWARD_HAS_DW 1

//bool useGui = false;

size_t approximateNumber = 1e6;
std::string& skipLasersString = RemissionCalibrationHelper::skipLasersString;

void printUsage(const char* progName) {
  std::cout << "\n\nUsage: "<<progName<<" [options]\n\n"
       << "Options:\n"
       << "------------------------------\n"
       << "-g                          Use GUI\n"
       << "-r <float>                  Resolution of the map (default: "<<RemissionCalibrationHelper::parameters.mapCellSize<<"m per cell)\n"
       << "-n <int>                    Minimum number of points for averaging (default: "
       <<                              RemissionCalibrationHelper::parameters.minNoOfMeasurementsForDataPointMean<<")\n"
       << "-s <string>                 Skip lasers. The listed lasers will be ignored. E.g.: \"1,{2,3,5,7},4-8,{4,6,101,103-107}\".\n"
       << "-d <directory>              Folder\n"
       << "-h                          This help\n"
       << "\n\n";
}

int main(int argc, char** argv)
{
  setlocale(LC_NUMERIC, "C");
  
  std::string folder = ".";
  
  // --------------------------------------
  // -----Parse Command Line Arguments-----
  // --------------------------------------
  for (char c; (c = getopt(argc, argv, "r:n:s:d:h")) != -1; ) {
    switch (c) {
      //case 'g':
        //useGui = true;
        //std::cout << "Using GUI.\n";
        //break;
      case 'r':
        RemissionCalibrationHelper::parameters.mapCellSize = atof(optarg);
        std::cout << "Setting map resolution to "<<RemissionCalibrationHelper::parameters.mapCellSize<<"m per cell.\n";
        break;
      case 'n':
        RemissionCalibrationHelper::parameters.minNoOfMeasurementsForDataPointMean = atoi(optarg);
        std::cout << "Will use a minimum of "
                  << RemissionCalibrationHelper::parameters.minNoOfMeasurementsForDataPointMean<<" "
                  << "points for averaging.\n";
        break;
      case 'd':
        folder = optarg;
        std::cout << "Folder \""<<folder<<"\" given.\n";
        break;
      case 's':
        skipLasersString = optarg;
        std::cout << "Skipping lasers: \""<<skipLasersString<<"\".\n";
        break;
      case 'h':
        printUsage(argv[0]);
        exit(0);
    }
  }
  
  //QApplication* application = NULL;
  //if (useGui) {
    //application = new QApplication(argc,argv);
    //RemissionCalibrationHelper::qApplication = application;
    //application->processEvents();
  //}
  
  std::string mapCellsFilename = folder+"/remissionCalibrationMapCells.dat";
  
  RemissionCalibrationHelper::CollectedPointsMap mapCells;
  if (RemissionCalibrationHelper::calculateMapCells(folder, mapCells, mapCellsFilename))
    return 0;
  else
    return 1;
}
