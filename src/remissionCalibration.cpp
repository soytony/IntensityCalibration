// This software is licenced under the Creative Commons
// (Attribution-NonCommercial-ShareAlike):
// http://creativecommons.org/licenses/by-nc-sa/4.0

// Author: Bastian Steder <steder@informatik.uni-freiburg.de>


//#include <qapplication.h>
#include <fstream>
#include "remissionCalibrationHelper.h"

//bool useGui = false;

RemissionCalibrationResult& calibration = RemissionCalibrationHelper::remissionCalibrationResult;
int& noOfRangeCells          = calibration.noOfCellsX,
   & noOfIncidenceAngleCells = calibration.noOfCellsY,
   & noOfIntensityCells      = calibration.noOfCellsZ;
float & outlierPercentage    = RemissionCalibrationHelper::parameters.outlierPercentage;
int mode = calibration.mode;

float fakeDatasetNoise = 0.05f;

void printUsage(const char* progName) {
  std::cout << "\n\nUsage: "<<progName<<" [options]\n\n"
       << "Options:\n"
       << "------------------------------\n"
       << "-g                          Use GUI\n"
       << "-d <directory>              Folder\n"
       << "-i <int>                    Number of G2O iterations (default: "<<RemissionCalibrationHelper::parameters.noOfIterationsForG2O<<")\n"
       << "-e <int>                    Number of EM iterations (default: "<<RemissionCalibrationHelper::parameters.noOfIterationsForEM<<")\n"
       << "-n                          Normalize the point clouds in the given folder (-d) with the given calibration file (-c)\n"
       << "-f                          Write the point clouds in the given folder (-d) with fake intensity values\n"
       << "-c <calibrationValuesFile>  Calibration values file\n"
       << "-x <int>                    Number of cells in x, determining in how many steps the range is divided "
       <<                              "(default "<<noOfRangeCells<<").\n"
       << "-y <int>                    Number of cells in y, determining in how many steps the incidence angle is divided "
       <<                              "(default "<<noOfIncidenceAngleCells<<")\n"
       << "-z <int>                    Number of cells in z, determining in how many steps the intensity is divided "
       <<                              "(default "<<noOfIntensityCells<<")\n"
       << "-m <1,2,3>                  mode - 1: full graph. 2: range&angle independency.\n"
       << "                                   3: range&angle dependency plus laser grouping.\n"
       << "                                   4: mode 3, followed by mode 2. (default: "<<mode<<")\n"
       << "-l <string>                 laser groups. Only applicable for mode 3. E.g.: \"1,{2,3,5,7},4-8,{4,6,101,103-107}\".\n"
       << "                                          All lasers not in the string get their own group.\n"
       << "-s <string>                 skip lasers. The listed lasers will be ignored. Same format as for -l.\n"
       << "-o <float>                  outlierPercentage (default: "<<outlierPercentage<<").\n"
       << "-h                          this help\n"
       << "\n\n";
}

int main(int argc, char** argv)
{
  setlocale(LC_NUMERIC, "C");

  std::string folder = ".";
  std::string calibrationValuesFileName = "";
  bool normalizePointClouds = false;
  bool writeFakePointClouds = false;
  std::string& laserGroupString = RemissionCalibrationHelper::laserGroupString;
  std::string& skipLasersString = RemissionCalibrationHelper::skipLasersString;
  
  // --------------------------------------
  // -----Parse Command Line Arguments-----
  // --------------------------------------
  for (char c; (c = getopt(argc, argv, "fnx:y:z:m:l:i:e:c:d:s:o:h")) != -1; ) {
    switch (c) {
      //case 'g':
        //useGui = true;
        //std::cout << "Using GUI.\n";
        //break;
      case 'd':
        folder = optarg;
        std::cout << "Folder \""<<folder<<"\" given.\n";
        break;
      case 'c':
        calibrationValuesFileName = optarg;
        std::cout << "Calibration values file \""<<calibrationValuesFileName<<"\" given.\n";
        break;
      case 'n':
        normalizePointClouds = true;
        break;
      case 'f':
        writeFakePointClouds = true;
        break;
      //case 'o':
        //std::cout << "Will use greater than one edges to prevent trivial optimization solutions.\n";
        //RemissionCalibrationHelper::parameters.useGreaterThanOneMethod = true;
        //break;
      case 'x':
        noOfRangeCells = atoi(optarg);
        std::cout << "No of cells in x (for range) set to "<<noOfRangeCells<<".\n";
        break;
      case 'y':
        noOfIncidenceAngleCells = atoi(optarg);
        std::cout << "No of cells in y (for incidence angle) set to "<<noOfIncidenceAngleCells<<".\n";
        break;
      case 'z':
        noOfIntensityCells = atoi(optarg);
        std::cout << "No of cells in z (for intensity) set to "<<noOfIntensityCells<<".\n";
        break;
      case 'i':
        RemissionCalibrationHelper::parameters.noOfIterationsForG2O = atoi(optarg);
        std::cout << "Settin number of G2O iterations to "<<RemissionCalibrationHelper::parameters.noOfIterationsForG2O<<".\n";
        break;
      case 'e':
        RemissionCalibrationHelper::parameters.noOfIterationsForEM = atoi(optarg);
        std::cout << "Settin number of EM iterations to "<<RemissionCalibrationHelper::parameters.noOfIterationsForEM<<".\n";
        break;
      case 'm':
        mode = atoi(optarg);
        if (mode==1)
          std::cout << "Mode 1: Full graph\n";
        else if (mode==2)
          std::cout << "Mode 2: Independency assumption for range and incidence angle\n";
        else if (mode==3)
          std::cout << "Mode 3: Independency assumption for range and incidence angle plus laser grouping\n";
        else if (mode==4)
          std::cout << "Mode 4: First mode 3, then mode 2\n";
        else {
          std::cerr << "\nInvalid mode selected ("<<mode<<").\n";
          return 1;
        }
        break;
      case 'l':
        laserGroupString = optarg;
        std::cout << "Laser group string: \""<<laserGroupString<<"\".\n";
        break;
      case 's':
        skipLasersString = optarg;
        std::cout << "Skipping lasers: \""<<skipLasersString<<"\".\n";
        break;
      case 'o':
        outlierPercentage = atof(optarg);
        std::cout << "Setting outlier percentage to "<<outlierPercentage<<".\n";
        break;
      case 'h':
        printUsage(argv[0]);
        exit(0);
    }
  }
  
  std::string mapCellsFilename = folder+"/remissionCalibrationMapCells.dat";
  
  //if (!calibrationValuesFileName.empty() && !normalizePointClouds) {
    //calibration.readCalibrationValuesFromFile(calibrationValuesFileName);
    //calibration.printCalibrationValues();
    //RemissionCalibrationHelper::CollectedPointsMap mapCells;
    //if (mapCells.readFromDisk(mapCellsFilename)) {
      //mapCells.setCellValuesToCalibratedMeans(calibration);
      //std::cout << "RMSE: "<<mapCells.getRMSE(calibration) << "\n";
      ////mapCells.removeOutlierCells(50.0, calibration);
      //std::cout << "Writing calibrated map cloud of size "<<mapCells.size()<<".\n";
      //mapCells.writeCalibratedMapCloudToDisk(calibration);
    //}
    //return 0;
  //}
  
  calibration.setNoOfCells(noOfRangeCells, noOfIncidenceAngleCells, noOfIntensityCells);

  if (mode==4)
    calibration.mode = 3;
  else
    calibration.mode = mode;
  
  RemissionCalibrationHelper::parameters.gridWeight = 1e4;
  RemissionCalibrationHelper::parameters.measurementWeight = 1.0;
  
  //QApplication* application = NULL;
  //if (useGui) {
    //application = new QApplication(argc,argv);
    //RemissionCalibrationHelper::qApplication = application;
    //application->processEvents();
  //}
   
  if (writeFakePointClouds) {
    std::cerr << "I will now write point clouds with fake intensities from the point clouds in folder \""
              << folder<<".\n";
    if (!RemissionCalibrationHelper::writeFakeDataset(folder, fakeDatasetNoise)) {
      std::cerr << "Something went wrong...\n";
      return 1;
    }
    return 0;
  }
  
  if (normalizePointClouds) {
    if (calibrationValuesFileName.empty()) {
      std::cerr << "I am supposed to normalize the intensities of the point clouds in folder \""
                << folder<<"\", but no calibration file is given. :-(\n";
      return 1;
    }
    std::cout << "I will now normalize the intensities of the point clouds in folder \""
              << folder<<"\" using the calibration file \""<<calibrationValuesFileName<<"\".\n";
    if (!RemissionCalibrationHelper::writeCalibratedClouds(folder, calibrationValuesFileName)) {
      std::cerr << "Something went wrong...\n";
      return 1;
    }
    return 0;
  }
  
  if (!calibrationValuesFileName.empty())
    calibration.readCalibrationValuesFromFile(calibrationValuesFileName);

  if (!RemissionCalibrationHelper::calculateCalibrationBasedOnMapCells(mapCellsFilename))
    std::cerr << "\n\nCalibration based on mapCells file \""<<mapCellsFilename<<"\" failed.\n\n";
  else {
    if (mode==4) {
      calibration.mode = 2;
      calibration.copyMode3ToMode2();
      RemissionCalibrationHelper::calculateCalibrationBasedOnMapCells(mapCellsFilename, true);
    }
  }
  return 0;
}
