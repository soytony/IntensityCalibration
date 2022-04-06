// This software is licenced under the Creative Commons
// (Attribution-NonCommercial-ShareAlike):
// http://creativecommons.org/licenses/by-nc-sa/4.0

// Author: Bastian Steder <steder@informatik.uni-freiburg.de>


#include "remissionCalibrationHelper.h"

std::string calibrationValuesFileName = "";
float errorStepSize = 0.1;
std::string& skipLasersString = RemissionCalibrationHelper::skipLasersString;

void printUsage(const char* progName) {
  std::cout << "\n\nUsage: "<<progName<<" [options]\n\n"
       << "Options:\n"
       << "------------------------------\n"
       << "-c <calibrationValuesFile>  Calibration values file\n"
       << "-e <float>                  Error step size (default: "<<errorStepSize<<")\n"
       << "-d <directory>              Folder\n"
       << "-s <string>                 skip lasers. The listed lasers will be ignored. E.g.: \"1,{2,3,5,7},4-8,{4,6,101,103-107}\".\n"
       << "-h                          this help\n"
       << "\n\n";
}

int main(int argc, char** argv)
{
  setlocale(LC_NUMERIC, "C");
  
  std::string folder = ".";
  
  // --------------------------------------
  // -----Parse Command Line Arguments-----
  // --------------------------------------
  for (char c; (c = getopt(argc, argv, "c:d:e:s:h")) != -1; ) {
    switch (c) {
      case 'c':
        calibrationValuesFileName = optarg;
        std::cout << "Calibration values file \""<<calibrationValuesFileName<<"\" given.\n";
        break;
      case 'e':
        errorStepSize = atof(optarg);
        std::cout << "Setting error step size to "<<errorStepSize<<".\n";
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
  
  std::string mapCellsFilename = folder+"/remissionCalibrationMapCells.dat";
  
  if (calibrationValuesFileName.empty()) {
    printUsage(argv[0]);
    std::cerr << "No calibration file given. :-(\n\n";
    return 1;
  }
  
  RemissionCalibrationResult& calibration = RemissionCalibrationHelper::remissionCalibrationResult;
  calibration.readCalibrationValuesFromFile(calibrationValuesFileName);
  //calibration.printCalibrationValues();
  std::cout << "\nRange cells: "    <<calibration.noOfCellsX<<"\n"
            <<   "Angle cells: "    <<calibration.noOfCellsY<<"\n"
            <<   "Intensity cells: "<<calibration.noOfCellsZ<<"\n\n";
  
  RemissionCalibrationHelper::CollectedPointsMap mapCells;
  if (!mapCells.readFromDisk(mapCellsFilename)) {
    printUsage(argv[0]);
    std::cerr << "Could not read map file \""<<mapCellsFilename<<"\". :-(\n\n";
    return 1;
  }
  
  std::cout << "No of data points in map: "<<mapCells.getNoOfDataPoints()<<"\n";
  
  std::ofstream file;
  std::stringstream ss;
  
  for (int step=0; step<2; ++step) {
    if (step==0) {
      std::cout << "\n---------------------------------\n"
                <<   "----------INITIAL GUESS----------\n"
                <<   "---------------------------------\n";
      RemissionCalibrationHelper::setCalibrationToInitialGuess(mapCells);
    }
    else {
      std::cout << "\n------------------------------\n"
                <<   "----------CALIBRATED----------\n"
                <<   "------------------------------\n";
      calibration.readCalibrationValuesFromFile(calibrationValuesFileName);
      
      mapCells.setCellValuesToCalibratedMeans(calibration);
      calibration.normalizeCalibrationFactors(1.0f / mapCells.getAverageCalibratedIntensity(calibration));

    }
    mapCells.setCellValuesToCalibratedMeans(calibration);
    std::cout << "Calibrated RMSE: "<<mapCells.getRMSE(calibration) << "\n";
    //mapCells.removeOutlierCells(50.0, calibration);
    
    std::map<int,int> errorDistribution = mapCells.getErrorDistribution(calibration, errorStepSize);
    for (std::map<int,int>::const_iterator it=errorDistribution.begin(); it!=errorDistribution.end(); ++it) {
      //std::cout << it->first*errorStepSize << ": "<<it->second<<"\n";
      //file << it->first*errorStepSize << " " << it->second << "\n";
    }
    //file.close();
    
    ss.str("");
    ss << "errorDistribution";
    if (step==0)
      ss << "_initialGuess.dat";
    else
      ss << "_calibratedMode"<<calibration.mode<<".dat";
    file.open(ss.str().c_str());
    std::map<int,float> errorPercentages = mapCells.getErrorPercentages(calibration);
    for (std::map<int,float>::const_iterator it=errorPercentages.begin(); it!=errorPercentages.end(); ++it) {
      //std::cout << it->first << "%: "<<it->second<<"\n";
      file << it->first << " "<<it->second<<"\n";
    }
    file.close();
    
    if (step==0)
      file.open("medianError_initial.txt");
    else
      file.open("medianError.txt");
    std::cout << "Median error: "<<errorPercentages[50] << "\n";
    file << round(errorPercentages[50]*1e3)/1e3 << "\n";
    file.close();
  }

  //if (calibration.mode==3) {
    //file.open("laserCalib.tex");
    //file.precision(3);
    //for (size_t laserGroupIdx=0; laserGroupIdx<calibration.calibrationValuesMode3.size(); ++laserGroupIdx) {
      //const RemissionCalibrationResult::LaserGroup& lg = calibration.calibrationValuesMode3[laserGroupIdx];
      //for (std::map<int, float>::const_iterator it=lg.calibrationValuesLaser.begin(); it!=lg.calibrationValuesLaser.end(); ++it)
        //file << (it==lg.calibrationValuesLaser.begin()&&laserGroupIdx==0 ? "" : " & ") << it->second;
    //}
    //file << "\n";
    //file.close();
  //}

  std::cout << "Average calibrated intensity: "<<mapCells.getAverageCalibratedIntensity(calibration)<<"\n";
  
  std::cout << "Writing gnuplot data files for the calibration.\n";
  calibration.saveGnuplotDataFiles();
  
  std::cout << "Writing calibrated map cloud of size "<<mapCells.size()<<".\n";
  mapCells.writeCalibratedMapCloudToDisk(calibration);
  
  return 0;
}
