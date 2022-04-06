// This software is licenced under the Creative Commons
// (Attribution-NonCommercial-ShareAlike):
// http://creativecommons.org/licenses/by-nc-sa/4.0

// Author: Bastian Steder <steder@informatik.uni-freiburg.de>


#include "remissionCalibrationResult.h"

#include <fstream>
#include <cstdlib>
#include "ais3dTools/basics/macros.h"

RemissionCalibrationResult::RemissionCalibrationResult() : mode(3) {
  setRangeStepFactor();  // Set default value
  setNoOfCells();        // Set default value
}

void RemissionCalibrationResult::updateDependingValues() {
  incidenceAngleToCellFactor = float(noOfCellsY-1)/(0.5f*M_PI);
  cellToIncidenceAngleFactor = (noOfCellsY>1 ? 1.0f/incidenceAngleToCellFactor : 0.0);
}

void RemissionCalibrationResult::setNoOfCells(int newNoOfCellsX, int newNoOfCellsY, int newNoOfCellsZ) {
  noOfCellsX = newNoOfCellsX;
  noOfCellsY = newNoOfCellsY;
  noOfCellsZ = newNoOfCellsZ;
  updateDependingValues();
}

void RemissionCalibrationResult::saveCalibrationValues(const std::string& filename) const {
  std::ofstream calibrationValuesFile(filename.c_str());
  calibrationValuesFile << mode<<" "<<rangeStepFactor<<" "<<noOfCellsX<<" "<<noOfCellsY<<" "<<noOfCellsZ<<"\n";
  if (mode==1) {  // full graph
    for (std::map<int, std::vector<float> >::const_iterator it=calibrationValuesMode1.begin(); it!=calibrationValuesMode1.end(); ++it) {
      int laser = it->first;
      const std::vector<float>& calibrationValues = it->second;
      calibrationValuesFile << laser << "\n"
                            << maxUncalibratedIntensities.find(laser)->second<<" "
                            << minMaxRanges.find(laser)->second.first<<" "<<minMaxRanges.find(laser)->second.second<<"\n";
      for (size_t calibrationValueIdx=0;  calibrationValueIdx<calibrationValues.size(); ++calibrationValueIdx) {
        calibrationValuesFile << " "<<calibrationValues[calibrationValueIdx];
      }
      calibrationValuesFile << "\n";
    }
  }
  else if (mode==2) {  // no full graph
    for (CalibIndependentPerLaser::const_iterator it=calibrationValuesMode2.begin(); it!=calibrationValuesMode2.end(); ++it) {
      int laser = it->first;
      const CalibIndependent& c = calibrationValuesMode2.find(laser)->second;
      calibrationValuesFile << laser << "\n"
                            << maxUncalibratedIntensities.find(laser)->second<<" "
                            << minMaxRanges.find(laser)->second.first<<" "<<minMaxRanges.find(laser)->second.second<<"\n";
      for (size_t calibrationValueIdx=0;  calibrationValueIdx<c.calibrationValuesRange.size(); ++calibrationValueIdx)
        calibrationValuesFile << " "<<c.calibrationValuesRange[calibrationValueIdx];
      calibrationValuesFile << "\n";
      for (size_t calibrationValueIdx=0;  calibrationValueIdx<c.calibrationValuesIncidenceAngle.size(); ++calibrationValueIdx)
        calibrationValuesFile << " "<<c.calibrationValuesIncidenceAngle[calibrationValueIdx];
      calibrationValuesFile << "\n";
      for (size_t calibrationValueIdx=0;  calibrationValueIdx<c.calibrationValuesIntensity.size(); ++calibrationValueIdx)
        calibrationValuesFile << " "<<c.calibrationValuesIntensity[calibrationValueIdx];
      calibrationValuesFile << "\n";
    }
  }
  else {  // no full graph with laser groups
    for (size_t laserGroupIdx=0; laserGroupIdx<calibrationValuesMode3.size(); ++laserGroupIdx) {
      calibrationValuesFile << "{\n";
      const LaserGroup& lg = calibrationValuesMode3[laserGroupIdx];
      for (size_t calibrationValueIdx=0;  calibrationValueIdx<lg.calibrationValuesRange.size(); ++calibrationValueIdx)
        calibrationValuesFile << " "<<lg.calibrationValuesRange[calibrationValueIdx];
      calibrationValuesFile << "\n";
      for (size_t calibrationValueIdx=0;  calibrationValueIdx<lg.calibrationValuesIncidenceAngle.size(); ++calibrationValueIdx)
        calibrationValuesFile << " "<<lg.calibrationValuesIncidenceAngle[calibrationValueIdx];
      calibrationValuesFile << "\n";
      for (size_t calibrationValueIdx=0;  calibrationValueIdx<lg.calibrationValuesIntensity.size(); ++calibrationValueIdx)
        calibrationValuesFile << " "<<lg.calibrationValuesIntensity[calibrationValueIdx];
      calibrationValuesFile << "\n";
      for (std::map<int, float>::const_iterator it=lg.calibrationValuesLaser.begin(); it!=lg.calibrationValuesLaser.end(); ++it)
        calibrationValuesFile << " "<<it->first<<" "<<it->second<<" "
                              << maxUncalibratedIntensities.find(it->first)->second<<" "
                              << minMaxRanges.find(it->first)->second.first<<" "<<minMaxRanges.find(it->first)->second.second<<"\n";
      calibrationValuesFile << "}\n";
    }
  }
  calibrationValuesFile.close();
}


bool RemissionCalibrationResult::readCalibrationValuesFromFile(const std::string& fileName) {
  calibrationValuesMode1.clear();
  calibrationValuesMode2.clear();
  calibrationValuesMode3.clear();
  maxUncalibratedIntensities.clear();
  minMaxRanges.clear();
  
  std::ifstream calibrationValuesFile(fileName.c_str());
  if (calibrationValuesFile) {
    calibrationValuesFile >> mode >> rangeStepFactor >> noOfCellsX >> noOfCellsY >> noOfCellsZ;
    updateDependingValues();
    
    if (mode==1) { // full graph
      while (true) {
        int laser;
        calibrationValuesFile >> laser;
        if (!calibrationValuesFile)
          break;
        //std::cout << PVARN(laser);
        float& maxUncalibratedIntensity = maxUncalibratedIntensities[laser];
        std::pair<float, float>& minMaxRange = minMaxRanges[laser];
        calibrationValuesFile >> maxUncalibratedIntensity >> minMaxRange.first >> minMaxRange.second;

        std::vector<float>& calibrationValues = calibrationValuesMode1[laser];
        calibrationValues.resize(noOfCellsX*noOfCellsY);
        for (size_t calibrationValueIdx=0; calibrationValueIdx<calibrationValues.size(); ++calibrationValueIdx)
          calibrationValuesFile >> calibrationValues[calibrationValueIdx];
      }
    }
    else if (mode==2) { // no full graph
      while (true) {
        int laser;
        calibrationValuesFile >> laser;
        if (!calibrationValuesFile)
          break;
        //std::cout << PVARN(laser);
        float& maxUncalibratedIntensity = maxUncalibratedIntensities[laser];
        std::pair<float, float>& minMaxRange = minMaxRanges[laser];
        calibrationValuesFile >> maxUncalibratedIntensity >> minMaxRange.first >> minMaxRange.second;
        
        CalibIndependent& c = calibrationValuesMode2[laser];
        c.calibrationValuesRange.resize(noOfCellsX);
        c.calibrationValuesIncidenceAngle.resize(noOfCellsY);
        c.calibrationValuesIntensity.resize(noOfCellsZ);
        for (size_t calibrationValueIdx=0;  calibrationValueIdx<c.calibrationValuesRange.size(); ++calibrationValueIdx)
          calibrationValuesFile >> c.calibrationValuesRange[calibrationValueIdx];
        for (size_t calibrationValueIdx=0;  calibrationValueIdx<c.calibrationValuesIncidenceAngle.size(); ++calibrationValueIdx)
          calibrationValuesFile >> c.calibrationValuesIncidenceAngle[calibrationValueIdx];
        for (size_t calibrationValueIdx=0;  calibrationValueIdx<c.calibrationValuesIntensity.size(); ++calibrationValueIdx)
          calibrationValuesFile >> c.calibrationValuesIntensity[calibrationValueIdx];
      }
    }
    else {  // no full graph with laser groups
      while (true) {
        std::string groupStart;
        calibrationValuesFile >> groupStart;
        if (groupStart!="{")
          break;
        calibrationValuesMode3.push_back(LaserGroup());
        LaserGroup& lg = calibrationValuesMode3.back();
        lg.calibrationValuesRange.resize(noOfCellsX);
        for (size_t calibrationValueIdx=0;  calibrationValueIdx<lg.calibrationValuesRange.size(); ++calibrationValueIdx)
          calibrationValuesFile >> lg.calibrationValuesRange[calibrationValueIdx];
        lg.calibrationValuesIncidenceAngle.resize(noOfCellsY);
        for (size_t calibrationValueIdx=0;  calibrationValueIdx<lg.calibrationValuesIncidenceAngle.size(); ++calibrationValueIdx)
          calibrationValuesFile >> lg.calibrationValuesIncidenceAngle[calibrationValueIdx];
        lg.calibrationValuesIntensity.resize(noOfCellsZ);
        for (size_t calibrationValueIdx=0;  calibrationValueIdx<lg.calibrationValuesIntensity.size(); ++calibrationValueIdx)
          calibrationValuesFile >> lg.calibrationValuesIntensity[calibrationValueIdx];
        int laser;
        float calibrationValue;
        while (calibrationValuesFile) {
          std::string token;
          calibrationValuesFile >> token;
          if (token == "}")
            break;
          laser = atoi(token.c_str());
          calibrationValuesFile >> calibrationValue
                                >> maxUncalibratedIntensities[laser]
                                >> minMaxRanges[laser].first >> minMaxRanges[laser].second;
          //std::cout << PVARC(laser)<<PVARN(calibrationValue);
          lg.calibrationValuesLaser[laser] = calibrationValue;
        }
      }
      //for (size_t laserGroupIdx=0; laserGroupIdx<calibrationValuesMode3.size(); ++laserGroupIdx) {
        //std::cout << "Laser group "<<laserGroupIdx+1<<"\n";
        //const LaserGroup& lg = calibrationValuesMode3[laserGroupIdx];
        //for (std::map<int, float>::const_iterator it=lg.calibrationValuesLaser.begin(); it!=lg.calibrationValuesLaser.end(); ++it) {
          //int laser = it->first;
          //float calibrationValue = it->second;
          //std::cout << PVARC(laser)<<PVARN(calibrationValue);
        //}
      //}
      calibrationValuesMode3.updateLaserIndices();
      updateMaxUncalibratedIntensitiesForLaserGroups();
    }
    calibrationValuesFile.close();
  }
  else
    return false;
  
  return true;
}

void RemissionCalibrationResult::normalizeCalibrationFactors(float normalization) {
  if (mode==1) { // full graph
    for (std::map<int, std::vector<float> >::iterator it=calibrationValuesMode1.begin(); it!=calibrationValuesMode1.end(); ++it) {
      for (size_t calibrationValueIdx=0;  calibrationValueIdx<it->second.size(); ++calibrationValueIdx)
        it->second[calibrationValueIdx] *= normalization;
    }
  }
  else if (mode==2) { // no full graph
    for (CalibIndependentPerLaser::iterator it=calibrationValuesMode2.begin(); it!=calibrationValuesMode2.end(); ++it) {
      std::vector<float>& calibrationValuesIncidenceAngle = it->second.calibrationValuesIncidenceAngle;
      std::vector<float>& calibrationValuesRange = it->second.calibrationValuesRange;
      float angleNormalizationZero = calibrationValuesIncidenceAngle[0];
      for (size_t calibrationValueIdx=0;  calibrationValueIdx<calibrationValuesIncidenceAngle.size(); ++calibrationValueIdx)
        calibrationValuesIncidenceAngle[calibrationValueIdx] /= angleNormalizationZero;
      for (size_t calibrationValueIdx=0;  calibrationValueIdx<calibrationValuesRange.size(); ++calibrationValueIdx)
        calibrationValuesRange[calibrationValueIdx] *= normalization*angleNormalizationZero;
    }
  }
  else {  // no full graph with laser groups
    for (size_t laserGroupIdx=0; laserGroupIdx<calibrationValuesMode3.size(); ++laserGroupIdx) {
      LaserGroup& lg = calibrationValuesMode3[laserGroupIdx];
      float rangeNormalization = lg.calibrationValuesRange[0];
      for (size_t cellIdx=0; cellIdx<lg.calibrationValuesRange.size(); ++cellIdx)
        lg.calibrationValuesRange[cellIdx] /= rangeNormalization;
      float angleNormalization = lg.calibrationValuesIncidenceAngle[0];
      for (size_t cellIdx=0; cellIdx<lg.calibrationValuesIncidenceAngle.size(); ++cellIdx)
        lg.calibrationValuesIncidenceAngle[cellIdx] /= angleNormalization;
      float intensityNormalization = lg.calibrationValuesIntensity[0];
      for (size_t cellIdx=0; cellIdx<lg.calibrationValuesIntensity.size(); ++cellIdx)
        lg.calibrationValuesIntensity[cellIdx] /= intensityNormalization;
      for (std::map<int, float>::iterator it=lg.calibrationValuesLaser.begin(); it!=lg.calibrationValuesLaser.end(); ++it)
        it->second *= normalization*rangeNormalization*angleNormalization*intensityNormalization;
    }
  }
}

void RemissionCalibrationResult::LaserGroups::updateLaserIndices() {
  laserIndices.clear();
  for (size_t laserGroupIdx=0; laserGroupIdx<size(); ++laserGroupIdx) {
    const LaserGroup& lg = at(laserGroupIdx);
    for (std::map<int, float>::const_iterator it=lg.calibrationValuesLaser.begin(); it!=lg.calibrationValuesLaser.end(); ++it)
      laserIndices[it->first] = laserGroupIdx;
  }
}

std::set<int> RemissionCalibrationResult::getAvailableLaserIds() const {
  std::set<int> ret;
  
  if (mode==1) {
    for (std::map<int, std::vector<float> >::const_iterator it=calibrationValuesMode1.begin(); it!=calibrationValuesMode1.end(); ++it) {
      ret.insert(it->first);
    }
  }
  else if (mode==2) {
    for (CalibIndependentPerLaser::const_iterator it=calibrationValuesMode2.begin(); it!=calibrationValuesMode2.end(); ++it)
      ret.insert(it->first);
  }
  else {
    for (size_t laserGroupIdx=0; laserGroupIdx<calibrationValuesMode3.size(); ++laserGroupIdx) {
      const RemissionCalibrationResult::LaserGroup& lg = calibrationValuesMode3[laserGroupIdx];
      for (std::map<int, float>::const_iterator it=lg.calibrationValuesLaser.begin(); it!=lg.calibrationValuesLaser.end(); ++it) 
        ret.insert(it->first);
    }
  }
  return ret;
}

void RemissionCalibrationResult::printCalibrationValues() const {
  float maxCalibratedRange;
  cellToRange(noOfCellsX-1, maxCalibratedRange);
  
  if (mode==1) {
    //std::ofstream file("calibration.dat");
    for (std::map<int, std::vector<float> >::const_iterator it=calibrationValuesMode1.begin(); it!=calibrationValuesMode1.end(); ++it) {
      int laser = it->first;
      float minRange = minMaxRanges.find(laser)->second.first,
            maxRange = minMaxRanges.find(laser)->second.second;
      const std::vector<float>& calibrationValues = it->second;
      std::cout << "\n-----Calibration values for laser "<<laser<<"-----\n";
      for (int cellY=0; cellY<noOfCellsY; ++cellY) {
        for (int cellX=0; cellX<noOfCellsX; ++cellX) {
          float range;
          cellToRange(cellX, range);
          if (maxCalibratedRange>minRange && (range<minRange || range>maxRange))
            continue;
          float incidenceAngle;
          cellToIncidenceAngle(cellY, incidenceAngle);
          float calibrationValue = calibrationValues[cellY*noOfCellsX+cellX];
          std::cout << PVARC(range)<<PVARA(incidenceAngle)<<": "<<calibrationValue<<"\n";
          //file << laser << " " << range << " " << incidenceAngle<<" " << calibrationValue << "\n";
        }
      }
    }
    //file.close();
  }
  else if (mode==2) {
    for (CalibIndependentPerLaser::const_iterator it=calibrationValuesMode2.begin(); it!=calibrationValuesMode2.end(); ++it) {
      int laser = it->first;
      float minRange = minMaxRanges.find(laser)->second.first,
            maxRange = minMaxRanges.find(laser)->second.second;
      std::cout << "\n-----Calibration values for laser "<<laser<<"-----\n";
      const RemissionCalibrationResult::CalibIndependent& c = it->second;
      for (int cellX=0; cellX<noOfCellsX; ++cellX) {
        float range;
        cellToRange(cellX, range);
        if (maxCalibratedRange>minRange && (range<minRange || range>maxRange))
          continue;
        float calibrationValue = c.calibrationValuesRange[cellX];
        std::cout << PVARC(range)<<PVARN(calibrationValue);
      }
      for (int cellY=0; cellY<noOfCellsY; ++cellY) {
        float incidenceAngle;
        cellToIncidenceAngle(cellY, incidenceAngle);
        float calibrationValue = c.calibrationValuesIncidenceAngle[cellY];
        std::cout << PVARAC(incidenceAngle)<<PVARN(calibrationValue);
      }
      for (int cellZ=0; cellZ<noOfCellsZ; ++cellZ) {
        float intensity;
        cellToIntensity(cellZ, laser, intensity);
        float calibrationValue = c.calibrationValuesIntensity[cellZ];
        std::cout << PVARC(intensity)<<PVARN(calibrationValue);
      }
    }
  }
  else {
    for (size_t laserGroupIdx=0; laserGroupIdx<calibrationValuesMode3.size(); ++laserGroupIdx) {
      const RemissionCalibrationResult::LaserGroup& lg = calibrationValuesMode3[laserGroupIdx];
      float minRange=1e10, maxRange=0;
      for (std::map<int, float>::const_iterator it=lg.calibrationValuesLaser.begin(); it!=lg.calibrationValuesLaser.end(); ++it) {
        minRange = std::min(minRange, minMaxRanges.find(it->first)->second.first);
        maxRange = std::max(maxRange, minMaxRanges.find(it->first)->second.second);
      }
      
      std::cout << "Laser group "<<laserGroupIdx+1<<":\n";
      for (int cellX=0; cellX<noOfCellsX; ++cellX) {
        float range;
        cellToRange(cellX, range);
        if (range<minRange || range>maxRange)
          continue;
        float calibrationValue = lg.calibrationValuesRange[cellX];
        std::cout << PVARC(range)<<PVARN(calibrationValue);
      }
      for (int cellY=0; cellY<noOfCellsY; ++cellY) {
        float incidenceAngle;
        cellToIncidenceAngle(cellY, incidenceAngle);
        float calibrationValue = lg.calibrationValuesIncidenceAngle[cellY];
        std::cout << PVARAC(incidenceAngle)<<PVARN(calibrationValue);
      }
      for (int cellZ=0; cellZ<noOfCellsZ; ++cellZ) {
        float intensity;
        cellToNormalizedIntensity(cellZ, intensity);
        intensity *= lg.maxUncalibratedIntensity;
        float calibrationValue = lg.calibrationValuesIntensity[cellZ];
        std::cout << PVARC(intensity)<<PVARN(calibrationValue);
      }
      for (std::map<int, float>::const_iterator it=lg.calibrationValuesLaser.begin(); it!=lg.calibrationValuesLaser.end(); ++it) {
        int laser = it->first;
        float calibrationValue = it->second;
        std::cout << PVARC(laser)<<PVARN(calibrationValue);
      }
    }
  }
  
  //  GNUPLOT:
  //    set palette model RGB defined (0 "red", 32 "blue", 63 "green")
  //    plot './calibration.dat' using 2:4:1 palette
}

void RemissionCalibrationResult::copyMode3ToMode2() {
  calibrationValuesMode2.clear();
  for (size_t laserGroupIdx=0; laserGroupIdx<calibrationValuesMode3.size(); ++laserGroupIdx) {
    const RemissionCalibrationResult::LaserGroup& lg = calibrationValuesMode3[laserGroupIdx];
    for (std::map<int, float>::const_iterator it=lg.calibrationValuesLaser.begin(); it!=lg.calibrationValuesLaser.end(); ++it) {
      int laser = it->first;
      float calibrationValue = it->second;
      RemissionCalibrationResult::CalibIndependent& c = calibrationValuesMode2[laser];
      for (int cellX=0; cellX<noOfCellsX; ++cellX)
        c.calibrationValuesRange.push_back(calibrationValue*lg.calibrationValuesRange[cellX]);
      c.calibrationValuesIncidenceAngle.insert(c.calibrationValuesIncidenceAngle.begin(),
                                               lg.calibrationValuesIncidenceAngle.begin(),lg.calibrationValuesIncidenceAngle.end());
      c.calibrationValuesIntensity.insert(c.calibrationValuesIntensity.begin(),
                                          lg.calibrationValuesIntensity.begin(),lg.calibrationValuesIntensity.end());
    }
  }
}

void RemissionCalibrationResult::saveGnuplotDataFiles() const {
  float maxCalibratedRange;
  cellToRange(noOfCellsX-1, maxCalibratedRange);
  
  if (mode==1) {
    std::ofstream file("calibration_mode1.dat");
    for (std::map<int, std::vector<float> >::const_iterator it=calibrationValuesMode1.begin(); it!=calibrationValuesMode1.end(); ++it) {
      int laser = it->first;
      float minRange = minMaxRanges.find(laser)->second.first,
            maxRange = minMaxRanges.find(laser)->second.second;
      const std::vector<float>& calibrationValues = it->second;
      for (int cellY=0; cellY<noOfCellsY; ++cellY) {
        for (int cellX=0; cellX<noOfCellsX; ++cellX) {
          float range;
          cellToRange(cellX, range);
          if (maxCalibratedRange>minRange && (range<minRange || range>maxRange))
            continue;
          float incidenceAngle;
          cellToIncidenceAngle(cellY, incidenceAngle);
          float calibrationValue = calibrationValues[cellY*noOfCellsX+cellX];
          file << laser << " " << range << " " << incidenceAngle<<" " << calibrationValue << "\n";
        }
      }
    }
    file.close();
    std::cout << "Wrote calibration values for mode 1. You can use a gnuplot script like this to visualize it:\n"
              << "set palette model RGB defined (0 \"red\", 32 \"blue\", 63 \"green\")\n"
              << "plot './calibration.dat' using 2:4:1 palette\n";
  }
  else if (mode==2) {
    std::ofstream fileRange("calibrationRange_mode2.dat"),
                  fileAngle("calibrationAngle_mode2.dat");
                  //fileIntensity("calibrationIntensity_mode2.dat");
    for (int cellX=0; cellX<noOfCellsX; ++cellX) {
      float range;
      cellToRange(cellX, range);
      fileRange << range;
      for (CalibIndependentPerLaser::const_iterator it=calibrationValuesMode2.begin(); it!=calibrationValuesMode2.end(); ++it) {
        float minRange = minMaxRanges.find(it->first)->second.first,
              maxRange = minMaxRanges.find(it->first)->second.second;
        if (maxCalibratedRange>minRange && (range<minRange || range>maxRange))
          fileRange << " 0";
        else
          fileRange << " "<<it->second.calibrationValuesRange[cellX];
      }
      fileRange << "\n";
    }
    for (int cellY=0; cellY<noOfCellsY; ++cellY) {
      float incidenceAngle;
      cellToIncidenceAngle(cellY, incidenceAngle);
      fileAngle << incidenceAngle;
      for (CalibIndependentPerLaser::const_iterator it=calibrationValuesMode2.begin(); it!=calibrationValuesMode2.end(); ++it)
        fileAngle << " "<<it->second.calibrationValuesIncidenceAngle[cellY];
      fileAngle << "\n";
    }
    //for (int cellZ=0; cellZ<noOfCellsZ; ++cellZ) {
      //float intensity;
      //cellToIntensity(cellZ, laser, intensity);
      //for (CalibIndependentPerLaser::const_iterator it=calibrationValuesMode2.begin(); it!=calibrationValuesMode2.end(); ++it)
        //fileIntensity << " "<<it->second.calibrationValuesIntensity[cellZ];
      //fileIntensity << "\n";
    //}
    fileRange.close();
    fileAngle.close();
    //fileIntensity.close();
  }
  else {
    std::ofstream fileRange("calibrationRange_mode3.dat"),
                  fileAngle("calibrationAngle_mode3.dat"),
                  //fileIntensity("calibrationIntensity_mode3.dat"),
                  fileLaser("calibrationLaser_mode3.dat");
    for (int cellX=0; cellX<noOfCellsX; ++cellX) {
      float range;
      cellToRange(cellX, range);
      fileRange << range;
      for (size_t laserGroupIdx=0; laserGroupIdx<calibrationValuesMode3.size(); ++laserGroupIdx) {
        const RemissionCalibrationResult::LaserGroup& lg = calibrationValuesMode3[laserGroupIdx];
        float minRange=1e10, maxRange=0;
        for (std::map<int, float>::const_iterator it=lg.calibrationValuesLaser.begin(); it!=lg.calibrationValuesLaser.end(); ++it) {
          minRange = std::min(minRange, minMaxRanges.find(it->first)->second.first);
          maxRange = std::max(maxRange, minMaxRanges.find(it->first)->second.second);
        }
        if (maxCalibratedRange>minRange && (range<minRange || range>maxRange))
          fileRange << " 0";
        else
          fileRange << " "<<lg.calibrationValuesRange[cellX];
      }
      fileRange << "\n";
    }
    for (int cellY=0; cellY<noOfCellsY; ++cellY) {
      float incidenceAngle;
      cellToIncidenceAngle(cellY, incidenceAngle);
      fileAngle << incidenceAngle;
      for (size_t laserGroupIdx=0; laserGroupIdx<calibrationValuesMode3.size(); ++laserGroupIdx)
        fileAngle << " "<<calibrationValuesMode3[laserGroupIdx].calibrationValuesIncidenceAngle[cellY];
      fileAngle << "\n";
    }
    //for (int cellZ=0; cellZ<noOfCellsZ; ++cellZ) {
      //float intensity;
      //cellToNormalizedIntensity(cellZ, intensity);
      //intensity *= lg.maxUncalibratedIntensity;
      //for (size_t laserGroupIdx=0; laserGroupIdx<calibrationValuesMode3.size(); ++laserGroupIdx) {
        //fileIntensity << " "<<calibrationValuesMode3[laserGroupIdx].calibrationValuesIntensity[cellZ];
      //fileIntensity << "\n";
    //}
    for (size_t laserGroupIdx=0; laserGroupIdx<calibrationValuesMode3.size(); ++laserGroupIdx) {
      const RemissionCalibrationResult::LaserGroup& lg = calibrationValuesMode3[laserGroupIdx];
      for (std::map<int, float>::const_iterator it=lg.calibrationValuesLaser.begin(); it!=lg.calibrationValuesLaser.end(); ++it) {
        int laser = it->first;
        float calibrationValue = it->second;
        fileLaser << laserGroupIdx<<" "<<laser<<" "<<calibrationValue<<"\n";
      }
    }
    fileRange.close();
    fileAngle.close();
    //fileIntensity.close();
    fileLaser.close();
  }
}

void RemissionCalibrationResult::setMaxUncalibratedIntensities(const std::map<int,float>& newMaxUncalibratedIntensities) {
  maxUncalibratedIntensities = newMaxUncalibratedIntensities;
  updateMaxUncalibratedIntensitiesForLaserGroups();
}

void RemissionCalibrationResult::updateMaxUncalibratedIntensitiesForLaserGroups() {
  //std::cout << __PRETTY_FUNCTION__<<" called.\n";
  if (mode!=3)
    return;
  for (size_t laserGroupIdx=0; laserGroupIdx<calibrationValuesMode3.size(); ++laserGroupIdx) {
    RemissionCalibrationResult::LaserGroup& lg = calibrationValuesMode3[laserGroupIdx];
    lg.maxUncalibratedIntensity = 0.0f;
    for (std::map<int, float>::const_iterator it=lg.calibrationValuesLaser.begin(); it!=lg.calibrationValuesLaser.end(); ++it)
      lg.maxUncalibratedIntensity = std::max(lg.maxUncalibratedIntensity, maxUncalibratedIntensities[it->first]);
    //std::cout << "Maximum uncalibrated intensity of laser group "<<laserGroupIdx<<" is "<<lg.maxUncalibratedIntensity<<".\n";
  }
}
