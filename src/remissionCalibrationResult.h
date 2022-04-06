// This software is licenced under the Creative Commons
// (Attribution-NonCommercial-ShareAlike):
// http://creativecommons.org/licenses/by-nc-sa/4.0

// Author: Bastian Steder <steder@informatik.uni-freiburg.de>


#ifndef REMISSION_CALIBRATION_RESULT_H
#define REMISSION_CALIBRATION_RESULT_H

#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <cmath>
//#include "aislib/stuff/macros.h"

class RemissionCalibrationResult {
  public:
    //-----STRUCTS & TYPEDEFS-----
    
    //-----CONSTRUCTOR & DESTRUCTOR-----
    RemissionCalibrationResult();
    
    //-----PUBLIC STATIC FUNCTIONS-----

    //-----PUBLIC FUNCTIONS-----
    void setNoOfCells(int newNoOfCellsX=50, int newNoOfCellsY=1, int newNoOfCellsZ=1);
    void saveCalibrationValues(const std::string& filename) const;
    void normalizeCalibrationFactors(float normalization);
    bool readCalibrationValuesFromFile(const std::string& fileName);

    void setMaxUncalibratedIntensities(const std::map<int,float>& newMaxUncalibratedIntensities);
    void updateMaxUncalibratedIntensitiesForLaserGroups();
    
    inline float getCalibratedIntensity(float range, float incidenceAngle, float intensity, int laser) const;
    
    inline float getCalibrationFactor(float cellX, float cellY, float cellZ, int laser) const;
    
    //! value>1.0 (range_0=0.0, range_(i+1)=rangeStepFactor*(range_i-1.0))
    void setRangeStepFactor(float newRangeStepFactor=1.1f) { rangeStepFactor = newRangeStepFactor; }
    
    inline void rangeToCell(float range, float& cell) const;
    inline void rangeToCell(float range, int&   cell) const;
    inline void cellToRange(float cell, float& range) const;
    
    inline void incidenceAngleToCell(float incidenceAngle, float& cell) const;
    inline void incidenceAngleToCell(float incidenceAngle, int&   cell) const;
    inline void cellToIncidenceAngle(float cell, float& incidenceAngle) const;
    
    inline void normalizedIntensityToCell(float intensity, float& cell) const;
    inline void normalizedIntensityToCell(float intensity, int& cell) const;
    inline void cellToNormalizedIntensity(float cell, float& intensity) const;
    inline void intensityToCell(float intensity, int laser, float& cell) const;
    inline void intensityToCell(float intensity, int laser, int&  cell) const;
    inline void cellToIntensity(float cell, int laser, float& intensity) const;
    
    inline void rangeAndIncidenceAngleToCell(float range, float incidenceAngle, float& cellX, float& cellY) const;
    inline void rangeAndIncidenceAngleToCell(float range, float incidenceAngle, int&   cellX, int&   cellY) const;
    
    std::set<int> getAvailableLaserIds() const;

    void printCalibrationValues() const;
    void saveGnuplotDataFiles() const;
    void copyMode3ToMode2();
    
    //-----PUBLIC MEMBERS-----
    typedef std::map<int, std::vector<float> > CalibPerLaser;
    
    struct CalibIndependent {
      std::vector<float> calibrationValuesRange,
                         calibrationValuesIncidenceAngle,
                         calibrationValuesIntensity;
    };
    typedef std::map<int, CalibIndependent> CalibIndependentPerLaser;
    
    struct LaserGroup {
      std::vector<float> calibrationValuesRange;
      std::vector<float> calibrationValuesIncidenceAngle;
      std::vector<float> calibrationValuesIntensity;
      std::map<int, float> calibrationValuesLaser;
      float maxUncalibratedIntensity;
    };
    struct LaserGroups : public std::vector<LaserGroup> {
      std::map<int,int> laserIndices;
      void updateLaserIndices();
    };
    
    CalibPerLaser calibrationValuesMode1;
    //calibrationValuesRangePerLaser, calibrationValuesIncidenceAnglePerLaser;
    CalibIndependentPerLaser calibrationValuesMode2;
    LaserGroups calibrationValuesMode3;
    int mode;  // Mode 1: Treat influence of lasers, angles, and ranges individually (full graph)
               // Mode 2: Treat lasers individually and influence of angles, ranges, and intensities as independent
               // Mode 3: Treat lasers in groups and angles and ranges as independent.
    int noOfCellsX, noOfCellsY, noOfCellsZ;  // x range, y angle, z intensity
    float incidenceAngleToCellFactor, cellToIncidenceAngleFactor;
    float rangeStepFactor;
    
    std::map<int,std::pair<float,float> > minMaxRanges;
    std::map<int,float> maxUncalibratedIntensities;

  protected:    
    void updateDependingValues();
};

#include "remissionCalibrationResult.hpp"

#endif
