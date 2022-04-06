// This software is licenced under the Creative Commons
// (Attribution-NonCommercial-ShareAlike):
// http://creativecommons.org/licenses/by-nc-sa/4.0

// Author: Bastian Steder <steder@informatik.uni-freiburg.de>


float RemissionCalibrationResult::getCalibrationFactor(float cellX, float cellY, float cellZ, int laser) const {
  float normalizationFactor = 1.0f;
  if (mode==1) {
    float leftX = floorf(cellX),
          topY  = floorf(cellY),
          d1 = cellX-leftX,
          d2 = cellY-topY;
    
      //std::cout << PVARC(cellX)<<PVARC(cellY);
    
    float topLeft, topRight, bottomLeft, bottomRight;
    
    std::map<int, std::vector<float> >::const_iterator it = calibrationValuesMode1.find(laser);
    const std::vector<float>& calibrationValues = it->second;
    size_t topLeftIdx = lrintf(topY*noOfCellsX+leftX);
    topLeft     = calibrationValues[topLeftIdx];
    topRight    = (leftX+1>=noOfCellsX ? topLeft : calibrationValues[topLeftIdx+1]);
    bottomLeft  = (topY+1>=noOfCellsY  ? topLeft : calibrationValues[topLeftIdx+noOfCellsX]);
    bottomRight = (leftX+1>=noOfCellsX ? bottomLeft :
                    (topY+1>=noOfCellsY  ? topRight : calibrationValues[topLeftIdx+noOfCellsX+1]));
    normalizationFactor = (1.0f-d1)*(1.0f-d2) * topLeft +
                          d1*(1.0f-d2)        * topRight +
                          (1.0f-d1)*d2        * bottomLeft +
                          d1*d2               * bottomRight;
  }
  else if (mode==2) { // no full graph
    const CalibIndependent& c = calibrationValuesMode2.find(laser)->second;
    int lowerCellX = cellX,
        lowerCellY = cellY,
        lowerCellZ = cellZ;
    float normalizationFactorX, normalizationFactorY, normalizationFactorZ;
    if (lowerCellX >= noOfCellsX-1)
      normalizationFactorX = c.calibrationValuesRange[noOfCellsX-1];
    else {
      float factor = cellX-lowerCellX;
      normalizationFactorX = (1.0f-factor)*c.calibrationValuesRange[lowerCellX] + factor*c.calibrationValuesRange[lowerCellX+1];
    }
    if (lowerCellY >= noOfCellsY-1)
      normalizationFactorY = c.calibrationValuesIncidenceAngle[noOfCellsY-1];
    else {
      float factor = cellY-lowerCellY;
      normalizationFactorY = (1.0f-factor)*c.calibrationValuesIncidenceAngle[lowerCellY] + factor*c.calibrationValuesIncidenceAngle[lowerCellY+1];
    }
    if (lowerCellZ >= noOfCellsZ-1)
      normalizationFactorZ = c.calibrationValuesIntensity[noOfCellsZ-1];
    else {
      float factor = cellZ-lowerCellZ;
      normalizationFactorZ = (1.0f-factor)*c.calibrationValuesIntensity[lowerCellZ] + factor*c.calibrationValuesIntensity[lowerCellZ+1];
    }
    normalizationFactor = normalizationFactorX*normalizationFactorY*normalizationFactorZ;
  }
  else if (mode==3) {  // no full graph with laser groups
    int laserGroupIdx = calibrationValuesMode3.laserIndices.find(laser)->second;
    const LaserGroup& lg = calibrationValuesMode3[laserGroupIdx];
    int lowerCellX = cellX,
        lowerCellY = cellY,
        lowerCellZ = cellZ;
    float normalizationFactorX, normalizationFactorY, normalizationFactorZ;
    if (lowerCellX >= noOfCellsX-1)
      normalizationFactorX = lg.calibrationValuesRange[noOfCellsX-1];
    else {
      float factor = cellX-lowerCellX;
      normalizationFactorX = (1.0f-factor)*lg.calibrationValuesRange[lowerCellX] + factor*lg.calibrationValuesRange[lowerCellX+1];
    }
    if (lowerCellY >= noOfCellsY-1)
      normalizationFactorY = lg.calibrationValuesIncidenceAngle[noOfCellsY-1];
    else {
      float factor = cellY-lowerCellY;
      normalizationFactorY = (1.0f-factor)*lg.calibrationValuesIncidenceAngle[lowerCellY] + factor*lg.calibrationValuesIncidenceAngle[lowerCellY+1];
    }
    if (lowerCellZ >= noOfCellsZ-1)
      normalizationFactorZ = lg.calibrationValuesIntensity[noOfCellsZ-1];
    else {
      float factor = cellZ-lowerCellZ;
      normalizationFactorZ = (1.0f-factor)*lg.calibrationValuesIntensity[lowerCellZ] + factor*lg.calibrationValuesIntensity[lowerCellZ+1];
    }
    normalizationFactor = lg.calibrationValuesLaser.find(laser)->second*normalizationFactorX*normalizationFactorY*normalizationFactorZ;
  }
  
  return normalizationFactor;
}

float RemissionCalibrationResult::getCalibratedIntensity(float range, float incidenceAngle, float intensity, int laser) const {
  float cellX, cellY, cellZ;
  rangeAndIncidenceAngleToCell(range, incidenceAngle, cellX, cellY);
  intensityToCell(intensity, laser, cellZ);
  float normalizationFactor = getCalibrationFactor(cellX, cellY, cellZ, laser);
  float normalizedIntensity = normalizationFactor*intensity;
  //std::cout << PVARC(cellX)<<PVARC(cellY)<<PVARC(range)<<PVARC(incidenceAngle)<<PVARC(laser)<<PVARC(intensity)<<PVARC(normalizationFactor)<<PVARC(normalizedIntensity);
  return normalizedIntensity;
}

void RemissionCalibrationResult::rangeToCell(float range, float& cell) const {
  float logBasis = std::log(rangeStepFactor);  // TODO use member variable
  cell = std::min(float(noOfCellsX-1), std::log(range+1.0f) / logBasis);
}
void RemissionCalibrationResult::rangeToCell(float range, int& cell) const {
  float logBasis = std::log(rangeStepFactor);  // TODO use member variable
  cell = std::min(noOfCellsX-1, int(lrintf(std::log(range+1.0f) / logBasis)));
}
void RemissionCalibrationResult::cellToRange(float cell, float& range) const {
  range = std::pow(rangeStepFactor, cell)-1.0f;
}

void RemissionCalibrationResult::incidenceAngleToCell(float incidenceAngle, float& cell) const {
  cell = incidenceAngleToCellFactor*incidenceAngle;
}
void RemissionCalibrationResult::incidenceAngleToCell(float incidenceAngle, int& cell) const {
  cell = lrintf(incidenceAngleToCellFactor*incidenceAngle);
}
void RemissionCalibrationResult::cellToIncidenceAngle(float cell, float& incidenceAngle) const {
  incidenceAngle = cell*cellToIncidenceAngleFactor;
}

void RemissionCalibrationResult::normalizedIntensityToCell(float intensity, float& cell) const {
  cell = std::min(float(noOfCellsZ-1), (noOfCellsZ-1)*intensity);
}
void RemissionCalibrationResult::normalizedIntensityToCell(float intensity, int& cell) const {
  float cellf;
  normalizedIntensityToCell(intensity, cellf);
  cell = lrintf(cellf);
}
void RemissionCalibrationResult::intensityToCell(float intensity, int laser, float& cell) const {
  float normalizedIntensity;
  if (mode==3) {
    const LaserGroup& lg = calibrationValuesMode3[calibrationValuesMode3.laserIndices.find(laser)->second];
    normalizedIntensity = intensity/lg.maxUncalibratedIntensity;
  }
  else
    normalizedIntensity = intensity/maxUncalibratedIntensities.find(laser)->second;
  normalizedIntensityToCell(normalizedIntensity, cell);
}
void RemissionCalibrationResult::intensityToCell(float intensity, int laser, int& cell) const {
  float cellf;
  intensityToCell(intensity, laser, cellf);
  cell = lrintf(cellf);
}
void RemissionCalibrationResult::cellToNormalizedIntensity(float cell, float& intensity) const {
  intensity = cell/std::max(1, noOfCellsZ-1);
}
void RemissionCalibrationResult::cellToIntensity(float cell, int laser, float& intensity) const {
  cellToNormalizedIntensity(cell, intensity);
  if (mode==3)
    intensity *= calibrationValuesMode3[calibrationValuesMode3.laserIndices.find(laser)->second].maxUncalibratedIntensity;
  else
    intensity *= maxUncalibratedIntensities.find(laser)->second;
}

void RemissionCalibrationResult::rangeAndIncidenceAngleToCell(float range, float incidenceAngle, float& cellX, float& cellY) const {
  rangeToCell(range, cellX);
  incidenceAngleToCell(incidenceAngle, cellY);
}

// the cell id that this range and incident angle falls into
void RemissionCalibrationResult::rangeAndIncidenceAngleToCell(float range, float incidenceAngle, int& cellX, int& cellY) const {
  rangeToCell(range, cellX);
  incidenceAngleToCell(incidenceAngle, cellY);
}
