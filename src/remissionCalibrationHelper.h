// This software is licenced under the Creative Commons
// (Attribution-NonCommercial-ShareAlike):
// http://creativecommons.org/licenses/by-nc-sa/4.0

// Author: Bastian Steder <steder@informatik.uni-freiburg.de>


#ifndef REMISSION_CALIBRATION_HELPER_H
#define REMISSION_CALIBRATION_HELPER_H

#include <iostream>
#include <vector>
#include <set>
#include <Eigen/Geometry>

#include "g2o/core/sparse_optimizer.h"

#include "remissionCalibrationResult.h"
#include "pointXYZIVPL.h"

// Forward declarations
class QApplication;
namespace pcl {
  template <typename PointType>
  class PointCloud;
}

class RemissionCalibrationHelper {
  public:
    //-----STRUCTS & TYPEDEFS-----
    typedef pcl::PointXYZIVPL PointType;
    typedef pcl::PointCloud<PointType> PointCloudType;
    
    struct Parameters {
      Parameters() : mapCellSize(0.5f), maxIntensityDeviation(0.02), minNoOfMeasurementsForDataPointMean(10),
                     rangeStepFactorForPointAveraging(1.1), noOfIncidenceAngleCellsForAveraging(10),
                     noOfIterationsForG2O(-1), noOfIterationsForEM(-1), useExpectationMaximizationMethod(true),
                     gridWeight(1.0f), measurementWeight(1.0f), useRobustKernel(false), outlierPercentage(1.0f),
                     autoStopRmseChange(1e-4), useMode3ResultToInitializeMode2(false) {}
      
      float mapCellSize;  //!< Cell size in meters, used for clustering points together
      float maxIntensityDeviation;  //!< Maximum allowed deviation for intensity values in a cell (for values normalized to [0,1])
      size_t minNoOfMeasurementsForDataPointMean;  //!< Minimum number of points to average before adding to map
      float rangeStepFactorForPointAveraging;      //!< Range factor to split up averaging
      float noOfIncidenceAngleCellsForAveraging;   //!< no of cells to split up averaging
      int noOfIterationsForG2O;                    //!< No of iterations performed by the optimizer. -1 means automatic.
      int noOfIterationsForEM;                     //!< No of iterations performed by the EM algorithm.
      bool useExpectationMaximizationMethod;       //!< Should the calibration procedure be based on EM?
      float gridWeight;                            //!< Weight of the edges defining the cell grid
      float measurementWeight;                     //!< Weight of the actual measurement edges
      bool useRobustKernel;                        //!< Use robust kernel for the optimization
      float outlierPercentage;                     //!< Percentage of edges with highest error to be ignored
      float autoStopRmseChange;                    //!< Maximum change in RMSE to stop the EM iterations
      bool useMode3ResultToInitializeMode2;        //!< Use result from mode 3 to initialize mode 2 optimization
    };
    
    //#pragma pack(1)
    struct CollectedPoint {
      CollectedPoint() : x(0),y(0),z(0),
                         //nX(0),nY(0),nZ(0),
                         range(0),incidenceAngle(0),intensity(0),laser(0) {} 
      
      void projectOntoCellMiddle(const RemissionCalibrationResult& calibration);
      
      float x, y, z;
      //float nX, nY, nZ;
      float range;
      float incidenceAngle;
      float intensity;
      int laser;
    };
    class CollectedPoints : public std::vector<CollectedPoint> {
      public:
        CollectedPoints() : std::vector<CollectedPoint>(), cellValue(-1), cellBad(false) {}
        float getAverageCalibratedIntensity(const RemissionCalibrationResult& calibration) const;
        pcl::PointXYZI getPointRepresentingCell() const;
        float getChi2(const RemissionCalibrationResult& calibration) const;
        float getRMSE(const RemissionCalibrationResult& calibration) const;
        float cellValue;
        bool cellBad;
    };
    class CollectedPointsMap : public std::vector<CollectedPoints> {
      public:
        std::set<int> getLaserIds() const;
        std::map<int,std::pair<float,float> > getMinMaxRanges() const;
        std::map<int,float> getMaxIntensities() const;
        float getMaxIntensity() const;
        std::map<int,float> getAverageIntensities() const;
        float getAverageIntensity() const;
        float getAverageCalibratedIntensity(const RemissionCalibrationResult& calibration) const;
        float getMaxCalibratedIntensity(const RemissionCalibrationResult& calibration) const;
        void calibrateIntensities(const RemissionCalibrationResult& calibration);
        void normalizeIntensities(const std::map<int,float>& normalizationValues, bool reciprocal=false);
        void normalizeIntensities(float normalizationValue);
        void setCellValuesToCalibratedMeans(const RemissionCalibrationResult& calibration);
        float getChi2(const RemissionCalibrationResult& calibration) const;
        size_t getNoOfDataPoints() const;
        float getRMSE(const RemissionCalibrationResult& calibration) const;
        void projectOntoCellMiddles(const RemissionCalibrationResult& calibration);
        void removeExtremeIntensities(float percentage);
        void removeOutliers(float percentage, const RemissionCalibrationResult& calibration);
        void removeOutlierCells(float percentage, const RemissionCalibrationResult& calibration);
        std::map<int,int> getErrorDistribution(const RemissionCalibrationResult& calibration, float stepSize=0.1) const;
        std::map<int, float> getErrorPercentages(const RemissionCalibrationResult& calibration) const;
        
        void writeMapCloudToDisk() const;
        void writeCalibratedMapCloudToDisk(const RemissionCalibrationResult& calibration) const;
        void writeToDisk(const std::string& filename) const;
        bool readFromDisk(const std::string& filename);
    };
    
    //-----CONSTRUCTOR & DESTRUCTOR-----
    RemissionCalibrationHelper();
    
    //-----PUBLIC STATIC FUNCTIONS-----
    static inline float getCosIncidenceAngle(const PointType& point, const Eigen::Vector3f& normal);
    static void runG2oOptimization(CollectedPointsMap& mapCells, int noOfIterations,
                                   bool performIntensityInfluenceCalibration=false);
    static void runG2oOptimizationMode1(const CollectedPointsMap& mapCells, g2o::SparseOptimizer& optimizer);
    static void runG2oOptimizationMode2(const CollectedPointsMap& mapCells, g2o::SparseOptimizer& optimizer,
                                        bool performIntensityInfluenceCalibration=false);
    static void runG2oOptimizationMode3(const CollectedPointsMap& mapCells, g2o::SparseOptimizer& optimizer,
                                        bool performIntensityInfluenceCalibration=false);
    static void polynomialApproximation();
    //static void showCalibrationValues();
    static void setCalibrationToInitialGuess(const CollectedPointsMap& mapCells);
    static bool calculateCalibrationBasedOnMapCells(const std::string& filename, bool skipInitialization=false);
    static bool calculateMapCells(const std::string& folder, CollectedPointsMap& mapCells, const std::string& mapCellsFilename="");
    static bool writeCalibratedClouds(const std::string& folder, const std::string& calibrationValuesFilename);
    static bool writeFakeDataset(const std::string& folder, float standardDeviation=0.0f);
    static std::vector<std::string> getScanFilenames(const std::string& folder);
    static void initializeCalibrationValues(const std::set<int>& laserIds);
    static void initializeCalibrationValues(const std::map<int,float>& initialValuesPerLaser);
    static void calculateInitialRMSE(float& maxNormalizationRMSE, float& avgNormalizationRMSE);
    static std::vector<std::set<int> > parseLaserGroupString(const std::string& inputString);
    static std::set<int> getSkipLaserSet();

    //-----PUBLIC MEMBERS-----
    static QApplication* qApplication;
    static RemissionCalibrationResult remissionCalibrationResult;
    static Parameters parameters;
    static std::string laserGroupString;
    static std::string skipLasersString;
};

#include "remissionCalibrationHelper.hpp"

#endif
