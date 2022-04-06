// This software is licenced under the Creative Commons
// (Attribution-NonCommercial-ShareAlike):
// http://creativecommons.org/licenses/by-nc-sa/4.0

// Author: Bastian Steder <steder@informatik.uni-freiburg.de>


#include "remissionCalibrationHelper.h"

#include <pcl/console/print.h>
#include <pcl/io/pcd_io.h>

#include <limits>
#include <fstream>
#include <dirent.h>

#define BACKWARD_HAS_DW 1

using namespace std;
//#include "aislib/stuff/string_tools.h"
//#include "aislib/stuff/misc.h"
//#include "aislib/stuff/timeutil.h"
#include "ais3dTools/basics/pclDownwardsCompatibility.h"
#include "ais3dTools/basics/hashMap3D/hash_map_3d.h"
#include "ais3dTools/basics/filesys_tools.h"
#include "ais3dTools/basics/transformationRepresentation.h"

//#include <pcl/common/polynomial_calculations.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/kdtree/impl/kdtree_flann.hpp>
#include <pcl/common/vector_average.h>

#include "g2o/core/optimization_algorithm_factory.h"
#include "g2o/core/robust_kernel_factory.h"
#include "g2o/core/sparse_optimizer_terminate_action.h"
#include "g2o/core/optimization_algorithm_levenberg.h"
G2O_USE_OPTIMIZATION_LIBRARY(cholmod);
G2O_USE_OPTIMIZATION_LIBRARY(csparse);
G2O_USE_OPTIMIZATION_LIBRARY(dense);

// tony added
#include <g2o/core/optimization_algorithm_gauss_newton.h>
#include <g2o/core/block_solver.h>
#include <g2o/solvers/dense/linear_solver_dense.h>

#include "g2oGraphElements.h"

#include <omp.h> 

RemissionCalibrationHelper::Parameters RemissionCalibrationHelper::parameters;
std::string RemissionCalibrationHelper::laserGroupString;
std::string RemissionCalibrationHelper::skipLasersString;
RemissionCalibrationResult RemissionCalibrationHelper::remissionCalibrationResult;
RemissionCalibrationResult& remissionCalibrationResult = RemissionCalibrationHelper::remissionCalibrationResult;

std::string g2oAlgorithmName = "gn_var";
// std::string g2oAlgorithmName = "lm_var";
//std::string g2oAlgorithmName = "gn_dense";
//std::string g2oAlgorithmName = "lm_dense";

double huberWidth = 1.0;

float initialValue = 1.0;
//const double stabilizingEdgeWeightFactor = 1.0;

static int& mode = remissionCalibrationResult.mode;

RemissionCalibrationResult::CalibPerLaser& calibrationValuesMode1 = remissionCalibrationResult.calibrationValuesMode1;
RemissionCalibrationResult::CalibIndependentPerLaser& calibrationValuesMode2 = remissionCalibrationResult.calibrationValuesMode2;
RemissionCalibrationResult::LaserGroups& calibrationValuesMode3 = remissionCalibrationResult.calibrationValuesMode3;

VerticesMode1 verticesMode1;
VerticesMode2 verticesMode2;
VerticesMode3 verticesMode3;
std::vector<RemissionLearnerVertex*> verticesForMapCells;

int& noOfCellsX = remissionCalibrationResult.noOfCellsX,
   & noOfCellsY = remissionCalibrationResult.noOfCellsY,
   & noOfCellsZ = remissionCalibrationResult.noOfCellsZ;
//int& polynomialDegree = remissionCalibrationResult.polynomialDegree;
//pcl::BivariatePolynomialT<double>& polynomial = remissionCalibrationResult.polynomial;

#include "ais3dTools/basics/sparsePointCloud.h"
typedef Ais3dTools::SparsePointCloudT<RemissionCalibrationHelper::CollectedPoints> PointDatabase;


struct PointDyn
{
  float dynProb;
};
POINT_CLOUD_REGISTER_POINT_STRUCT (PointDyn,
    (float, dynProb, dynProb)
)
typedef pcl::PointCloud<PointDyn> PointCloudDyn;


G2oPostIterationAction g2oPostIterationAction;  // Object to update GUI during g2o-iterations and copy over optimized values

static double getLocalTime() {
  struct timeval ts; 
  gettimeofday(&ts,0);
  return ts.tv_sec + ts.tv_usec*1e-6;
}

RemissionCalibrationHelper::RemissionCalibrationHelper() {
}

typedef RemissionCalibrationHelper::PointType PointType;
typedef std::vector<PointType, Eigen::aligned_allocator<PointType> > PointListType;
typedef RemissionCalibrationHelper::PointCloudType PointCloudType;
typedef pcl::KdTreeFLANN<PointType> KdTreeType;
typedef std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f> > VectorOfEigenVector3f;

static const Eigen::Vector3f invalidNormal(0,0,0);

void estimateResolution(const PointCloudType& cloud, float& resX, float& resY) {
  //std::cout << PVARC(cloud.width)<<PVARN(cloud.height);
  float angleSumX=0, angleSumY=0;
  int noOfValuesX=0, noOfValuesY=0;
  for (size_t y=1; y<cloud.height; ++y) {
    for (size_t x=1; x<cloud.width; ++x) {
      const PointType& p = cloud[y*cloud.width+x];
      if (!std::isfinite(p.x))
        continue;
      Eigen::Vector3f sensorPos(p.vp_x, p.vp_y, p.vp_z);
      Eigen::Vector3f direction1 = (p.getVector3fMap()-sensorPos).normalized();
      const PointType& neighborX = cloud[y*cloud.width+(x-1)];
      const PointType& neighborY = cloud[(y-1)*cloud.width+x];
      if (std::isfinite(neighborX.x)) {
        Eigen::Vector3f direction2 = (neighborX.getVector3fMap()-sensorPos).normalized();
        float angle = acosf(direction1.dot(direction2));
        if (std::isfinite(angle)) {
          angleSumX += angle;
          ++noOfValuesX;
        }
        else {
          //std::cerr << PVARC(p)<<PVARC(neighborX)<<PVARN(direction1.dot(direction2));
        }
      }
      if (std::isfinite(neighborY.x)) {
        Eigen::Vector3f direction2 = (neighborY.getVector3fMap()-sensorPos).normalized();
        float angle = acosf(direction1.dot(direction2));
        if (std::isfinite(angle)) {
          angleSumY += angle;
          ++noOfValuesY;
        }
        else {
          //std::cerr << PVARC(p)<<PVARC(neighborY)<<PVARN(direction1.dot(direction2));
        }
      }
    }
  }
  resX = angleSumX/noOfValuesX;
  resY = angleSumY/noOfValuesY;
  //std::cout << PVARAC(resX)<<PVARAN(resY);
}

// 生成每个离散的cell中最多只有一个点的稀疏点云
template <typename CloudType>
void createSparsePointCloud(const CloudType& cloud, CloudType& sparseCloud, float resolution) {
  sparseCloud.clear();
  Ais3dTools::HashMap3D<bool> hashMap3D(resolution);
  for (size_t pointIdx=0; pointIdx<cloud.size(); ++pointIdx) {
    const typename CloudType::PointType& p = cloud[pointIdx];
    if (!std::isfinite(p.x))
      continue;
    Eigen::Vector3i cellIndex = hashMap3D.worldToGrid(p.getVector3fMap());
    bool& cellValue = hashMap3D.findOrCreate(cellIndex, false);
    // 如果该cell中已经有元素，就不再向其中增加点
    if (cellValue) {
      continue;
    }
    sparseCloud.push_back(p);
    // 每个cell里存了一个bool，用于表示该cell里有没有元素
    cellValue = true;
  }
}

// 通过PCA计算一个点的surface normal，使用kdTree中存储的数据得到附近的点
bool getNormal(const PointType& point, const KdTreeType& kdTree, Eigen::Vector3f& normal, float maxDistance, float* quality=NULL) {
  normal = invalidNormal;
  
  if (quality != NULL)
    *quality = 0.0f;
  
  const PointCloudType& cloud = *const_cast<KdTreeType*>(&kdTree)->getInputCloud();  // HACK for old PCL versions 
  
  if (!std::isfinite(point.x))
    return false;

  Eigen::Vector3f sensorPos(point.vp_x, point.vp_y, point.vp_z);
  Eigen::Vector3f viewingDirection = (sensorPos-point.getVector3fMap()).normalized();
  
  std::vector<int> neighborIndices;
  std::vector<float> neighborSquaredDistances;
  //neighborIndices.clear();
  //neighborSquaredDistances.clear();
  kdTree.radiusSearch(point, maxDistance, neighborIndices, neighborSquaredDistances);
  pcl::VectorAverage3f vectorAverage;
  for (size_t neighborIndicesIdx=0; neighborIndicesIdx<neighborIndices.size(); ++neighborIndicesIdx) {
    size_t neighborIdx = neighborIndices[neighborIndicesIdx];
    const PointType& neighbor = cloud[neighborIdx];
    vectorAverage.add(neighbor.getVector3fMap());
  }
  
  if (vectorAverage.getNoOfSamples()<3)
    return false;
  
  //vectorAverage.getEigenVector1(normal);
  Eigen::Vector3f eigenValues, eigenVector2, eigenVector3;
  vectorAverage.doPCA(eigenValues, normal, eigenVector2, eigenVector3);
  
  if (quality != NULL)
    *quality = (1.0f - eigenValues[0]/eigenValues[1]) * (eigenValues[1]/eigenValues[2]);   // Quality regarding planarity
  
  // 确保normal和入射光线在同一个象限内
  if (normal.dot(viewingDirection) < 0.0f)
    normal *= -1.0f;
  
  //std::cout << "("<<normal[0]<<","<<normal[1]<<","<<normal[2]<<")["<<vectorAverage.getNoOfSamples()<<"] "<<std::flush;
  
  return true;
}

bool getNormal(const PointCloudType& cloud, size_t index, Eigen::Vector3f& normal, bool newCloud) {
  const float minStepAngle = DEG2RAD(1.0f);
  normal = invalidNormal;
  
  static int stepSizeX=1, stepSizeY=1;
  if (newCloud) {
    float resX, resY;
    estimateResolution(cloud, resX, resY);
    stepSizeX = std::max(1, int(lrintf(minStepAngle/resX)));
    stepSizeY = std::max(1, int(lrintf(minStepAngle/resY)));
    //std::cout << "Cloud has rectangular structure with approximated resolution of "<<RAD2DEG(resX)<<"deg x "<<RAD2DEG(resY)<<"deg.\n";
    //std::cout << PVARC(stepSizeX)<<PVARN(stepSizeY);
  }
  
  if (cloud.height <= 1) {
    std::cerr << __PRETTY_FUNCTION__<<": Error, point cloud has to be organized (width and height > 1).\n";
    return false;
  }
  
  //normal = Eigen::Vector3f(1.0f, 0.0f, 0.0f);  return true;  // HACK!
  
  const PointType point = cloud[index];
  
  if (!std::isfinite(point.x))
    return false;
  
  Eigen::Vector3f p1 = point.getVector3fMap();
  Eigen::Vector3f sensorPos(point.vp_x, point.vp_y, point.vp_z);
  int width = cloud.width,
      height=cloud.height;
  int y=index/width,
      x=index-y*width;
  
  Eigen::Vector3f viewingDirection = (sensorPos-p1).normalized();
  
  pcl::VectorAverage3f vectorAverage;
  for (int y2=std::max(0, y-stepSizeY); y2<=std::min(y+stepSizeY,height-1); ++y2) {
    for (int x2=std::max(0, x-stepSizeX); x2<=std::min(x+stepSizeX,width-1); ++x2) {
      //std::cout << PVARC(x)<<PVARC(y)<<PVARC(x2)<<PVARN(y2);
      const PointType& p2 = cloud.points[y2*width + x2];
      if (!std::isfinite(p2.x))
        continue;
      vectorAverage.add(p2.getVector3fMap());
    }
  }
  //std::cout << PVARN(vectorAverage.getNoOfSamples());
  if (vectorAverage.getNoOfSamples()<3)
    return false;
  vectorAverage.getEigenVector1(normal);
  
  //int triangles[18] = {-1,-1, 0,-1, 1,1, 1,0, 1,1, 0,1, -1,1, -1,0, -1,-1};
  
  ////std::vector<float> cosIncidenceAngles;
  //for (int triangleIdx=0; triangleIdx<8; ++triangleIdx) {
    //int x2=x+stepSizeX*triangles[2*triangleIdx],   y2=y+stepSizeY*triangles[2*triangleIdx+1],
        //x3=x+stepSizeX*triangles[2*triangleIdx+2], y3=y+stepSizeY*triangles[2*triangleIdx+3];
    //if (x2<0 || x2>=width || y2<0|| y2>=height || x3<0 || x3>=width || y3<0 || y3>=height)
      //continue;
    //Eigen::Vector3f p2 = cloud.points[y2*width + x2].getVector3fMap(),
                    //p3 = cloud.points[y3*width + x3].getVector3fMap();
    //if (!isfinite(p2[0]) || !isfinite(p3[0]))
      //continue;
    //Eigen::Vector3f alternativeNormal = (p2-p1).cross(p3-p1).normalized();
    //if (std::abs(alternativeNormal.dot(normal))<0.8)
      //return false;
    ////float cosIncidenceAngle = std::abs(viewingDirection.dot(normal));
    ////cosIncidenceAngles.push_back(cosIncidenceAngle);
  //}
  ////if (cosIncidenceAngles.empty())
    ////return 0.0f;
  ////std::sort(cosIncidenceAngles.begin(), cosIncidenceAngles.end());
  ////cosIncidenceAngleToUse = cosIncidenceAngles.back();
  
  if (normal.dot(viewingDirection) < 0.0f)
    normal *= -1.0f;
  
  return true;
}


void getNormals(const boost::shared_ptr<PointCloudType> cloud, VectorOfEigenVector3f& normals,
                boost::shared_ptr<KdTreeType> kdTree=boost::shared_ptr<KdTreeType>(), float maxDistance=0.1f)
{
  normals.clear();
  normals.resize(cloud->size());
  
  if (cloud->height > 1) {  // Organized cloud
    for (size_t pointIdx=0; pointIdx<cloud->size(); ++pointIdx)
      getNormal(*cloud, pointIdx, normals[pointIdx], pointIdx==0);
  }
  else {  // Unorganized cloud
    std::cout << "Cloud is unordered (does not have rectangular structure) - Will use kdTree to extract normals.\n";
    if (kdTree == NULL) {
      kdTree.reset(new KdTreeType);
      kdTree->setInputCloud(cloud);
    }
    for (size_t pointIdx=0; pointIdx<cloud->size(); ++pointIdx)
      getNormal(cloud->points[pointIdx], *kdTree, normals[pointIdx], maxDistance);
  
  }
}

struct NormalWithQuality {
  NormalWithQuality() : quality(0.0f) {}
  Eigen::Vector3f normal;
  float quality;
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

typedef Ais3dTools::HashMap3D<std::shared_ptr<NormalWithQuality> > NormalsMap;

// 对点云将采样，并计算其各cell中点的surface normal
template <typename CloudType>
void createSparseNormalsMap(const CloudType& cloud, NormalsMap& normalsMap) {
  normalsMap.clear();
  float resolution = normalsMap.getGridSize();
  pcl::KdTreeFLANN<PointType> kdTree;
  boost::shared_ptr<PointCloudType> sparseCloudPtr(new PointCloudType);
  PointCloudType& sparseCloud = *sparseCloudPtr;
  createSparsePointCloud(cloud, *sparseCloudPtr, resolution);
  kdTree.setInputCloud(sparseCloudPtr);
  
  for (size_t pointIdx=0; pointIdx<sparseCloud.size(); ++pointIdx) {
    const typename CloudType::PointType& p = sparseCloud[pointIdx];
    std::shared_ptr<NormalWithQuality>& normalWithQuality = normalsMap.cellFromPoint(p.getVector3fMap());
    normalWithQuality.reset(new NormalWithQuality);
    getNormal(p, kdTree, normalWithQuality->normal, 3.0f*resolution, &normalWithQuality->quality);
  }
}

//void addStabilizingEdgeForAllVertices(g2o::SparseOptimizer& optimizer) {
  //double stabilizingEdgeWeight = stabilizingEdgeWeightFactor*optimizer.edges().size();
  //for (g2o::SparseOptimizer::VertexIDMap::iterator it=optimizer.vertices().begin(); it!=optimizer.vertices().end(); ++it) { 
    //g2o::HyperGraph::Vertex* vertex = it->second;
    //RemissionLearnerUnaryEdge* stabilizingEdge = new RemissionLearnerUnaryEdge;
    //stabilizingEdge->information()(0,0) = stabilizingEdgeWeight;
    //stabilizingEdge->vertices()[0] = vertex;
    //optimizer.addEdge(stabilizingEdge);
  //}
//}

//void addGlobalStabilizingEdge(g2o::SparseOptimizer& optimizer) {
  //double stabilizingEdgeWeight = stabilizingEdgeWeightFactor*optimizer.edges().size();
  //RemissionLearnerMultiEdge_n* stabilizingEdge = new RemissionLearnerMultiEdge_n(optimizer.vertices().size());
  //stabilizingEdge->information()(0,0) = stabilizingEdgeWeight;
  //size_t counter = 0;
  //for (g2o::SparseOptimizer::VertexIDMap::iterator it=optimizer.vertices().begin(); it!=optimizer.vertices().end(); ++it) { 
    //stabilizingEdge->vertices()[counter++] = it->second;
  //}
  //optimizer.addEdge(stabilizingEdge);
//}

void RemissionCalibrationHelper::runG2oOptimization(CollectedPointsMap& mapCells, int noOfIterations,
                                                    bool performIntensityInfluenceCalibration)
{
  g2o::SparseOptimizer optimizer;
  //optimizer.setVerbose(true);

  // g2o::OptimizationAlgorithmProperty solverProperty;
  // g2o::OptimizationAlgorithmFactory* solverFactory = g2o::OptimizationAlgorithmFactory::instance();
  // g2o::OptimizationAlgorithm* g2oAlgorithm = solverFactory->construct(g2oAlgorithmName, solverProperty);
  // optimizer.setAlgorithm(g2oAlgorithm);


  // tony added
  // since solver factory is not working, select solver directly instead
  typedef g2o::BlockSolver<g2o::BlockSolverTraits<1, 1>> BlockSolverType;
  typedef g2o::LinearSolverDense<BlockSolverType::PoseMatrixType> LinearSolverType;

  auto solver = new g2o::OptimizationAlgorithmGaussNewton(
      g2o::make_unique<BlockSolverType>(g2o::make_unique<LinearSolverType>()));
  optimizer.setAlgorithm(solver);

  //if (g2oAlgorithmName=="lm_var" || g2oAlgorithmName=="lm_dense")
    //static_cast<g2o::OptimizationAlgorithmLevenberg*>(g2oAlgorithm)->setUserLambdaInit(1e5);
  
  g2o::SparseOptimizerTerminateAction* terminateAction = NULL;
  if (noOfIterations < 0) {
    //std::cout << "# setup termination criterion based on the gain of the iteration" << endl;
    noOfIterations = 1e3; //std::numeric_limits<int>::max();
    terminateAction = new g2o::SparseOptimizerTerminateAction;
    optimizer.addPostIterationAction(terminateAction);
  }
  
  if (mode==1) { // full graph
    runG2oOptimizationMode1(mapCells, optimizer);
  }
  else if (mode==2) { // no full graph
    runG2oOptimizationMode2(mapCells, optimizer, performIntensityInfluenceCalibration);
  }
  else {  // no full graph with laser groups
    runG2oOptimizationMode3(mapCells, optimizer, performIntensityInfluenceCalibration);
  }
  
  g2oPostIterationAction.mapCells = &mapCells;
  g2oPostIterationAction.verticesMode1 = &verticesMode1;
  g2oPostIterationAction.verticesMode2 = &verticesMode2;
  g2oPostIterationAction.verticesMode3 = &verticesMode3;
  g2oPostIterationAction.verticesForMapCells = &verticesForMapCells;
  
  g2oPostIterationAction.remissionCalibrationResult = &remissionCalibrationResult;
  optimizer.addPostIterationAction(&g2oPostIterationAction);
  optimizer.initializeOptimization();
  
  if (parameters.useRobustKernel) {
    g2o::AbstractRobustKernelCreator* creator = g2o::RobustKernelFactory::instance()->creator("Cauchy");
    //g2o::AbstractRobustKernelCreator* creator = g2o::RobustKernelFactory::instance()->creator("Saturated");
    if (!creator) {
      std::cerr << "Not a a valid robust kernel" << endl;
    }
    else {
      for (g2o::SparseOptimizer::EdgeSet::const_iterator it = optimizer.edges().begin(); it != optimizer.edges().end(); ++it) {
        g2o::OptimizableGraph::Edge* e = static_cast<g2o::OptimizableGraph::Edge*>(*it);
        e->setRobustKernel(creator->construct());
        e->robustKernel()->setDelta(huberWidth);
      }
    }
  }
  
  //optimizer.computeActiveErrors();
  //double chi2 = optimizer.chi2();
  //std::cout << "G2O Chi2 " << FIXED(chi2) << " (RMSE "<<sqrt(chi2/optimizer.edges().size())<<")\n";
  
  optimizer.optimize(noOfIterations);
  
  //optimizer.computeActiveErrors();
  //chi2 = optimizer.chi2();
  //std::cout << "G2O Chi2 " << FIXED(chi2) << " (RMSE "<<sqrt(chi2/optimizer.edges().size())<<")\n";
  
  //remissionCalibrationResult.printCalibrationValues();
}

void RemissionCalibrationHelper::runG2oOptimizationMode1(const CollectedPointsMap& mapCells, g2o::SparseOptimizer& optimizer) {
  // Add vertices and grid structure
  int vertexId = 0;
  verticesMode1.clear();
  for (std::map<int, std::vector<float> >::iterator it=calibrationValuesMode1.begin(); it!=calibrationValuesMode1.end(); ++it) {
    int laser = it->first;
    //std::cout << PVARN(laser);
    std::vector<float>& calibrationValues = it->second;
    std::vector<RemissionLearnerVertex*>& vertices = verticesMode1[laser];
    
    for (int cellY=0; cellY<noOfCellsY; ++cellY) {
      for (int cellX=0; cellX<noOfCellsX; ++cellX) {
        int index = cellY*noOfCellsX + cellX;
        RemissionLearnerVertex* vertex = new RemissionLearnerVertex;
        vertex->setEstimate(calibrationValues[index]);
        vertex->setId(vertexId++);
        vertices.push_back(vertex);
        optimizer.addVertex(vertex);
        if (cellX>0) {
          optimizer.addEdge(new RemissionCalibrationGridEdge(vertices[cellY*noOfCellsX+cellX-1], vertex, parameters.gridWeight));
        }
        if (cellX==0 && cellY>0) {
          optimizer.addEdge(new RemissionCalibrationGridEdge(vertices[(cellY-1)*noOfCellsX+cellX], vertex, parameters.gridWeight));
        }
      }
    }
  }
  
  // Add constraints from actual measurements
  verticesForMapCells.clear();
  for (size_t mapCellIdx=0; mapCellIdx<mapCells.size(); ++mapCellIdx) {
    const CollectedPoints& mapCell = mapCells[mapCellIdx];
    RemissionLearnerVertex* mapCellVertex = new RemissionLearnerVertex;
    mapCellVertex->setEstimate(mapCell.cellValue);
    mapCellVertex->setId(vertexId++);
    verticesForMapCells.push_back(mapCellVertex);
    optimizer.addVertex(mapCellVertex);
    mapCellVertex->setFixed(true);
    for (size_t dataPointIdx=0; dataPointIdx<mapCell.size(); ++dataPointIdx) {
      const CollectedPoint& cp = mapCell[dataPointIdx];
      int cellX, cellY;
      remissionCalibrationResult.rangeAndIncidenceAngleToCell(cp.range, cp.incidenceAngle, cellX, cellY);
      RemissionLearnerVertex* calibrationValueVertex = verticesMode1[cp.laser][cellY*noOfCellsX+cellX];
      
      RemissionCalibrationEdge_2* edge = new RemissionCalibrationEdge_2;
      edge->setMeasurement(cp.intensity);
      edge->information()(0,0) = parameters.measurementWeight;
      edge->vertices()[0] = calibrationValueVertex;
      edge->vertices()[1] = mapCellVertex;
      optimizer.addEdge(edge);
    }
  }
}

void RemissionCalibrationHelper::runG2oOptimizationMode2(const CollectedPointsMap& mapCells, g2o::SparseOptimizer& optimizer,
                                                         bool performIntensityInfluenceCalibration)
{
  // Add vertices and grid structure
  int vertexId = 0;
  verticesMode2.clear();
  for (RemissionCalibrationResult::CalibIndependentPerLaser::iterator it=calibrationValuesMode2.begin(); it!=calibrationValuesMode2.end(); ++it) {
    int laser = it->first;
    RemissionCalibrationResult::CalibIndependent& c = calibrationValuesMode2[laser];
    CalibIndependentVertices& v = verticesMode2[laser];
    // traverse all range for one laser
    for (int cellX=0; cellX<noOfCellsX; ++cellX) {
      RemissionLearnerVertex* vertex = new RemissionLearnerVertex;
      vertex->setEstimate(c.calibrationValuesRange[cellX]);
      vertex->setId(vertexId++);
      v.verticesRange.push_back(vertex);
      optimizer.addVertex(vertex);
      if (cellX > 0) {
        optimizer.addEdge(new RemissionCalibrationGridEdge(v.verticesRange[cellX-1], vertex, parameters.gridWeight));
      }
      if (performIntensityInfluenceCalibration)
        vertex->setFixed(true);
    }
    
    // traverse all incident angle for one laser
    for (int cellY=0; cellY<noOfCellsY; ++cellY) {
      RemissionLearnerVertex* vertex = new RemissionLearnerVertex;
      vertex->setEstimate(c.calibrationValuesIncidenceAngle[cellY]);
      vertex->setId(vertexId++);
      v.verticesIncidenceAngle.push_back(vertex);
      optimizer.addVertex(vertex);
      if (cellY > 0) {
        optimizer.addEdge(new RemissionCalibrationGridEdge(v.verticesIncidenceAngle[cellY-1], vertex, parameters.gridWeight));
      }
      if (performIntensityInfluenceCalibration)
        vertex->setFixed(true);
    }
    v.verticesIncidenceAngle[0]->setFixed(true);  // First incidence Angle normalization factor should always be 1.0
    
    for (int cellZ=0; cellZ<noOfCellsZ; ++cellZ) {
      RemissionLearnerVertex* vertex = new RemissionLearnerVertex;
      vertex->setEstimate(c.calibrationValuesIntensity[cellZ]);
      vertex->setId(vertexId++);
      v.verticesIntensity.push_back(vertex);
      optimizer.addVertex(vertex);
      if (cellZ > 0) {
        optimizer.addEdge(new RemissionCalibrationGridEdge(v.verticesIntensity[cellZ-1], vertex, parameters.gridWeight));
      }
      if (!performIntensityInfluenceCalibration)
        vertex->setFixed(true);
    }
    v.verticesIntensity[0]->setFixed(true);  // First intensity normalization factor should always be 1.0
  }
  
  // Add constraints from actual measurements
  verticesForMapCells.clear();
  for (size_t mapCellIdx=0; mapCellIdx<mapCells.size(); ++mapCellIdx) {
    const CollectedPoints& mapCell = mapCells[mapCellIdx];
    RemissionLearnerVertex* mapCellVertex = new RemissionLearnerVertex;
    mapCellVertex->setEstimate(mapCell.cellValue);
    mapCellVertex->setId(vertexId++);
    verticesForMapCells.push_back(mapCellVertex);
    optimizer.addVertex(mapCellVertex);
    mapCellVertex->setFixed(true);
    for (size_t dataPointIdx=0; dataPointIdx<mapCell.size(); ++dataPointIdx) {
      const CollectedPoint& cp = mapCell[dataPointIdx];
      CalibIndependentVertices& v = verticesMode2[cp.laser];
      int cellX, cellY, cellZ;
      remissionCalibrationResult.rangeAndIncidenceAngleToCell(cp.range, cp.incidenceAngle, cellX, cellY);
      remissionCalibrationResult.intensityToCell(cp.intensity, cp.laser, cellZ);
      RemissionLearnerVertex* rangeVertex     = v.verticesRange[cellX],
                            * angleVertex     = v.verticesIncidenceAngle[cellY],
                            * intensityVertex = v.verticesIntensity[cellZ];
      RemissionCalibrationEdgeMode2* edge = new RemissionCalibrationEdgeMode2;
      edge->setMeasurement(cp.intensity);
      edge->information()(0,0) = parameters.measurementWeight;
      edge->vertices()[0] = rangeVertex;
      edge->vertices()[1] = angleVertex;
      edge->vertices()[2] = intensityVertex;
      edge->vertices()[3] = mapCellVertex;
      optimizer.addEdge(edge);
    }
  }
}

void RemissionCalibrationHelper::runG2oOptimizationMode3(const CollectedPointsMap& mapCells, g2o::SparseOptimizer& optimizer,
                                                         bool performIntensityInfluenceCalibration)
{
  // Add vertices and grid structure
  int vertexId = 0;
  verticesMode3.clear();
  verticesMode3.resize(calibrationValuesMode3.size());
  for (size_t laserGroupIdx=0; laserGroupIdx<calibrationValuesMode3.size(); ++laserGroupIdx) {
    RemissionCalibrationResult::LaserGroup& lg = calibrationValuesMode3[laserGroupIdx];
    LaserGroupVertices& lgv = verticesMode3[laserGroupIdx];
    
    for (int cellX=0; cellX<noOfCellsX; ++cellX) {
      RemissionLearnerVertex* vertex = new RemissionLearnerVertex;
      vertex->setEstimate(lg.calibrationValuesRange[cellX]);
      vertex->setId(vertexId++);
      lgv.verticesRange.push_back(vertex);
      optimizer.addVertex(vertex);
      if (cellX > 0) {
        optimizer.addEdge(new RemissionCalibrationGridEdge(lgv.verticesRange[cellX-1], vertex, parameters.gridWeight));
      }
      if (performIntensityInfluenceCalibration)
        vertex->setFixed(true);
    }
    //lgv.verticesRange[0]->setFixed(true);  // First incidence angle normalization factor should always be 1.0
    
    for (int cellY=0; cellY<noOfCellsY; ++cellY) {
      RemissionLearnerVertex* vertex = new RemissionLearnerVertex;
      vertex->setEstimate(lg.calibrationValuesIncidenceAngle[cellY]);
      vertex->setId(vertexId++);
      lgv.verticesIncidenceAngle.push_back(vertex);
      optimizer.addVertex(vertex);
      if (cellY > 0) {
        optimizer.addEdge(new RemissionCalibrationGridEdge(lgv.verticesIncidenceAngle[cellY-1], vertex, parameters.gridWeight));
      }
      if (performIntensityInfluenceCalibration)
        vertex->setFixed(true);
    }
    //lgv.verticesIncidenceAngle[0]->setFixed(true);  // First incidence angle normalization factor should always be 1.0
    
    for (int cellZ=0; cellZ<noOfCellsZ; ++cellZ) {
      RemissionLearnerVertex* vertex = new RemissionLearnerVertex;
      vertex->setEstimate(lg.calibrationValuesIntensity[cellZ]);
      vertex->setId(vertexId++);
      lgv.verticesIntensity.push_back(vertex);
      optimizer.addVertex(vertex);
      if (cellZ > 0) {
        optimizer.addEdge(new RemissionCalibrationGridEdge(lgv.verticesIntensity[cellZ-1], vertex, parameters.gridWeight));
      }
      if (!performIntensityInfluenceCalibration)
        vertex->setFixed(true);
    }
    //lgv.verticesIntensity[0]->setFixed(true);  // First intensity normalization factor should always be 1.0
    
    for (std::map<int, float>::const_iterator it=lg.calibrationValuesLaser.begin(); it!=lg.calibrationValuesLaser.end(); ++it) {
      RemissionLearnerVertex* vertex = new RemissionLearnerVertex;
      vertex->setEstimate(it->second);
      vertex->setId(vertexId++);
      lgv.verticesLaser[it->first] = vertex;
      optimizer.addVertex(vertex);
      if (performIntensityInfluenceCalibration)
        vertex->setFixed(true);
    }
  }
  
  // Add constraints from actual measurements
  verticesForMapCells.clear();
  for (size_t mapCellIdx=0; mapCellIdx<mapCells.size(); ++mapCellIdx) {
    const CollectedPoints& mapCell = mapCells[mapCellIdx];
    RemissionLearnerVertex* mapCellVertex = new RemissionLearnerVertex;
    mapCellVertex->setEstimate(mapCell.cellValue);
    mapCellVertex->setId(vertexId++);
    verticesForMapCells.push_back(mapCellVertex);
    optimizer.addVertex(mapCellVertex);
    mapCellVertex->setFixed(true);
    for (size_t dataPointIdx=0; dataPointIdx<mapCell.size(); ++dataPointIdx) {
      const CollectedPoint& cp = mapCell[dataPointIdx];
      int cellX, cellY, cellZ;
      remissionCalibrationResult.rangeAndIncidenceAngleToCell(cp.range, cp.incidenceAngle, cellX, cellY);
      remissionCalibrationResult.intensityToCell(cp.intensity, cp.laser, cellZ);
      int laserGroupIdx = calibrationValuesMode3.laserIndices[cp.laser];
      //std::cout << laserGroupIdx<<", cellZ="<<cellZ<<".\n";
      LaserGroupVertices& lgv = verticesMode3[laserGroupIdx];
      RemissionLearnerVertex* calibrationValueRangeVertex     = lgv.verticesRange[cellX],
                            * calibrationValueAngleVertex     = lgv.verticesIncidenceAngle[cellY],
                            * calibrationValueIntensityVertex = lgv.verticesIntensity[cellZ],
                            * calibrationValueLaserVertex     = lgv.verticesLaser[cp.laser];
      RemissionCalibrationEdgeMode3* edge = new RemissionCalibrationEdgeMode3;
      edge->setMeasurement(cp.intensity);
      edge->information()(0,0) = parameters.measurementWeight;
      edge->vertices()[0] = calibrationValueRangeVertex;
      edge->vertices()[1] = calibrationValueAngleVertex;
      edge->vertices()[2] = calibrationValueIntensityVertex;
      edge->vertices()[3] = calibrationValueLaserVertex;
      edge->vertices()[4] = mapCellVertex;
      optimizer.addEdge(edge);
    }
  }
}


void RemissionCalibrationHelper::polynomialApproximation() {
  //std::vector<Eigen::Vector3d> samplePoints;
  //for (size_t calibrationValueIdx=0;  calibrationValueIdx<calibrationValues.size(); ++calibrationValueIdx) {
    //double range = 1.0f/(std::max(0.4999,double(calibrationValueIdx%noOfCellsX))*cellToRangeReciprocalFactor),
           //incidenceAngle = double(calibrationValueIdx/noOfCellsX)*cellToIncidenceAngleFactor;
    //Eigen::Vector3d sample(std::pow(range,2), cosf(incidenceAngle), calibrationValues[calibrationValueIdx]);
    //samplePoints.push_back(sample);
  //}
  //pcl::PolynomialCalculationsT<double> polynomialCalculations;
  //if (!polynomialCalculations.bivariatePolynomialApproximation(samplePoints, polynomialDegree, polynomial))
    //std::cout << "Could not approximate polynomial.\n";
  //else
    //std::cout << "Resulting polynomial is "<<polynomial<<"\n";
  
  //calibrationValuesPolynomial.clear();
  //for (size_t calibrationValueIdx=0;  calibrationValueIdx<calibrationValues.size(); ++calibrationValueIdx) {
    //double range = 1.0f/(std::max(0.4999,double(calibrationValueIdx%noOfCellsX))*cellToRangeReciprocalFactor),
           //incidenceAngle = double(calibrationValueIdx/noOfCellsX)*cellToIncidenceAngleFactor;
    //calibrationValuesPolynomial.push_back(polynomial.getValue(std::pow(range,2), cosf(incidenceAngle)));
  //}
}

//void RemissionCalibrationHelper::showCalibrationValues() {
  //if (qApplication==NULL || calibrationValuesMode1.empty())
    //return;
  
  //static std::map<int, Ais3dTools::ImageWidget*> hitsWidgets, calibrationValuesWidgets;
  //for (std::map<int, std::vector<float> >::const_iterator it=calibrationValuesMode1.begin(); it!=calibrationValuesMode1.end(); ++it) {
    //int laser = it->first;
    //const std::vector<float>& calibrationValues = it->second;
    //const std::vector<float>& hits = hitsPerLaser[laser];
    //if (hitsWidgets.find(laser) == hitsWidgets.end()) {
      //hitsWidgets[laser] = new Ais3dTools::ImageWidget;
      //calibrationValuesWidgets[laser] = new Ais3dTools::ImageWidget;
    //}
    //Ais3dTools::ImageWidget& hitsWidget = *hitsWidgets[laser];
    //Ais3dTools::ImageWidget& calibrationValuesWidget = *calibrationValuesWidgets[laser];
    //std::stringstream windowTitle;
    //windowTitle << "Hits for laser "<<laser<<"\n";
    //hitsWidget.setWindowTitle(windowTitle.str().c_str());
    //hitsWidget.useFastScaling = true;
    //hitsWidget.setRealImage(&hits[0], noOfCellsX, noOfCellsY);
    
    //windowTitle.str("");
    //windowTitle << "Calibration values for laser "<<laser<<"\n";
    //calibrationValuesWidget.setWindowTitle(windowTitle.str().c_str());
    //calibrationValuesWidget.useFastScaling = true;
    //calibrationValuesWidget.showRealImagesInGrayScales = true;
    //calibrationValuesWidget.setRealImage(&calibrationValues[0], noOfCellsX, noOfCellsY);

    ////static Ais3dTools::ImageWidget* calibrationValuesPolynomialWidget = new Ais3dTools::ImageWidget;
    ////calibrationValuesPolynomialWidget->useFastScaling = true;
    ////calibrationValuesPolynomialWidget->setWindowTitle("calibration values polynomial");
    ////calibrationValuesPolynomialWidget->setRealImage(&calibrationValuesPolynomial[0], noOfCellsX, noOfCellsY, true, 0.0f, 20.0f);
  //}
  
  //while (true) {
    //qApplication->processEvents();
    //usleep(10000);
  //}
//}


//void RemissionCalibrationHelper::processDataPoints(const std::string& dataFilename, const std::string& calibrationValuesFilename) {
  //bool useGui = qApplication!=NULL;
  
  //Ais3dTools::ALUGLViewer* viewerPtr = NULL;
  //Vis::RemissionCalibrationResult visRemissionCalibrationResult(&remissionCalibrationResult); 
  //if (useGui) {
    //setlocale(LC_NUMERIC, "C");
    //viewerPtr = new Ais3dTools::ALUGLViewer;
    //viewerPtr->show();
    //viewerPtr->updateGL();
    //viewerPtr->add(&visRemissionCalibrationResult);
    //viewerPtr->setBackgroundColor(QColor(200, 200, 255));
    //viewerPtr->updateGL();
    //qApplication->processEvents();
    //g2oPostIterationAction.qApplication = qApplication;
    //g2oPostIterationAction.viewer = viewerPtr;
  //}
  
  //dataPoints.clear();
  //if (!dataFilename.empty()) {
    //bool useFakePoints = false;
    ////useFakePoints = true;  // HACK
    //if (!useFakePoints) {
      //std::set<int> skipLaser;
      //std::vector<std::set<int> > parsedLaserGroups = parseLaserGroupString(skipLasersString);
      //for (size_t laserGroupIdx=0; laserGroupIdx<parsedLaserGroups.size(); ++laserGroupIdx) {
        //const std::set<int>& lasers = parsedLaserGroups[laserGroupIdx];
        //for (std::set<int>::const_iterator laserIt=lasers.begin(); laserIt!=lasers.end(); ++laserIt)
          //skipLaser.insert(*laserIt);
      //}
      
      //std::ifstream file(dataFilename.c_str());
      //while (file && !file.eof()) {
        //RemissionCalibrationResult::DataPoint dataPoint;
        //file.read((char*)&dataPoint, sizeof(dataPoint));
        //if (file.eof() || dataPoint.intensity1<1e-4 || dataPoint.intensity2<1e-4)
          //continue;

        //if (skipLaser.find(dataPoint.laser1)!=skipLaser.end() || skipLaser.find(dataPoint.laser2)!=skipLaser.end())
          //continue;
        
        //dataPoints.push_back(dataPoint);
        ////dataPoints.back().incidenceAngle1 = dataPoints.back().intensity1;  // HACK!!!
        ////dataPoints.back().incidenceAngle2 = dataPoints.back().intensity2;  // HACK!!!
      //}
      //file.close();
      ////std::cout << PVARN(dataPoints.size());
    //}
    //else {
      ////for (float range1=1.0; range1<=100; range1+=10.0) {
        ////for (float angle1=0.0; angle1<DEG2RAD(90.0f); angle1+=DEG2RAD(10.0)) {
          ////for (float range2=1.0; range2<=100; range2+=10.0) {
            ////for (float angle2=0.0; angle2<DEG2RAD(90.0f); angle2+=DEG2RAD(10.0)) {
              ////RemissionCalibrationResult::DataPoint dataPoint;
              ////dataPoint.laser1 = dataPoint.laser2 = 0;
              ////dataPoint.range1 = range1;
              ////dataPoint.range2 = range2;
              ////dataPoint.incidenceAngle1 = angle1;
              ////dataPoint.incidenceAngle2 = angle2;
              ////dataPoint.intensity1 = 1.0/dataPoint.range1;
              ////dataPoint.intensity1 /= 1.0 + dataPoint.incidenceAngle2;
              ////dataPoint.intensity2 = 1.0/dataPoint.range2;
              ////dataPoint.intensity2 /= 1.0 + dataPoint.incidenceAngle2;
              ////dataPoints.push_back(dataPoint);
            ////}
          ////}
        ////}
      ////}
      //float start = 0.1, step=0.1, end=10.0;
      //for (float range1=start; range1<=end; range1+=step) {
        //for (float range2=start; range2<=end; range2+=step) {
          //RemissionCalibrationResult::DataPoint dataPoint;
          //dataPoint.laser1 = dataPoint.laser2 = 0;
          //dataPoint.range1 = range1;
          //dataPoint.range2 = range2;
          //dataPoint.incidenceAngle1 = 0;
          //dataPoint.incidenceAngle2 = 0;
          //dataPoint.intensity1 = 1.0/dataPoint.range1;
          //dataPoint.intensity2 = 1.0/dataPoint.range2;
          //dataPoints.push_back(dataPoint);
        //}
      //}

      ////RemissionCalibrationResult::DataPoint dataPoint;
      ////dataPoint.laser1 = dataPoint.laser2 = 0;
      ////dataPoint.range1 = dataPoint.range2 = 0.1;
      ////dataPoint.incidenceAngle1 = dataPoint.incidenceAngle2 = 0.0;
      ////dataPoint.intensity1 = dataPoint.intensity2 = 1.0;
      ////dataPoints.push_back(dataPoint);
    //}
    //std::cout << PVARN(dataPoints.size());

    ////std::map<int, std::pair<float,float> > minMaxRanges = getMinMaxRanges(dataPoints);
    
    //float maxNormalizationRMSE, avgNormalizationRMSE;
    //calculateInitialRMSE(maxNormalizationRMSE, avgNormalizationRMSE);
    //std::cout << "Initial RMSE (max normalization): "<<maxNormalizationRMSE<<", Initial RMSE (avg normalization): "<<avgNormalizationRMSE<<"\n";
    
    //std::set<int> laserIds = remissionCalibrationResult.getLaserIds();
    
    //std::map<int, float> averageIntensities;
    //for (std::set<int>::const_iterator it=laserIds.begin(); it!=laserIds.end(); ++it)
      //averageIntensities[*it] = getAverageIntensity(*it, dataPoints);
    
    //if (!calibrationValuesFilename.empty()) {
      //remissionCalibrationResult.readCalibrationValuesFromFile(calibrationValuesFilename);
      ////remissionCalibrationResult.printCalibrationValues();
      //std::cout << PVARN(getRMSE());
      //if (mode==1) { // full graph
        //for (std::map<int, std::vector<float> >::iterator it=calibrationValuesMode1.begin(); it!=calibrationValuesMode1.end(); ++it)
          //for (size_t calibrationValueIdx=0;  calibrationValueIdx<it->second.size(); ++calibrationValueIdx)
            //it->second[calibrationValueIdx] *= averageIntensities[it->first];
      //}
      //else if (mode==2) { // no full graph
        //for (std::map<int, std::vector<float> >::iterator it=calibrationValuesRangePerLaser.begin(); it!=calibrationValuesRangePerLaser.end(); ++it)
          //for (size_t calibrationValueIdx=0;  calibrationValueIdx<it->second.size(); ++calibrationValueIdx)
            //it->second[calibrationValueIdx] *= averageIntensities[it->first];
      //}
      //else {  // no full graph with laser groups
        //for (size_t laserGroupIdx=0; laserGroupIdx<calibrationValuesMode3.size(); ++laserGroupIdx) {
          //RemissionCalibrationResult::LaserGroup& lg = calibrationValuesMode3[laserGroupIdx];
          //for (std::map<int, float>::iterator it=lg.calibrationValuesLaser.begin(); it!=lg.calibrationValuesLaser.end(); ++it)
            //it->second *= averageIntensities[it->first];
        //}
      //}
    //}
    //else {
      //initializeCalibrationValues(remissionCalibrationResult.getLaserIds());
    //}

    
    //for (size_t dataPointIdx=0; dataPointIdx<dataPoints.size(); ++dataPointIdx) {
      //RemissionCalibrationResult::DataPoint& p = dataPoints[dataPointIdx];
      //p.intensity1 /= averageIntensities[p.laser1];
      //p.intensity2 /= averageIntensities[p.laser2];
    //}
    
    ////std::cout << "Initial RMSE: "<<getRMSE()<<".\n";
    //if (useGui) {
      //for (int guiUpdateCounter=0; guiUpdateCounter<200; ++guiUpdateCounter) {
        //qApplication->processEvents();
        //viewerPtr->updateGL();
        //usleep(10000);
      //}
    //}
    
    ////float minRange, maxRange;
    ////getMinMaxRange(dataPoints, minRange, maxRange);
    ////std::cout << "Setting minimum range to "<<minRange<<"m and maximum range to "<<maxRange<<"m.\n";
    
    //int noOfRuns = 3;
    //noOfRuns = 1;
    //for (int run=0; run<noOfRuns && parameters.noOfIterationsForG2O!=0; ++run) {
      //if (run==0) {
        //runG2oOptimization(&dataPoints, NULL, parameters.noOfIterationsForG2O);
      //}
      //else {
        //// Step 1: Project onto cell middle
        //std::cout << "Projecting points onto cell middles.\n";
        //RemissionCalibrationResult::DataPoints learningPoints1;
        //for (size_t dataPointIdx=0; dataPointIdx<dataPoints.size(); ++dataPointIdx) {
          //const RemissionCalibrationResult::DataPoint& dataPoint = dataPoints[dataPointIdx];
          //RemissionCalibrationResult::DataPoint newDataPoint = remissionCalibrationResult.projectOntoCellMiddle(dataPoint);
          //learningPoints1.push_back(newDataPoint);
        //}
        
        //// Step 2: Remove 1% of the edges with the highest errors
        //bool removeOutliers = run>=2;
        ////removeOutliers = false;
        
        //RemissionCalibrationResult::DataPoints learningPoints2;
        //RemissionCalibrationResult::DataPoints& learningPoints = (removeOutliers ? learningPoints2 : learningPoints1);
        
        //if (removeOutliers) {
          //std::cout << "Removing outliers.\n";
          //std::vector<float> errors;
          //for (size_t dataPointIdx=0; dataPointIdx<learningPoints1.size(); ++dataPointIdx) {
            //const RemissionCalibrationResult::DataPoint& dataPoint = learningPoints1[dataPointIdx];
            //float value1 = remissionCalibrationResult.getCalibratedIntensity(dataPoint.range1, dataPoint.incidenceAngle1,
                                                                             //dataPoint.intensity1, dataPoint.laser1),
                  //value2 = remissionCalibrationResult.getCalibratedIntensity(dataPoint.range2, dataPoint.incidenceAngle2,
                                                                             //dataPoint.intensity2, dataPoint.laser2);
            //float error = RemissionLearnerEdge::computeError(value1, value2);
            //errors.push_back(error);
          //}
          //std::sort(errors.begin(), errors.end());
          //float keepFactor = 1.0f - 0.01f*(run-1);
          //float maxError = errors[errors.size()*keepFactor];
          //for (size_t dataPointIdx=0; dataPointIdx<learningPoints1.size(); ++dataPointIdx) {
            //const RemissionCalibrationResult::DataPoint& dataPoint = learningPoints1[dataPointIdx];
            //float value1 = remissionCalibrationResult.getCalibratedIntensity(dataPoint.range1, dataPoint.incidenceAngle1,
                                                                             //dataPoint.intensity1, dataPoint.laser1),
                  //value2 = remissionCalibrationResult.getCalibratedIntensity(dataPoint.range2, dataPoint.incidenceAngle2,
                                                                             //dataPoint.intensity2, dataPoint.laser2);
            //float error = RemissionLearnerEdge::computeError(value1, value2);
            //if (error < maxError)
              //learningPoints2.push_back(dataPoint);
          //}
        //}
        //runG2oOptimization(&learningPoints, NULL, parameters.noOfIterationsForG2O);
      //}
    //}
    
    //// Reverse intensity normalization
    //for (size_t dataPointIdx=0; dataPointIdx<dataPoints.size(); ++dataPointIdx) {
      //RemissionCalibrationResult::DataPoint& p = dataPoints[dataPointIdx];
      //p.intensity1 *= averageIntensities[p.laser1];
      //p.intensity2 *= averageIntensities[p.laser2];
    //}
    
    //if (mode==1) { // full graph
      //for (std::map<int, std::vector<float> >::iterator it=calibrationValuesMode1.begin(); it!=calibrationValuesMode1.end(); ++it)
        //for (size_t calibrationValueIdx=0;  calibrationValueIdx<it->second.size(); ++calibrationValueIdx)
          //it->second[calibrationValueIdx] /= averageIntensities[it->first];
    //}
    //else if (mode==2) { // no full graph
      //for (std::map<int, std::vector<float> >::iterator it=calibrationValuesRangePerLaser.begin(); it!=calibrationValuesRangePerLaser.end(); ++it)
        //for (size_t calibrationValueIdx=0;  calibrationValueIdx<it->second.size(); ++calibrationValueIdx)
          //it->second[calibrationValueIdx] /= averageIntensities[it->first];
    //}
    //else {  // no full graph with laser groups
      //for (size_t laserGroupIdx=0; laserGroupIdx<calibrationValuesMode3.size(); ++laserGroupIdx) {
        //RemissionCalibrationResult::LaserGroup& lg = calibrationValuesMode3[laserGroupIdx];
        //for (std::map<int, float>::iterator it=lg.calibrationValuesLaser.begin(); it!=lg.calibrationValuesLaser.end(); ++it)
          //it->second /= averageIntensities[it->first];
      //}
    //}
    ////RemissionCalibrationResult bla = remissionCalibrationResult;
    //remissionCalibrationResult.normalizeCalibrationFactors();
    
    ////for (size_t dataPointIdx=0; dataPointIdx<dataPoints.size(); dataPointIdx+=1) {
      ////const RemissionCalibrationResult::DataPoint& dataPoint = dataPoints[dataPointIdx];
      ////double value1 = remissionCalibrationResult.getCalibratedIntensity(dataPoint.range1, dataPoint.incidenceAngle1,
                                                                        ////dataPoint.intensity1, dataPoint.laser1);
      ////double value2 = bla.getCalibratedIntensity(dataPoint.range1, dataPoint.incidenceAngle1, dataPoint.intensity1, dataPoint.laser1);
      ////std::cout << value1/value2<<" "<<std::flush;
    ////}
    
    //remissionCalibrationResult.printCalibrationValues();
    //std::cout << "Root mean squared error is "<<getRMSE()<<"\n";
    
    //remissionCalibrationResult.saveCalibrationValues("calibrationValues.txt");
    
    ////remissionCalibrationResult.printCalibrationValues();
    
    ////polynomialApproximation();
    ////double squaredError = getMSE(),
           ////squaredErrorPolynomial = getMSEPolynomial();
    ////std::cout << PVARC(squaredError)<<PVARC(squaredErrorPolynomial);
    ////printCalibratedPoints();
    ////calculateHits();
    ////showCalibrationValues();
    
    ////std::vector<SamplePointT<float> > samplePoints;
    ////for (size_t dataPointIdx=0; dataPointIdx<dataPoints.size(); ++dataPointIdx) {
      ////const RemissionCalibrationResult::DataPoint& dataPoint = dataPoints[dataPointIdx];
      ////SamplePointT<float> samplePoint;
      ////samplePoint.x1 = dataPoint.squaredRange1;
      ////samplePoint.y1 = 1.0f/dataPoint.cosIncidenceAngle1;
      ////samplePoint.z1 = dataPoint.intensity1;
      ////samplePoint.x2 = dataPoint.squaredRange2;
      ////samplePoint.y2 = 1.0f/dataPoint.cosIncidenceAngle2;
      ////samplePoint.z2 = dataPoint.intensity2;
      ////samplePoints.push_back(samplePoint);
    ////}
    
    ////BivariatePolynomialT<float> polynomial;
    ////int degree = 3,
        ////fixedParameter = BivariatePolynomialT<float>::getNoOfParametersFromDegree(degree)-(degree+3);
    ////float fixedValue = 0.1;
    ////if (!bivariatePolynomialApproximation(samplePoints, degree, fixedParameter, fixedValue, polynomial)) {
      ////std::cout << "Could not calculate polynomial. :-(\n";
      ////return 1;
    ////}
    
    ////cout << "Resulting polynomial is "<<polynomial<<"\n";
  //}
  ////visRemissionCalibrationResult.showEdges = true;
  //while (useGui) {
    //if (useGui && omp_get_thread_num()==0) {
      //bool updateGui = false;
      //AISLIB_DO_EVERY(0.5, updateGui=true;);
      //if (updateGui) {
        //qApplication->processEvents();
        //viewerPtr->updateGL();
      //}
    //}
    //usleep(1000);
  //}
//}

std::vector<std::string> RemissionCalibrationHelper::getScanFilenames(const std::string& folder) {
  std::vector<std::string> ret;
  
  DIR* directory;
  if((directory  = opendir(folder.c_str())) == NULL) {
    std::cerr << "Could not open directory.\n";
    return ret;
  }
  struct dirent *dirEntry;
  while ((dirEntry = readdir(directory)) != NULL) {
    if (dirEntry->d_type == DT_REG)  // Only regular files
    {
      std::string file_name = dirEntry->d_name;
      // Ignore non-pcd files
      if (file_name.size() <= 4 || file_name.substr(file_name.size()-4, 4)!=".pcd")
        continue;
      // Ignore far ranges files
      if (file_name.size()>15 && file_name.substr(file_name.size()-15, 11)=="_far_ranges")
        continue;
      ret.push_back(dirEntry->d_name);
    }
  }
  closedir(directory);
  std::sort(ret.begin(), ret.end());
  
  return ret;
}

bool readCloud(const std::string& filename, PointCloudType& cloud, bool transformPoints) {
  std::cout << "Now reading \""<<filename<<"\".\n";
  
  //pcl::io::loadPCDFile(pointCloudFilename, cloud);
  pcl::PCLPointCloud2 pointCloudData;
  if (pcl::io::loadPCDFile(filename, pointCloudData, cloud.sensor_origin_, cloud.sensor_orientation_) == -1) {
    cerr << "Could not open file \""<<filename<<"\".\n";
    return false;
  }
  
  fromPCLPointCloud2(pointCloudData, cloud);

  if (pcl::getFieldIndex(pointCloudData, "dynProb") >= 0) {
    std::cout << "Point cloud contains information about dynamic objects.\n";
    PointCloudDyn cloudDyn;
    fromPCLPointCloud2(pointCloudData, cloudDyn);
    for (size_t pointIdx=0; pointIdx<cloud.size(); ++pointIdx) {
      if (cloudDyn[pointIdx].dynProb >= 0.3)
        cloud[pointIdx].x = cloud[pointIdx].y = cloud[pointIdx].z = std::numeric_limits<float>::quiet_NaN();
    }
  }
  else {
    fromPCLPointCloud2(pointCloudData, cloud);
  }
  
  if (pcl::getFieldIndex(pointCloudData, "vp_x") < 0)
  {
    float vp_x = cloud.sensor_origin_.x(),
          vp_y = cloud.sensor_origin_.y(),
          vp_z = cloud.sensor_origin_.z();
    std::cout << "Points do not include viewpoint information.\n"
              << "Will use the whole clouds viewpoint information ("<<vp_x<<","<<vp_y<<","<<vp_z<<").\n";
    for (size_t pointIdx=0; pointIdx<cloud.points.size(); ++pointIdx) {
      PointType& point = cloud.points[pointIdx];
      point.vp_x = vp_x;
      point.vp_y = vp_y;
      point.vp_z = vp_z;
    }
  }
  
  if (pcl::getFieldIndex(pointCloudData, "laser") < 0)
  {
    std::cout << "Points do not include a laser id. Will assume id 0.\n";
    for (size_t pointIdx=0; pointIdx<cloud.points.size(); ++pointIdx) {
      PointType& point = cloud.points[pointIdx];
      point.laser = 0;
    }
  }
  //else {
    //std::cout << "\n";
    //for (size_t pointIdx=0; pointIdx<cloud.points.size(); ++pointIdx)
      //std::cout << cloud.points[pointIdx].laser<<", ";
    //std::cout << "\n";
  //}

  if (!transformPoints)
    return true;
  
  std::string baseFilename = Ais3dTools::getPureFilename(filename);
  std::string infoFilename = baseFilename+"_info.dat";
  std::ifstream infoFile(infoFilename.c_str());
  
  std::stringstream line;
  while (infoFile && !infoFile.eof()) {
    std::string lineString;
    std::getline(infoFile, lineString);
    line.str(lineString);
    std::string token;
    line >> token;
    if (token=="SLAM:")
      break;
    line.str("");
  }
  infoFile.close();
  if (line.str().empty()) {
    std::cout << "Could not find SLAM pose.\n";
    return false;
  }
  float x, y, z, roll, pitch, yaw;
  line >> x >> y >> z >> roll >> pitch >> yaw;
  
  //std::cout << PVARC(x)<<PVARC(y)<<PVARC(z)<<PVARAC(roll)<<PVARAC(pitch)<<PVARAN(yaw);
  
  Eigen::Affine3f robotPose;
  Ais3dTools::TransformationRepresentation::getMatrixFromTranslationAndEuler(
      x, y, z, roll, pitch, yaw, robotPose);
  
  for (size_t pointIdx=0; pointIdx<cloud.size(); ++pointIdx) {
    PointType& p = cloud[pointIdx];
    if (!std::isfinite(p.x))
      continue;
    p.getVector3fMap() = robotPose * p.getVector3fMap();
    Eigen::Vector3f sensorPos(p.vp_x, p.vp_y, p.vp_z);
    sensorPos = robotPose * sensorPos;
    p.vp_x=sensorPos.x(), p.vp_y=sensorPos.y(), p.vp_z=sensorPos.z();
  }
  
  return true;
}

bool RemissionCalibrationHelper::calculateMapCells(
  const std::string& folder, CollectedPointsMap& mapCells, const std::string& mapCellsFilename) {
  //bool useGui = qApplication!=NULL;
  
  std::vector<std::string> scan_file_names = getScanFilenames(folder);

  if (scan_file_names.empty()) {
    std::cerr << "Could not open directory or find any *.pcd files in the folder.\n";
    return false;
  }
  
  //Ais3dTools::ALUGLViewer* viewerPtr = NULL;
  
  //Ais3dTools::PclPointCloudObjectT<PointCloudType> visPointCloud;
  //visPointCloud.setDrawColor(0.0, 0.0, 1.0);
  //pcl::PointCloud<pcl::PointXYZ> mapCloud;
  //Ais3dTools::PclPointCloudObjectT<pcl::PointCloud<pcl::PointXYZ> > visMapCloud(&mapCloud);
  //visMapCloud.setDrawColor(0.0, 1.0, 0.0);
  //double lastGuiUpdate = 0.0;
  
  //if (useGui) {
    //setlocale(LC_NUMERIC, "C");
    //viewerPtr = new Ais3dTools::ALUGLViewer;
    //viewerPtr->show();
    //viewerPtr->updateGL();
    //viewerPtr->add(&visPointCloud);
    //viewerPtr->add(&visMapCloud);
    //viewerPtr->setBackgroundColor(Qt::white);
  //}
  PointDatabase pointDatabase(parameters.mapCellSize);
  
  // 要跳过的一些laser
  std::set<int> skipLaser = getSkipLaserSet();
  
  std::cout << "Step 1: Go through dataset to collect some properties...\n";
  float maxRange = 0.0f;
  std::map<int, float> maxIntensities;
  for (size_t pcdIdx=0; pcdIdx<scan_file_names.size(); ++pcdIdx) {
    std::string pointCloudFilename = folder+"/"+scan_file_names[pcdIdx];
    PointCloudType cloud;
    if (!readCloud(pointCloudFilename, cloud, true))
      continue;
    std::cout << "Read cloud "<<pcdIdx+1<<"/"<<scan_file_names.size()<<" with "<<cloud.size()<<" points.\n";
    //std::cout << "Sensor origin is "<<cloud.sensor_origin_[0]<<","<<cloud.sensor_origin_[1]<<","<<cloud.sensor_origin_[2]<<".\n";
    for (size_t pointIdx=0; pointIdx<cloud.size(); ++pointIdx) {
      const PointType& p = cloud[pointIdx];
      if (!std::isfinite(p.x))
        continue;
      if (skipLaser.find(p.laser)!=skipLaser.end())
        continue;
      if (maxIntensities.find(p.laser)==maxIntensities.end())
        maxIntensities[p.laser] = p.intensity;
      else
        maxIntensities[p.laser] = std::max(maxIntensities[p.laser], p.intensity);
      float range = (p.getVector3fMap()-Eigen::Vector3f(p.vp_x, p.vp_y, p.vp_z)).norm();
      //if (range<0.1 || (p.laser>=100 && range>100))
        //std::cerr << "Weird range measurement ("<<range<<"m for Laser "<<p.laser<<": "
                  //<< "("<<p.x<<","<<p.y<<","<<p.z<<") with viewpoint ("<<p.vp_x<<","<<p.vp_y<<","<<p.vp_z<<")\n";
      maxRange = std::max(maxRange, range);
      //if (useGui && getLocalTime()-lastGuiUpdate>0.1) {
        //lastGuiUpdate = getLocalTime();
        //visPointCloud.setPointCloud(&cloud);
        //qApplication->processEvents();
        //viewerPtr->updateGL();
      //}
    }
  }
  
  for (std::map<int, float>::const_iterator it=maxIntensities.begin(); it!=maxIntensities.end(); ++it)
    std::cout << "Laser "<<it->first<<": "<<it->second<<", ";
  std::cout << "\n";
  
  RemissionCalibrationResult cellHelper;
  cellHelper.setNoOfCells(100, parameters.noOfIncidenceAngleCellsForAveraging);
  cellHelper.setRangeStepFactor(parameters.rangeStepFactorForPointAveraging);
  
  std::cout << "Step 2: Collect averaged cell points...\n";
  for (size_t pcdIdx=0; pcdIdx<scan_file_names.size(); ++pcdIdx) {
    std::string baseFilename = Ais3dTools::getPureFilename(scan_file_names[pcdIdx]);
    std::string pointCloudFilename = folder+"/"+scan_file_names[pcdIdx];
    boost::shared_ptr<PointCloudType> cloudPtr(new PointCloudType);
    PointCloudType& cloud = *cloudPtr;
    if (!readCloud(pointCloudFilename, cloud, true))
      continue;
    std::cout << "Read cloud "<<pcdIdx+1<<"/"<<scan_file_names.size()<<" with "<<cloud.size()<<" points.\n";
    
    //std::cout << "Updating map to extract datapoints for calibration procedure.\n";

    NormalsMap normalsMap(0.2*parameters.mapCellSize);
    createSparseNormalsMap(cloud, normalsMap);
    
    //bool cloudIsUnorganized = cloud.height<=1;
    
    //boost::shared_ptr<pcl::KdTreeFLANN<PointType> > kdTreePtr(new pcl::KdTreeFLANN<PointType>);
    //boost::shared_ptr<PointCloudType> sparseCloudPtr(new PointCloudType);
    //if (cloudIsUnorganized) {
      //float sparseResolution = 0.2 * parameters.mapCellSize;
      //createSparsePointCloud(cloud, *sparseCloudPtr, sparseResolution);
      //kdTreePtr->setInputCloud(sparseCloudPtr);
    //}
    
    std::map<int, std::map<int, std::map<int, PointDatabase> > > clusters; // laserId, rangeCell, incidenceAngleCell, cells of points
    
    //if (useGui)
      //visPointCloud.setPointCloud(&cloud);
    
    for (size_t pointIdx=0; pointIdx<cloud.size(); ++pointIdx) {
      const PointType& p = cloud[pointIdx];
      if (!std::isfinite(p.x))
        continue;
      
      if (p.intensity < 1e-6)
        continue;

      if (skipLaser.find(p.laser)!=skipLaser.end())
        continue;
      
      // 获取预先计算的，每个点对应cell内的surface normal
      NormalWithQuality normalWithQuality = *normalsMap.cellFromPoint(p.getVector3fMap());
      //Eigen::Vector3f normal;
      //if (cloudIsUnorganized)  // Unorganized cloud
        //getNormal(p, *kdTreePtr, normal, 0.5f*parameters.mapCellSize);
      //else  // Organized cloud
        //getNormal(cloud, pointIdx, normal, pointIdx==0);
      
      // TODO: PORQUE hacer esto?
      if (normalWithQuality.normal.squaredNorm()<0.9)  // Invalid normal?
        continue;
      
      //if (normalWithQuality.quality < 0.3)
        //continue;
      
      float cosIncidenceAngle = getCosIncidenceAngle(cloud[pointIdx], normalWithQuality.normal),
            incidenceAngle = acosf(cosIncidenceAngle);
      
      float range = (p.getVector3fMap()-Eigen::Vector3f(p.vp_x, p.vp_y, p.vp_z)).norm();
      int rangeCell;
      cellHelper.rangeToCell(range, rangeCell);

      int incidenceAngleCell;
      cellHelper.incidenceAngleToCell(incidenceAngle, incidenceAngleCell);
      
      std::map<int, PointDatabase>& pointsPerAngle = clusters[p.laser][rangeCell];
      std::map<int, PointDatabase>::iterator pointsIt = pointsPerAngle.find(incidenceAngleCell);
      if (pointsIt == pointsPerAngle.end()) {
        pointsPerAngle[incidenceAngleCell] = PointDatabase(parameters.mapCellSize);
        pointsIt = pointsPerAngle.find(incidenceAngleCell);
      }
      PointDatabase& points = pointsIt->second;
      
      int cellX, cellY, cellZ;
      points.getIndicesFromCoordinates(p.x, p.y, p.z, cellX, cellY, cellZ);
      
      CollectedPoints& pointsInCell = points.points[cellX][cellY][cellZ];
      
      CollectedPoint cp;
      cp.x=p.x, cp.y=p.y, cp.z=p.z;
      //cp.nX=normal.x(), cp.nY=normal.y(), cp.nZ=normal.z();
      cp.range = range;
      cp.incidenceAngle = incidenceAngle;
      cp.intensity = p.intensity;
      cp.laser = p.laser;
      
      pointsInCell.push_back(cp);
      
      //if (useGui && getLocalTime()-lastGuiUpdate>0.1) {
        //lastGuiUpdate = getLocalTime();
        //qApplication->processEvents();
        //viewerPtr->updateGL();
      //}
    }
    
    //PointCloudType cloud2;
    for (std::map<int, std::map<int, std::map<int, PointDatabase> > >::const_iterator it=clusters.begin(); it!=clusters.end(); ++it) {
      int laser = it->first;
      float maxIntensity = maxIntensities[laser];
      float maxIntensityDeviation = parameters.maxIntensityDeviation * maxIntensity;
      //std::cout << "Laser: "<<laser<<"\n";
      for (std::map<int, std::map<int, PointDatabase> >::const_iterator it2=it->second.begin(); it2!=it->second.end(); ++it2) {
        //float range = cellHelper.cellToRange(it2->first);
        //std::cout << "Range: "<<range<<"m\n";
        for (std::map<int, PointDatabase>::const_iterator it3=it2->second.begin(); it3!=it2->second.end(); ++it3) {
          //float incidenceAngle = cellHelper.cellToIncidenceAngle(it3->first);
          //std::cout << "Incidence Angle: "<<incidenceAngle<<"m\n";
          const PointDatabase& points = it3->second;
          for (PointDatabase::ConstIteratorX itX=points.points.begin(); itX!=points.points.end(); ++itX) {
            for (PointDatabase::ConstIteratorY itY=itX->second.begin(); itY!=itX->second.end(); ++itY) {
              for (PointDatabase::ConstIteratorZ itZ=itY->second.begin(); itZ!=itY->second.end(); ++itZ) {
                float x, y, z;
                points.getCoordinatesFromIndices(itX->first, itY->first, itZ->first, x, y, z);
                //std::cout << "("<<x<<","<<y<<","<<z<<"): "<<itZ->second.size()<<" points\n";
                const CollectedPoints& pointsInCell = itZ->second;
                
                if (pointsInCell.size() < parameters.minNoOfMeasurementsForDataPointMean)
                  continue;
                
                CollectedPoint averagePoint;
                averagePoint.laser = laser;
                for (size_t pointsInCellIdx=0; pointsInCellIdx<pointsInCell.size(); ++pointsInCellIdx) {
                  const CollectedPoint& cp = pointsInCell[pointsInCellIdx];
                  averagePoint.x+=cp.x; averagePoint.y+=cp.y; averagePoint.z+=cp.z;
                  //averagePoint.nX+=cp.nX; averagePoint.nY+=cp.nY; averagePoint.nZ+=cp.nZ;
                  averagePoint.range+=cp.range;
                  averagePoint.incidenceAngle+=cp.incidenceAngle;
                  averagePoint.intensity+=cp.intensity;
                }
                
                averagePoint.x/=pointsInCell.size(); averagePoint.y/=pointsInCell.size(); averagePoint.z/=pointsInCell.size();
                //averagePoint.nX/=pointsInCell.size(); averagePoint.nY/=pointsInCell.size(); averagePoint.nZ/=pointsInCell.size();
                averagePoint.range/=pointsInCell.size();
                averagePoint.incidenceAngle/=pointsInCell.size();
                averagePoint.intensity/=pointsInCell.size();
                
                //if (averagePoint.range<0.1 || (averagePoint.laser>=100 && averagePoint.range>100))
                  //std::cerr << "Weird average range ("<<averagePoint.range<<"m for Laser "<<averagePoint.laser<<".\n";
                
                float intensityDeviation = 0.0f;
                for (size_t pointsInCellIdx=0; pointsInCellIdx<pointsInCell.size(); ++pointsInCellIdx) {
                  const CollectedPoint& cp = pointsInCell[pointsInCellIdx];
                  intensityDeviation += std::pow(averagePoint.intensity-cp.intensity, 2);
                }
                intensityDeviation = sqrtf(intensityDeviation/pointsInCell.size());
                //std::cout << PVARC(intensityDeviation);

                if (intensityDeviation > maxIntensityDeviation)
                  continue;
                
                //PointType p2;
                //p2.x=averagePoint.x, p2.y=averagePoint.y, p2.z=averagePoint.z;
                //p2.intensity=averagePoint.intensity, p2.laser=averagePoint.laser;
                //cloud2.push_back(p2);
                
                int cellX, cellY, cellZ;
                pointDatabase.getIndicesFromCoordinates(averagePoint.x, averagePoint.y, averagePoint.z, cellX, cellY, cellZ);
                CollectedPoints& pointsInMapCell = pointDatabase.points[cellX][cellY][cellZ];
                
                //if (intensityDeviation > maxIntensityDeviation) {
                  //pointsInMapCell.cellBad = true;
                  ////std::cout << "Markig cell as bad.\n";
                //}
                
                pointsInMapCell.push_back(averagePoint);
                
                //if (useGui) {
                  //pcl::PointXYZ mapCloudPoint;
                  //mapCloudPoint.x=averagePoint.x, mapCloudPoint.y=averagePoint.y, mapCloudPoint.z=averagePoint.z;
                  //mapCloud.push_back(mapCloudPoint);
                  //if (getLocalTime()-lastGuiUpdate>0.1) {
                    //lastGuiUpdate = getLocalTime();
                    //qApplication->processEvents();
                    //viewerPtr->updateGL();
                  //}
                //}
              }
            }
          }
        }
      }
    }
  }
  
  mapCells.clear();
  int goodCellsCounter=0, badCellsCounter=0;
  for (PointDatabase::IteratorX itX=pointDatabase.points.begin(); itX!=pointDatabase.points.end(); ++itX) {
    for (PointDatabase::IteratorY itY=itX->second.begin(); itY!=itX->second.end(); ++itY) {
      for (PointDatabase::IteratorZ itZ=itY->second.begin(); itZ!=itY->second.end(); ++itZ) {
        CollectedPoints& collectedPoints = itZ->second;
        if (collectedPoints.empty())
          continue;
        if (collectedPoints.cellBad) {
          ++badCellsCounter;
          //std::cout << "Skipping a bad cell.\n";
        }
        else {
          ++goodCellsCounter;
          mapCells.push_back(collectedPoints);
        }
        // tony added: FIXME: Segmentayion fault
        // itY->second.erase(itZ);  // Reduce necessary amount of RAM by not keeping already handled cells
      }
    }
  }
  std::cout << badCellsCounter<<" bad cells, "<<goodCellsCounter<<" good cells.\n";
  
  mapCells.removeExtremeIntensities(5.0f);
  
  if (!mapCellsFilename.empty())
    mapCells.writeToDisk(mapCellsFilename);
  std::cout << "\nWrote "<<mapCells.getNoOfDataPoints()<<" data points in "<<mapCells.size()<<" cells to file \""<<mapCellsFilename<<"\".\n";
  
  //mapCells.writeMapCloudToDisk();
  
  return true;
}

void RemissionCalibrationHelper::setCalibrationToInitialGuess(const CollectedPointsMap& mapCells) {
  remissionCalibrationResult.minMaxRanges = mapCells.getMinMaxRanges();
  remissionCalibrationResult.setMaxUncalibratedIntensities(mapCells.getMaxIntensities());
  std::map<int,float> avgIntensities = mapCells.getAverageIntensities();
  std::map<int,float> initialGuess;
  for (std::map<int,float>::iterator it=avgIntensities.begin(); it!=avgIntensities.end(); ++it)
    initialGuess[it->first] = 1.0f/(it->second);
  //std::cout << "\nInitial guess:\n";
  //for (std::map<int,float>::const_iterator it=initialGuess.begin(); it!=initialGuess.end(); ++it)
    //std::cout << it->first << ": " << it->second <<"\n";
  //std::cout << "\n";
  initializeCalibrationValues(initialGuess);
}

std::set<int> RemissionCalibrationHelper::getSkipLaserSet() {
  std::set<int> skipLaser;
  std::vector<std::set<int> > parsedLaserGroups = parseLaserGroupString(skipLasersString);
  for (size_t laserGroupIdx=0; laserGroupIdx<parsedLaserGroups.size(); ++laserGroupIdx) {
    const std::set<int>& lasers = parsedLaserGroups[laserGroupIdx];
    for (std::set<int>::const_iterator laserIt=lasers.begin(); laserIt!=lasers.end(); ++laserIt)
      skipLaser.insert(*laserIt);
  }
  return skipLaser;
}

float gaussianNoise() {
  const static int q = 15;
  const static double c1 = (1 << q) - 1;
  const static double c2 = ((int)(c1 / 3)) + 1;
  const static double c3 = 1.f / c1;
  return c3 * (2.0*c2*(double(rand())/(double(RAND_MAX)+1) +
                       double(rand())/(double(RAND_MAX)+1) +
                       double(rand())/(double(RAND_MAX)+1)) -
               3.0*(c2-1.0));
}

bool RemissionCalibrationHelper::writeFakeDataset(const std::string& folder, float standardDeviation) {
  std::vector<std::string> scan_file_names = getScanFilenames(folder);
  std::set<int> skipLaser = getSkipLaserSet();
  
  for (size_t pcdIdx=0; pcdIdx<scan_file_names.size(); ++pcdIdx) {
    std::string baseFilename = Ais3dTools::getPureFilename(scan_file_names[pcdIdx]);
    std::string pointCloudFilename = folder+"/"+scan_file_names[pcdIdx];
    
    pcl::PCLPointCloud2 pointCloudData;
    boost::shared_ptr<PointCloudType> cloudPtr(new PointCloudType);
    PointCloudType& cloud = *cloudPtr;
    if (!readCloud(pointCloudFilename, cloud, false))
      continue;
    
    NormalsMap normalsMap(0.2*parameters.mapCellSize);
    createSparseNormalsMap(cloud, normalsMap);
    
    PointCloudType cloud2 = cloud;
    for (size_t pointIdx=0; pointIdx<cloud2.size(); ++pointIdx) {
      PointType& point = cloud2[pointIdx];
      
      if (skipLaser.find(point.laser)!=skipLaser.end())
        point.x=point.y=point.z = std::numeric_limits<float>::quiet_NaN();
      
      if (!std::isfinite(point.x))
        continue;
      Eigen::Vector3f p = point.getVector3fMap(),
                      vp(point.vp_x, point.vp_y, point.vp_z);
      float range = (p-vp).norm();
      
      NormalWithQuality normalWithQuality;
      float cosIncidenceAngle = 1.0f;
      normalWithQuality = *normalsMap.cellFromPoint(p);
      if (normalWithQuality.normal.squaredNorm()<0.9)  // Invalid normal
        cosIncidenceAngle = 1.0f;
      else
        cosIncidenceAngle = getCosIncidenceAngle(cloud2[pointIdx], normalWithQuality.normal);
      float incidenceAngle = acosf(cosIncidenceAngle);
      
      point.intensity = 0.1*(point.laser+1) / ((0.1f*range+1.0f)*std::pow(incidenceAngle+1.0f, 2));
      //point.intensity = point.laser;
      //point.intensity = range;
      //point.intensity = incidenceAngle;
      //std::cout << PVARAC(incidenceAngle);
      
      if (standardDeviation > 1e-6) {
        point.intensity *= 1.0f+standardDeviation*gaussianNoise();
      }
    }
    std::string cloud2Filename = baseFilename+".pcd";
    pcl::io::savePCDFile(cloud2Filename, cloud2, true);
    std::cout << "Wrote file \""<<cloud2Filename<<"\".\n";
  }
  
  return true;
}


bool RemissionCalibrationHelper::calculateCalibrationBasedOnMapCells(const std::string& filename, bool skipInitialization) {
  double startTime = getLocalTime();
  
  CollectedPointsMap mapCells;
  if (!mapCells.readFromDisk(filename))
    return false;
  std::cout << "Performing optimization with "<<mapCells.getNoOfDataPoints()<<" measurements.\n";
  
  //bool useGui = qApplication!=NULL;
  //Ais3dTools::ALUGLViewer* viewerPtr = NULL;
  //PointCloudType points1, points2;
  //Ais3dTools::PclPointCloudObjectT<PointCloudType> visPointCloud;
  //visPointCloud.setDrawColor(0.0, 1.0, 0.0);
  
  //if (useGui) {
    //viewerPtr = new Ais3dTools::ALUGLViewer;
    //viewerPtr->show();
    //viewerPtr->updateGL();
    //viewerPtr->setBackgroundColor(Qt::white);
    ////viewerPtr->add(&visPointCloud);
  //}
  
  if (!skipInitialization)
    setCalibrationToInitialGuess(mapCells);
  
  std::stringstream calibrationFileName;
  calibrationFileName << "calibrationValues_mode"<< remissionCalibrationResult.mode << ".txt";
  
  if (parameters.useExpectationMaximizationMethod) {
    //if (useGui) {
      //qApplication->processEvents();
    //}
    int noOfIterations = parameters.noOfIterationsForEM;
    bool autoStop = false;
    if (noOfIterations < 0) {
      autoStop = true;
      noOfIterations = 100;
    }
    float lastRMSE = 1e10;
    for (int iteration=0; iteration<=noOfIterations; ++iteration) {
      CollectedPointsMap tmpMapCells = mapCells;
      tmpMapCells.setCellValuesToCalibratedMeans(remissionCalibrationResult);
      if (iteration>0 && parameters.noOfIterationsForG2O!=0) {
        // 最后一轮才进行InfluenceCalibration
        bool performIntensityInfluenceCalibration = noOfCellsZ>1 && iteration==(noOfIterations-1);
                                                    //iteration%1==0 &&
                                                    //iteration >= lrintf((3.0/4.0)*float(noOfIterations));
        bool removeOutliers = iteration > noOfIterations/2;
        //removeOutliers = false;  // HACK
        bool projectOntoCellMiddles = iteration > noOfIterations/2;
        //projectOntoCellMiddles = false;  // HACK
        
        if (projectOntoCellMiddles) {
          tmpMapCells.projectOntoCellMiddles(remissionCalibrationResult);
          tmpMapCells.setCellValuesToCalibratedMeans(remissionCalibrationResult);
        }
        
        if (removeOutliers) {
          //std::cout << PVARN(tmpMapCells.getRMSE(remissionCalibrationResult));
          //tmpMapCells.removeOutliers(parameters.outlierPercentage, remissionCalibrationResult);
          tmpMapCells.removeOutlierCells(parameters.outlierPercentage, remissionCalibrationResult);
          tmpMapCells.setCellValuesToCalibratedMeans(remissionCalibrationResult);
          //std::cout << PVARN(tmpMapCells.getRMSE(remissionCalibrationResult));
        }
        //std::cout << PVARN(tmpMapCells.getNoOfDataPoints());
        
        if (performIntensityInfluenceCalibration) {
          std::cout << "Running G2O for intensity influence calibration\n";
          runG2oOptimization(tmpMapCells, parameters.noOfIterationsForG2O, true);
          tmpMapCells.setCellValuesToCalibratedMeans(remissionCalibrationResult);
          std::cout << "\n";
        }
        
        std::cout << "Running G2O";
        if (projectOntoCellMiddles)
          std::cout << ", projecting data points to cell middles";
        if (removeOutliers)
          std::cout << ", removing "<<parameters.outlierPercentage<<"\% outliers";
        std::cout << ".\n";
        runG2oOptimization(tmpMapCells, parameters.noOfIterationsForG2O, false);
        tmpMapCells.setCellValuesToCalibratedMeans(remissionCalibrationResult);
        std::cout << "\n";
      }
      remissionCalibrationResult.normalizeCalibrationFactors(1.0f / tmpMapCells.getAverageCalibratedIntensity(remissionCalibrationResult));
      //remissionCalibrationResult.normalizeCalibrationFactors(1.0f / tmpMapCells.getMaxCalibratedIntensity(remissionCalibrationResult));
      tmpMapCells.setCellValuesToCalibratedMeans(remissionCalibrationResult);
      float currentRMSE = tmpMapCells.getRMSE(remissionCalibrationResult);
      std::cout << "RMSE of iteration "<<iteration<<" is "<<currentRMSE<<".\n\n";
      //remissionCalibrationResult.printCalibrationValues();
      std::cout << "Saving current calibration result to \""<<calibrationFileName.str()<<"\"\n";
      remissionCalibrationResult.saveCalibrationValues(calibrationFileName.str());
      
      if (autoStop && std::abs(lastRMSE-currentRMSE)<parameters.autoStopRmseChange) {
        autoStop = false;
        noOfIterations = iteration + 5;
        std::cout << "\n-----\nChange in RMSE below "<<parameters.autoStopRmseChange<<". "
                  << "Will stop after iteration "<<noOfIterations<<".\n-----\n\n";
      }
      lastRMSE = currentRMSE;
    }
  }
  mapCells.setCellValuesToCalibratedMeans(remissionCalibrationResult);
  remissionCalibrationResult.printCalibrationValues();
  std::cout << "Final RMSE is "<<mapCells.getRMSE(remissionCalibrationResult)<<".\n";
  
  double runtime = getLocalTime()-startTime;
  std::cout << "Runtime: "<<int(runtime)/60<<"min, "<<int(runtime)%60<<"sec.\n";
  return true;
}

bool RemissionCalibrationHelper::writeCalibratedClouds(const std::string& folder, const std::string& calibrationValuesFilename) {
  std::vector<std::string> scan_file_names = getScanFilenames(folder);
  
  //bool useGui = qApplication!=NULL;
  //Ais3dTools::ALUGLViewer* viewerPtr = NULL;
  //PointCloudType points1, points2;
  //Ais3dTools::PclPointCloudObjectT<PointCloudType> visPointCloud;
  //visPointCloud.setDrawColor(0.0, 1.0, 0.0);
  
  //if (useGui) {
    //setlocale(LC_NUMERIC, "C");
    //viewerPtr = new Ais3dTools::ALUGLViewer;
    //viewerPtr->show();
    //viewerPtr->updateGL();
    //viewerPtr->add(&visPointCloud);
    //viewerPtr->setBackgroundColor(Qt::white);
  //}
  
  if (!remissionCalibrationResult.readCalibrationValuesFromFile(calibrationValuesFilename)) {
    std::cout << "Could not read calibration values file.\n";
    return false;
  }

  //remissionCalibrationResult.printCalibrationValues();
  
  std::set<int> availableLaserIds = remissionCalibrationResult.getAvailableLaserIds();
  
  std::set<int> skipLaser = getSkipLaserSet();
  
  for (size_t pcdIdx=0; pcdIdx<scan_file_names.size(); ++pcdIdx) {
    std::string baseFilename = Ais3dTools::getPureFilename(scan_file_names[pcdIdx]);
    std::string pointCloudFilename = folder+"/"+scan_file_names[pcdIdx];
    
    pcl::PCLPointCloud2 pointCloudData;
    boost::shared_ptr<PointCloudType> cloudPtr(new PointCloudType);
    PointCloudType& cloud = *cloudPtr;
    if (!readCloud(pointCloudFilename, cloud, false))
      continue;
    
    //std::cout << "Read cloud "<<pcdIdx<<" of size "<<cloud.size()<<".\n";
    
    NormalsMap normalsMap(0.2*parameters.mapCellSize);
    createSparseNormalsMap(cloud, normalsMap);
    
    //bool cloudIsUnorganized = cloud.height<=1;
    //boost::shared_ptr<pcl::KdTreeFLANN<PointType> > kdTreePtr(new pcl::KdTreeFLANN<PointType>);
    //boost::shared_ptr<PointCloudType> sparseCloudPtr(new PointCloudType);
    //if (cloudIsUnorganized && remissionCalibrationResult.noOfCellsY>1) {
      //float sparseResolution = 0.2 * parameters.mapCellSize;
      //createSparsePointCloud(cloud, *sparseCloudPtr, sparseResolution);
      //kdTreePtr->setInputCloud(sparseCloudPtr);
    //}
    
    PointCloudType cloud2 = cloud;
    for (size_t pointIdx=0; pointIdx<cloud2.size(); ++pointIdx) {
      PointType& point = cloud2[pointIdx];
      
      if (availableLaserIds.find(point.laser)==availableLaserIds.end() || skipLaser.find(point.laser)!=skipLaser.end())
        point.x=point.y=point.z = std::numeric_limits<float>::quiet_NaN();
      
      if (!std::isfinite(point.x))
        continue;
      Eigen::Vector3f p = point.getVector3fMap(),
                      vp(point.vp_x, point.vp_y, point.vp_z);
      float range = (p-vp).norm();
      
      NormalWithQuality normalWithQuality;
      float cosIncidenceAngle = 1.0f;
      if (remissionCalibrationResult.noOfCellsY > 1) {
        normalWithQuality = *normalsMap.cellFromPoint(p);
        //if (cloudIsUnorganized)  // Unorganized cloud
          //getNormal(point, *kdTreePtr, normal, 0.5f*parameters.mapCellSize);
        //else  // Organized cloud
          //getNormal(cloud, pointIdx, normal, pointIdx==0);
        if (normalWithQuality.normal.squaredNorm()<0.9)  // Invalid normal
          cosIncidenceAngle = 1.0f;
        else
          cosIncidenceAngle = getCosIncidenceAngle(cloud2[pointIdx], normalWithQuality.normal);
      }
      float incidenceAngle = acosf(cosIncidenceAngle);
      //std::cout << PVARAC(incidenceAngle);
      
      point.intensity = remissionCalibrationResult.getCalibratedIntensity(range, incidenceAngle, point.intensity, point.laser);
      //point.intensity = normalWithQuality.quality;  // HACK
      //point.intensity = incidenceAngle;  // HACK
      //point.intensity = point.laser;  // HACK

      //std::cout << point.intensity << " " << std::flush;

      //point.intensity = polynomial.getValue(squaredRange, cosIncidenceAngle) * point.intensity;
    }
    //std::string cloud2Filename = baseFilename+"_normalizedIntensities.pcd";
    std::string cloud2Filename = baseFilename+".pcd";
    pcl::io::savePCDFile(cloud2Filename, cloud2, true);
    std::cout << "Wrote file \""<<cloud2Filename<<"\".\n";
  }
  
  return true;
}


//template <typename real>
//struct SamplePointT {
  //real x1, y1, z1,
       //x2, y2, z2;
//};

//template <typename real>
//bool bivariatePolynomialApproximation(
    //const std::vector<SamplePointT<real> >& samplePoints, unsigned int polynomial_degree, unsigned int fixedParameter, float fixedValue, BivariatePolynomialT<real>& ret)
//{
  //MEASURE_FUNCTION_TIME;
  
  //unsigned int parameters_size = BivariatePolynomialT<real>::getNoOfParametersFromDegree(polynomial_degree);
  ////cout << PVARN (parameters_size);

  ////cout << "Searching for the "<<parameters_size<<" parameters for the bivariate polynom of degree "
  ////     << polynomial_degree<<" using "<<samplePoints.size ()<<" points.\n";
  
  //if (parameters_size-1 > samplePoints.size ()) // Too many parameters for this number of equations (points)?
  //{
    //return false;    
    //// Reduce degree of polynomial
    ////polynomial_degree = (unsigned int) (0.5f* (sqrtf (8*samplePoints.size ()+1) - 3));
    ////parameters_size = BivariatePolynomialT<real>::getNoOfParametersFromDegree (polynomial_degree);
    ////cout << "Not enough points, so degree of polynomial was decreased to "<<polynomial_degree
    ////     << " ("<<samplePoints.size ()<<" points => "<<parameters_size<<" parameters)\n";
  //}
  
  //ret.setDegree (polynomial_degree);
  
  ////double coeffStuffStartTime=-get_time ();
  ////Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> A (parameters_size, parameters_size);
  ////A.setZero();
  ////Eigen::Matrix<real, Eigen::Dynamic, 1> b (parameters_size);
  ////b.setZero();
  //Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> A (samplePoints.size(), parameters_size-1);
  //Eigen::Matrix<real, Eigen::Dynamic, 1> b (samplePoints.size());
  //std::vector<real> tmpC1, tmpC2;
  //tmpC1.resize(parameters_size);
  //tmpC2.resize(parameters_size);
  //real* tmpCEndPtr1 = &tmpC1.back(),
      //* tmpCEndPtr2 = &tmpC2.back();
  //for (size_t samplePointIdx=0; samplePointIdx<samplePoints.size(); ++samplePointIdx) {
    //const SamplePointT<real>& samplePoint = samplePoints[samplePointIdx];
    //real currentX1=samplePoint.x1, currentY1=samplePoint.y1, currentZ1=samplePoint.z1,
         //currentX2=samplePoint.x2, currentY2=samplePoint.y2, currentZ2=samplePoint.z2;
    //unsigned int posInC = parameters_size;
    //real* tmpCPtr1 = tmpCEndPtr1,
        //* tmpCPtr2 = tmpCEndPtr2;
    //real tmpX1=1, tmpX2=1;
    //for (unsigned int xDegree=0; xDegree<=polynomial_degree; ++xDegree) {
      //real tmpY1=1, tmpY2=1;
      //for (unsigned int yDegree=0; yDegree<=polynomial_degree-xDegree; ++yDegree) {
        //--posInC;
        //*(tmpCPtr1--) = tmpX1*tmpY1;
        //*(tmpCPtr2--) = tmpX2*tmpY2;
        ////cout << "x="<<currentX<<", y="<<currentY<<", Pos "<<posInC--<<": "<<tmpX<<"*"<<tmpY<<"="<<tmpC[posInC]<<"\n";
        //tmpY1 *= currentY1;
        //tmpY2 *= currentY2;
      //}
      //tmpX1 *= currentX1;
      //tmpX2 *= currentX2;
    //}
    
    //unsigned int column=0;
    //for (unsigned int parameterIdx=0; parameterIdx<parameters_size; ++parameterIdx) {
      //if (parameterIdx!=fixedParameter)
        //A(samplePointIdx, column++) = currentZ1*tmpC1[parameterIdx]-currentZ2*tmpC2[parameterIdx];
    //}
    //b[samplePointIdx] = fixedValue*(tmpC2[fixedParameter]*currentZ2-tmpC1[fixedParameter]*currentZ1);
    
    ////real* APtr = &A(0,0);
    ////real* bPtr = &b[0];
    ////real* tmpCPtr1=tmpC;
    ////for (unsigned int i=0; i<parameters_size; ++i)
    ////{
      ///[> (bPtr++) += currentZ * *tmpCPtr1;
      
      ////real* tmpCPtr2=tmpC;
      ////for (unsigned int j=0; j<parameters_size; ++j)
      ////{
        ///[> (APtr++) += *tmpCPtr1 * * (tmpCPtr2++);
      ////}
      
      ////++tmpCPtr1;
    ////}
    ////A += DMatrix<real>::outProd (tmpC);
    ////b += currentZ * tmpC;
  //}
  ////cout << "Calculating matrix A and vector b (size "<<b.size ()<<") from "<<samplePoints.size ()<<" points took "
       ////<< (coeffStuffStartTime+get_time ())*1000<<"ms using constant memory.\n";
    ////cout << PVARC (A)<<PVARN (b);
  
  //Eigen::Matrix<real, Eigen::Dynamic, 1> parameters;
  //parameters =(A.transpose()*A).inverse()*A.transpose() * b;
  
  //unsigned int parameterIdx=0;
  //for (unsigned int i=0; i<parameters_size-1; i++) {
    //if (i==fixedParameter)
      //ret.parameters[parameterIdx++] = fixedValue;
    //ret.parameters[parameterIdx++] = parameters[i];
  //}
  
  //BivariatePolynomialT<real> ret_ = ret;
  //for (unsigned int i=0; i<parameters_size; i++)
    //ret_.parameters[i] = 0.0f;
  //ret_.parameters[fixedParameter] = fixedValue;
  
  //cout << "Resulting polynomial is "<<ret<<"\n";
  //float squaredError=0.0f, otherSquaredError=0.0f;
  //for (size_t samplePointIdx=0; samplePointIdx<samplePoints.size(); ++samplePointIdx) {
    //const SamplePointT<real>& samplePoint = samplePoints[samplePointIdx];
    ////if (std::abs(sqrtf(samplePoint.x1)-sqrtf(samplePoint.x2))<1)
      ////continue;
    //float z1=ret.getValue(samplePoint.x1, samplePoint.y1)*samplePoint.z1,
          //z2=ret.getValue(samplePoint.x2, samplePoint.y2)*samplePoint.z2;
    //float z1_=ret_.getValue(samplePoint.x1, samplePoint.y1)*samplePoint.z1,
          //z2_=ret_.getValue(samplePoint.x2, samplePoint.y2)*samplePoint.z2;
    //std::cout << PVARC(samplePoint.x1) << PVARC(samplePoint.x2)
              //<< PVARC(samplePoint.y1) << PVARC(samplePoint.y2)
              //<< PVARC(samplePoint.z1) << PVARC(samplePoint.z2)
              //<< PVARC(z1) << PVARC(z2)
              //<< PVARC(z1_) << PVARN(z2_);
    //squaredError += std::pow(z2-z1, 2);
    //otherSquaredError += std::pow(z2_-z1_, 2);
  //}
  //std::cout << PVARC(squaredError)<<PVARN(otherSquaredError);
  
  //return true;
//}

std::vector<std::set<int> > RemissionCalibrationHelper::parseLaserGroupString(const std::string& inputString) {
  //std::cout << "Parsing laser groups string \""<<inputString<<"\".\n";
  std::vector<std::set<int> > ret;
  
  stringstream ss(inputString);
  std::string groupString;
  bool inBrackets = false;
  while(std::getline(ss, groupString, ',')) {
    //std::cout << PVARN(groupString) << '\n';
    if (!inBrackets)
      ret.push_back(std::set<int>());
    std::set<int>& lg = ret.back();
    //int laserGroupIdx = ret.size()-1;
    while (!groupString.empty()) {
      if (groupString[0]=='{') {
        if (inBrackets) {
          std::cerr << "\n"<<__PRETTY_FUNCTION__ << ": Recursive brackets not allowed.\n\n";
          exit(1);
        }
        inBrackets = true;
        groupString.erase(groupString.begin());
        continue;
      }
      if (groupString[0]=='}') {
        if (!inBrackets) {
          std::cerr << "\n"<<__PRETTY_FUNCTION__ << ": Bracket closed but not opened.\n\n";
          exit(1);
        }
        inBrackets = false;
        groupString.erase(groupString.begin());
        continue;
      }
      size_t numberEnd = groupString.find_first_not_of("0123456789", 0);
      if (numberEnd==0) {
        std::cerr << "\n"<<__PRETTY_FUNCTION__ << ": Expecting a number here \""<<groupString<<"\".\n\n";
        exit(1);
      }
      int laserId = atoi(groupString.substr(0,numberEnd).c_str());
      groupString.erase(groupString.begin(),(numberEnd==std::string::npos ? groupString.end() : groupString.begin()+numberEnd));
      if (!groupString.empty() && groupString[0]=='-') {
        groupString.erase(groupString.begin());
        numberEnd = groupString.find_first_not_of("0123456789", 0);
        if (numberEnd==0) {
          std::cerr << "\n"<<__PRETTY_FUNCTION__ << ": Expecting a number after '-'.\n\n";
          exit(1);
        }
        int laserId2 = atoi(groupString.substr(0,numberEnd).c_str());
        groupString.erase(groupString.begin(),(numberEnd==std::string::npos ? groupString.end() : groupString.begin()+numberEnd));
        //std::cout << "Adding laser ids from "<<laserId<<" to "<<laserId2<<" to laser group "<<laserGroupIdx<<".\n";
        for (int laserIdToAdd=laserId; laserIdToAdd<=laserId2; ++laserIdToAdd) {
          lg.insert(laserIdToAdd);
        }
      }
      else {
        //std::cout << "Adding laser id "<<laserId<<" to laser group "<<laserGroupIdx<<".\n";
        lg.insert(laserId);
      }
    }
  }
  return ret;
}

void RemissionCalibrationHelper::initializeCalibrationValues(const std::set<int>& laserIds) {
  std::map<int,float> initialValuesPerLaser;
  for (std::set<int>::const_iterator it=laserIds.begin(); it!=laserIds.end(); ++it)
    initialValuesPerLaser[*it] = initialValue;
  initializeCalibrationValues(initialValuesPerLaser);
}

void RemissionCalibrationHelper::initializeCalibrationValues(const std::map<int,float>& initialValuesPerLaser) {
  calibrationValuesMode1.clear();
  calibrationValuesMode2.clear();
  calibrationValuesMode3.clear();
  if (mode==1) {
    for (std::map<int,float>::const_iterator it=initialValuesPerLaser.begin(); it!=initialValuesPerLaser.end(); ++it) {
      calibrationValuesMode1[it->first].resize(noOfCellsX*noOfCellsY, it->second);
    }
  }
  else if (mode==2) {
    for (std::map<int,float>::const_iterator it=initialValuesPerLaser.begin(); it!=initialValuesPerLaser.end(); ++it) {
      calibrationValuesMode2[it->first].calibrationValuesRange.resize(noOfCellsX, it->second);
      calibrationValuesMode2[it->first].calibrationValuesIncidenceAngle.resize(noOfCellsY, initialValue);
      calibrationValuesMode2[it->first].calibrationValuesIntensity.resize(noOfCellsZ, initialValue);
    }
  }
  else {
    std::vector<std::set<int> > parsedLaserGroups = parseLaserGroupString(laserGroupString);
    calibrationValuesMode3.resize(parsedLaserGroups.size());
    std::set<int> remainingLaserIds;
    for (std::map<int,float>::const_iterator it=initialValuesPerLaser.begin(); it!=initialValuesPerLaser.end(); ++it)
      remainingLaserIds.insert(it->first);
    
    for (size_t laserGroupIdx=0; laserGroupIdx<parsedLaserGroups.size(); ++laserGroupIdx) {
      RemissionCalibrationResult::LaserGroup& lg = calibrationValuesMode3[laserGroupIdx];
      const std::set<int>& lasers = parsedLaserGroups[laserGroupIdx];
      for (std::set<int>::const_iterator laserIt=lasers.begin(); laserIt!=lasers.end(); ++laserIt) {
        std::map<int,float>::const_iterator initialValueIt = initialValuesPerLaser.find(*laserIt);
        if (initialValueIt==initialValuesPerLaser.end()) {
          std::cerr << "Laser "<<*laserIt<<" given in laser group string but not available in dataset.\n";
          lg.calibrationValuesLaser[*laserIt] = initialValue;
        }
        else
          lg.calibrationValuesLaser[*laserIt] = initialValueIt->second;
        remainingLaserIds.erase(*laserIt);
      }
    }
    
    if (!remainingLaserIds.empty()) {
      calibrationValuesMode3.push_back(RemissionCalibrationResult::LaserGroup());
      for (std::set<int>::const_iterator it=remainingLaserIds.begin(); it!=remainingLaserIds.end(); ++it) {
        calibrationValuesMode3.back().calibrationValuesLaser[*it] = initialValuesPerLaser.find(*it)->second;
      }
    }
    
    for (size_t laserGroupIdx=0; laserGroupIdx<calibrationValuesMode3.size(); ++laserGroupIdx) {
      RemissionCalibrationResult::LaserGroup& lg = calibrationValuesMode3[laserGroupIdx];
      lg.calibrationValuesRange.resize(noOfCellsX, initialValue);
      lg.calibrationValuesIncidenceAngle.resize(noOfCellsY, initialValue);
      lg.calibrationValuesIntensity.resize(noOfCellsZ, initialValue);
    }
    
    calibrationValuesMode3.updateLaserIndices();
    remissionCalibrationResult.updateMaxUncalibratedIntensitiesForLaserGroups();
  }
  //remissionCalibrationResult.normalizeCalibrationFactors();
  //remissionCalibrationResult.printCalibrationValues();
}

void RemissionCalibrationHelper::CollectedPointsMap::writeCalibratedMapCloudToDisk(const RemissionCalibrationResult& calibration) const {
  pcl::PointCloud<pcl::PointXYZI> mapCloud;
  for (size_t mapCellIdx=0; mapCellIdx<size(); ++mapCellIdx) {
    pcl::PointXYZI p = at(mapCellIdx).getPointRepresentingCell();
    p.intensity = at(mapCellIdx).getAverageCalibratedIntensity(calibration);
    mapCloud.push_back(p);
  }
  mapCloud.width=mapCloud.size(), mapCloud.height=1;
  pcl::io::savePCDFile("mapCloudCalibrated.pcd", mapCloud, true);
  
  mapCloud.clear();
  for (size_t mapCellIdx=0; mapCellIdx<size(); ++mapCellIdx) {
    pcl::PointXYZI p = at(mapCellIdx).getPointRepresentingCell();
    p.intensity = at(mapCellIdx).getRMSE(calibration);
    mapCloud.push_back(p);
  }
  mapCloud.width=mapCloud.size(), mapCloud.height=1;
  pcl::io::savePCDFile("mapCloudErrors.pcd", mapCloud, true);
}


void RemissionCalibrationHelper::CollectedPointsMap::writeMapCloudToDisk() const {
  pcl::PointCloud<pcl::PointXYZI> mapCloud;
  for (size_t mapCellIdx=0; mapCellIdx<size(); ++mapCellIdx) {
    const CollectedPoints& collectedPoints = at(mapCellIdx);
    for (size_t pointIdx=0; pointIdx<collectedPoints.size(); ++pointIdx) {
      const CollectedPoint& cp = collectedPoints[pointIdx];
      pcl::PointXYZI p;
      p.x = cp.x;
      p.y = cp.y;
      p.z = cp.z;
      p.intensity = cp.intensity;
      mapCloud.points.push_back(p);
    }
  }
  mapCloud.width = mapCloud.size();
  mapCloud.height = 1;
  pcl::io::savePCDFile("mapCloud.pcd", mapCloud, true);
}

void RemissionCalibrationHelper::CollectedPointsMap::writeToDisk(const std::string& filename) const {
  std::ofstream file(filename.c_str());
  size_t noOfMapCells = size();
  file.write(reinterpret_cast<const char*>(&noOfMapCells), sizeof(noOfMapCells));
  for (size_t mapCellIdx=0; mapCellIdx<size(); ++mapCellIdx) {
    const CollectedPoints& collectedPoints = at(mapCellIdx);
    size_t noOfPointsInCell = collectedPoints.size();
    file.write(reinterpret_cast<const char*>(&noOfPointsInCell), sizeof(noOfPointsInCell));
    file.write(reinterpret_cast<const char*>(&collectedPoints[0]), noOfPointsInCell*sizeof(collectedPoints[0]));
  }
}

bool RemissionCalibrationHelper::CollectedPointsMap::readFromDisk(const std::string& filename) {
  std::ifstream file(filename.c_str());
  if (!file)
    return false;
  
  std::set<int> skipLaser = getSkipLaserSet();

  size_t noOfMapCells;
  file.read(reinterpret_cast<char*>(&noOfMapCells), sizeof(noOfMapCells));
  resize(noOfMapCells);
  for (size_t mapCellIdx=0; mapCellIdx<size(); ++mapCellIdx) {
    CollectedPoints& collectedPoints = at(mapCellIdx);
    size_t noOfPointsInCell;
    file.read(reinterpret_cast<char*>(&noOfPointsInCell), sizeof(noOfPointsInCell));
    collectedPoints.resize(noOfPointsInCell);
    file.read(reinterpret_cast<char*>(&collectedPoints[0]), noOfPointsInCell*sizeof(collectedPoints[0]));
    
    if (!skipLaser.empty()) {
      for (CollectedPoints::iterator it=collectedPoints.begin(); it!=collectedPoints.end();) {
        if (skipLaser.find(it->laser)==skipLaser.end())
          ++it;
        else
          it = collectedPoints.erase(it);
      }
    }
  }
  return true;
}

std::set<int> RemissionCalibrationHelper::CollectedPointsMap::getLaserIds() const {
  std::set<int> ret;
  for (size_t mapCellIdx=0; mapCellIdx<size(); ++mapCellIdx) {
    const CollectedPoints& collectedPoints = at(mapCellIdx);
    for (size_t pointIdx=0; pointIdx<collectedPoints.size(); ++pointIdx) {
      const CollectedPoint& cp = collectedPoints[pointIdx];
      ret.insert(cp.laser);
    }
  }
  //std::cout << "\nFound laser ids ";
  //for (std::set<int>::const_iterator it=ret.begin(); it!=ret.end(); ++it) {
    //std::cout << *it<<" ";
  //}
  //std::cout << " in map.\n";
  return ret;
}

std::map<int,std::pair<float,float> > RemissionCalibrationHelper::CollectedPointsMap::getMinMaxRanges() const {
  std::map<int,std::pair<float,float> > ret;
  for (size_t mapCellIdx=0; mapCellIdx<size(); ++mapCellIdx) {
    const CollectedPoints& collectedPoints = at(mapCellIdx);
    for (size_t pointIdx=0; pointIdx<collectedPoints.size(); ++pointIdx) {
      const CollectedPoint& cp = collectedPoints[pointIdx];
      if (ret.find(cp.laser)==ret.end())
        ret[cp.laser] = std::make_pair(cp.range,cp.range);
      else {
        ret[cp.laser].first  = std::min(ret[cp.laser].first, cp.range);
        ret[cp.laser].second = std::max(ret[cp.laser].second, cp.range);
      }
    }
  }
  //std::cout << "\nMinimum and maximum ranges:\n";
  //for (std::map<int,std::pair<float,float> >::const_iterator it=ret.begin(); it!=ret.end(); ++it)
    //std::cout << it->first<<": "<<it->second.first<<"-"<<it->second.second<<"\n";
  return ret;
}

std::map<int,float> RemissionCalibrationHelper::CollectedPointsMap::getMaxIntensities() const {
  std::map<int,float> ret;
  for (size_t mapCellIdx=0; mapCellIdx<size(); ++mapCellIdx) {
    const CollectedPoints& collectedPoints = at(mapCellIdx);
    for (size_t pointIdx=0; pointIdx<collectedPoints.size(); ++pointIdx) {
      const CollectedPoint& cp = collectedPoints[pointIdx];
      if (ret.find(cp.laser)==ret.end())
        ret[cp.laser] = cp.intensity;
      else
        ret[cp.laser] = std::max(ret[cp.laser], cp.intensity);
    }
  }
  //std::cout << "\nMaximum intensities:\n";
  //for (std::map<int,float>::const_iterator it=ret.begin(); it!=ret.end(); ++it)
    //std::cout << it->first<<": "<<it->second<<"\n";
  return ret;
}

float RemissionCalibrationHelper::CollectedPointsMap::getMaxIntensity() const {
  float ret = 0.0f;
  for (size_t mapCellIdx=0; mapCellIdx<size(); ++mapCellIdx) {
    const CollectedPoints& collectedPoints = at(mapCellIdx);
    for (size_t pointIdx=0; pointIdx<collectedPoints.size(); ++pointIdx) {
      const CollectedPoint& cp = collectedPoints[pointIdx];
      ret = std::max(ret, cp.intensity);
    }
  }
  //std::cout << "\nMaximum intensity is "<<ret<<".\n";
  return ret;
}

std::map<int,float> RemissionCalibrationHelper::CollectedPointsMap::getAverageIntensities() const {
  std::map<int,float> ret;
  std::map<int,float> weights;
  for (size_t mapCellIdx=0; mapCellIdx<size(); ++mapCellIdx) {
    const CollectedPoints& collectedPoints = at(mapCellIdx);
    for (size_t pointIdx=0; pointIdx<collectedPoints.size(); ++pointIdx) {
      const CollectedPoint& cp = collectedPoints[pointIdx];
      if (ret.find(cp.laser)==ret.end()) {
        ret[cp.laser] = cp.intensity;
        weights[cp.laser] = 1.0f;
      }
      else {
        ret[cp.laser] += cp.intensity;
        weights[cp.laser] += 1.0f;
      }
    }
  }
  for (std::map<int,float>::iterator it=ret.begin(); it!=ret.end(); ++it)
    it->second /= weights[it->first];
  
  //std::cout << "\nAverage intensities:\n";
  //for (std::map<int,float>::const_iterator it=ret.begin(); it!=ret.end(); ++it)
    //std::cout << it->first<<": "<<it->second<<"\n";
  
  return ret;
}

float RemissionCalibrationHelper::CollectedPointsMap::getAverageIntensity() const {
  float ret=0.0f, weightSum=0.0f;
  for (size_t mapCellIdx=0; mapCellIdx<size(); ++mapCellIdx) {
    const CollectedPoints& collectedPoints = at(mapCellIdx);
    for (size_t pointIdx=0; pointIdx<collectedPoints.size(); ++pointIdx) {
      const CollectedPoint& cp = collectedPoints[pointIdx];
      ret += cp.intensity;
      weightSum += 1.0f;
    }
  }
  ret /= weightSum;
  //std::cout << "\nAverage intensity is "<<ret<<".\n";
  return ret;
}

float RemissionCalibrationHelper::CollectedPointsMap::getAverageCalibratedIntensity(const RemissionCalibrationResult& calibration) const {
  float ret=0.0f, weightSum=0.0f;
  for (size_t mapCellIdx=0; mapCellIdx<size(); ++mapCellIdx) {
    const CollectedPoints& collectedPoints = at(mapCellIdx);
    for (size_t pointIdx=0; pointIdx<collectedPoints.size(); ++pointIdx) {
      const CollectedPoint& cp = collectedPoints[pointIdx];
      ret += calibration.getCalibratedIntensity(cp.range, cp.incidenceAngle, cp.intensity, cp.laser);
      weightSum += 1.0f;
    }
  }
  ret /= weightSum;
  //std::cout << "\nAverage calibrated intensity is "<<ret<<".\n";
  return ret;
}

void RemissionCalibrationHelper::CollectedPointsMap::removeExtremeIntensities(float percentage) {
  std::map<int, std::vector<float> > intensitiesPerLaser;
  for (size_t mapCellIdx=0; mapCellIdx<size(); ++mapCellIdx) {
    const CollectedPoints& collectedPoints = at(mapCellIdx);
    for (size_t pointIdx=0; pointIdx<collectedPoints.size(); ++pointIdx) {
      const CollectedPoint& cp = collectedPoints[pointIdx];
      intensitiesPerLaser[cp.laser].push_back(cp.intensity);
    }
  }
  std::map<int, std::pair<float,float> > minMaxAllowedIntensities;
  for (std::map<int, std::vector<float> >::iterator it=intensitiesPerLaser.begin(); it!=intensitiesPerLaser.end(); ++it) {
    int laser = it->first;
    std::vector<float>& intensities = it->second;
    std::sort(intensities.begin(), intensities.end());
    minMaxAllowedIntensities[laser].first  = intensities[lrintf(float(intensities.size()-1)*percentage/100.0f)];
    minMaxAllowedIntensities[laser].second = intensities[lrintf(float(intensities.size()-1)*(1.0f-percentage/100.0f))];
    //std::cout << "Cropping intensities for laser "<<laser<<": ["
              //<< minMaxAllowedIntensities[laser].first<<","<<minMaxAllowedIntensities[laser].second<<"].\n";
  }
  for (size_t mapCellIdx=0; mapCellIdx<size(); ++mapCellIdx) {
    CollectedPoints& collectedPoints = at(mapCellIdx);
    for (CollectedPoints::iterator it=collectedPoints.begin(); it!=collectedPoints.end();) {
      if (it->intensity<minMaxAllowedIntensities[it->laser].first || it->intensity>minMaxAllowedIntensities[it->laser].second)
        it = collectedPoints.erase(it);
      else
        ++it;
    }
  }
}

void RemissionCalibrationHelper::CollectedPointsMap::removeOutliers(float percentage, const RemissionCalibrationResult& calibration) {
  std::vector<float> errors;
  for (size_t mapCellIdx=0; mapCellIdx<size(); ++mapCellIdx) {
    const CollectedPoints& collectedPoints = at(mapCellIdx);
    for (size_t pointIdx=0; pointIdx<collectedPoints.size(); ++pointIdx) {
      const CollectedPoint& cp = collectedPoints[pointIdx];
      errors.push_back(std::abs(collectedPoints.cellValue-calibration.getCalibratedIntensity(cp.range, cp.incidenceAngle, cp.intensity, cp.laser)));
    }
  }
  std::sort(errors.begin(), errors.end());
  float maxAllowedError = errors[lrintf(float(errors.size()-1)*(1.0f-percentage/100.0f))];
  //std::cout << "Removing all measurements with an error above "<<maxAllowedError<<".\n";
  for (size_t mapCellIdx=0; mapCellIdx<size(); ++mapCellIdx) {
    CollectedPoints& collectedPoints = at(mapCellIdx);
    for (CollectedPoints::iterator it=collectedPoints.begin(); it!=collectedPoints.end();) {
      float error = std::abs(collectedPoints.cellValue-calibration.getCalibratedIntensity(it->range, it->incidenceAngle, it->intensity, it->laser));
      if (error > maxAllowedError)
        it = collectedPoints.erase(it);
      else
        ++it;
    }
  }
}

std::map<int, float> RemissionCalibrationHelper::CollectedPointsMap::getErrorPercentages(const RemissionCalibrationResult& calibration) const {
  std::map<int, float> ret;
  std::vector<float> errors;
  for (size_t mapCellIdx=0; mapCellIdx<size(); ++mapCellIdx) {
    const CollectedPoints& collectedPoints = at(mapCellIdx);
    for (size_t pointIdx=0; pointIdx<collectedPoints.size(); ++pointIdx) {
      const CollectedPoint& cp = collectedPoints[pointIdx];
      errors.push_back(std::abs(collectedPoints.cellValue-calibration.getCalibratedIntensity(cp.range, cp.incidenceAngle, cp.intensity, cp.laser)));
    }
  }
  std::sort(errors.begin(), errors.end());
  for (int percentage=1; percentage<=100; ++percentage)
    ret[percentage] = errors[lrintf(float(errors.size()-1)*(float(percentage)/100.0f))];
  return ret;
}

void RemissionCalibrationHelper::CollectedPointsMap::removeOutlierCells(float percentage, const RemissionCalibrationResult& calibration) {
  
  std::vector<float> errors;
  for (size_t mapCellIdx=0; mapCellIdx<size(); ++mapCellIdx)
    errors.push_back(at(mapCellIdx).getRMSE(calibration));
  std::sort(errors.begin(), errors.end());
  
  //for (int percent=0; percent<=100; ++percent)
    //std::cout << percent<<"%: "<<errors[lrintf(float(errors.size()-1)*(percent/100.0f))]<<"\n";
  
  float maxAllowedError = errors[lrintf(float(errors.size()-1)*(1.0f-percentage/100.0f))];
  //std::cout << "Removing all measurements with an error above "<<maxAllowedError<<".\n";
  
  RemissionCalibrationHelper::CollectedPointsMap tmp;
  for (size_t mapCellIdx=0; mapCellIdx<size(); ++mapCellIdx) {
    if (at(mapCellIdx).getRMSE(calibration) <= maxAllowedError)
      tmp.push_back(at(mapCellIdx));
  }
  std::swap(tmp, *this);
  
  //errors.clear();
  //for (size_t mapCellIdx=0; mapCellIdx<size(); ++mapCellIdx)
    //errors.push_back(at(mapCellIdx).getRMSE(calibration));
  //std::sort(errors.begin(), errors.end());
  //for (int percent=0; percent<=100; ++percent)
    //std::cout << percent<<"%: "<<errors[lrintf(float(errors.size()-1)*(percent/100.0f))]<<"\n";
}

float RemissionCalibrationHelper::CollectedPointsMap::getMaxCalibratedIntensity(const RemissionCalibrationResult& calibration) const {
  float ret=0.0f;
  for (size_t mapCellIdx=0; mapCellIdx<size(); ++mapCellIdx) {
    const CollectedPoints& collectedPoints = at(mapCellIdx);
    for (size_t pointIdx=0; pointIdx<collectedPoints.size(); ++pointIdx) {
      const CollectedPoint& cp = collectedPoints[pointIdx];
      ret = std::max(ret, calibration.getCalibratedIntensity(cp.range, cp.incidenceAngle, cp.intensity, cp.laser));
    }
  }
  //std::cout << "\nMaximum calibrated intensity is "<<ret<<".\n";
  return ret;
}

void RemissionCalibrationHelper::CollectedPointsMap::calibrateIntensities(const RemissionCalibrationResult& calibration) {
  std::set<int> availableLaserIds = calibration.getAvailableLaserIds();
  for (size_t mapCellIdx=0; mapCellIdx<size(); ++mapCellIdx) {
    CollectedPoints& collectedPoints = at(mapCellIdx);
    for (size_t pointIdx=0; pointIdx<collectedPoints.size(); ++pointIdx) {
      CollectedPoint& cp = collectedPoints[pointIdx];
      if (availableLaserIds.find(cp.laser) != availableLaserIds.end()) {
        cp.intensity = calibration.getCalibratedIntensity(cp.range, cp.incidenceAngle, cp.intensity, cp.laser);
      }
      else {
        //...
      }
    }
  }
}

void RemissionCalibrationHelper::CollectedPointsMap::normalizeIntensities(const std::map<int,float>& normalizationValues, bool reciprocal) {
  std::map<int,float> normalizationValues2 = normalizationValues;
  if (reciprocal) {
    for (std::map<int,float>::iterator it=normalizationValues2.begin(); it!=normalizationValues2.end(); ++it)
      it->second = 1.0f/it->second;
  }
  
  for (size_t mapCellIdx=0; mapCellIdx<size(); ++mapCellIdx) {
    CollectedPoints& collectedPoints = at(mapCellIdx);
    for (size_t pointIdx=0; pointIdx<collectedPoints.size(); ++pointIdx) {
      CollectedPoint& cp = collectedPoints[pointIdx];
      cp.intensity *= normalizationValues2[cp.laser];
    }
  }
}

void RemissionCalibrationHelper::CollectedPointsMap::normalizeIntensities(float normalizationValue) {
  for (size_t mapCellIdx=0; mapCellIdx<size(); ++mapCellIdx) {
    CollectedPoints& collectedPoints = at(mapCellIdx);
    for (size_t pointIdx=0; pointIdx<collectedPoints.size(); ++pointIdx) {
      CollectedPoint& cp = collectedPoints[pointIdx];
      cp.intensity *= normalizationValue;
    }
  }
}

float RemissionCalibrationHelper::CollectedPoints::getAverageCalibratedIntensity(const RemissionCalibrationResult& calibration) const {
  float ret = 0.0f;
  for (size_t pointIdx=0; pointIdx<size(); ++pointIdx) {
    const CollectedPoint& cp = at(pointIdx);
    ret += calibration.getCalibratedIntensity(cp.range, cp.incidenceAngle, cp.intensity, cp.laser);
  }
  ret /= size();
  return ret;
}

pcl::PointXYZI RemissionCalibrationHelper::CollectedPoints::getPointRepresentingCell() const {
  pcl::PointXYZI p;
  p.intensity = cellValue;
  p.x=p.y=p.z = 0;
  for (size_t pointIdx=0; pointIdx<size(); ++pointIdx) {
    const CollectedPoint& cp = at(pointIdx);
    p.getVector3fMap() += Eigen::Vector3f(cp.x, cp.y, cp.z);
  }
  p.getVector3fMap() /= size();
  return p;
}

float RemissionCalibrationHelper::CollectedPoints::getChi2(const RemissionCalibrationResult& calibration) const {
  float ret = 0.0f;
  for (size_t pointIdx=0; pointIdx<size(); ++pointIdx) {
    const CollectedPoint& cp = at(pointIdx);
    ret += std::pow(cellValue-calibration.getCalibratedIntensity(cp.range, cp.incidenceAngle, cp.intensity, cp.laser), 2);
  }
  return ret;
}

float RemissionCalibrationHelper::CollectedPoints::getRMSE(const RemissionCalibrationResult& calibration) const {
  if (empty()) {
    //std::cout << "Cell empty...\n";
    return 0.0f;
  }
  return sqrtf(getChi2(calibration)/size());
}

float RemissionCalibrationHelper::CollectedPointsMap::getChi2(const RemissionCalibrationResult& calibration) const {
  float ret = 0.0f;
  for (size_t mapCellIdx=0; mapCellIdx<size(); ++mapCellIdx) {
    ret += at(mapCellIdx).getChi2(calibration);
  }
  return ret;
}

size_t RemissionCalibrationHelper::CollectedPointsMap::getNoOfDataPoints() const {
  size_t ret = 0;
  for (size_t mapCellIdx=0; mapCellIdx<size(); ++mapCellIdx)
    ret += at(mapCellIdx).size();
  return ret;
}

float RemissionCalibrationHelper::CollectedPointsMap::getRMSE(const RemissionCalibrationResult& calibration) const {
  return sqrtf(getChi2(calibration)/getNoOfDataPoints());
}

void RemissionCalibrationHelper::CollectedPointsMap::setCellValuesToCalibratedMeans(const RemissionCalibrationResult& calibration) {
  for (size_t mapCellIdx=0; mapCellIdx<size(); ++mapCellIdx) {
    at(mapCellIdx).cellValue = at(mapCellIdx).getAverageCalibratedIntensity(calibration);
  }
}

void RemissionCalibrationHelper::CollectedPoint::projectOntoCellMiddle(const RemissionCalibrationResult& calibration) {
  float normalizedIntensity = calibration.getCalibratedIntensity(range, incidenceAngle, intensity, laser);
  int cellX, cellY, cellZ;
  calibration.rangeAndIncidenceAngleToCell(range, incidenceAngle, cellX, cellY);
  calibration.intensityToCell(intensity, laser, cellZ);
  calibration.cellToRange(cellX, range);
  calibration.cellToIncidenceAngle(cellY, incidenceAngle);
  //calibration.cellToIntensity(cellZ, laser, intensity);
  float normalizationFactor = calibration.getCalibrationFactor(cellX, cellY, cellZ, laser);
  intensity = normalizedIntensity / normalizationFactor;
}

void RemissionCalibrationHelper::CollectedPointsMap::projectOntoCellMiddles(const RemissionCalibrationResult& calibration) {
  for (size_t mapCellIdx=0; mapCellIdx<size(); ++mapCellIdx) {
    CollectedPoints& collectedPoints = at(mapCellIdx);
    for (size_t pointIdx=0; pointIdx<collectedPoints.size(); ++pointIdx)
      collectedPoints[pointIdx].projectOntoCellMiddle(calibration);;
  }
}

std::map<int,int> RemissionCalibrationHelper::CollectedPointsMap::getErrorDistribution(const RemissionCalibrationResult& calibration,
                                                                                       float stepSize) const
{
  std::map<int,int> ret;
  for (size_t mapCellIdx=0; mapCellIdx<size(); ++mapCellIdx) {
    const CollectedPoints& collectedPoints = at(mapCellIdx);
    for (size_t pointIdx=0; pointIdx<collectedPoints.size(); ++pointIdx) {
      const CollectedPoint& cp = collectedPoints[pointIdx];
      float error = std::abs(calibration.getCalibratedIntensity(cp.range, cp.incidenceAngle, cp.intensity, cp.laser)-collectedPoints.cellValue);
      int histogramPosition = lrintf(error/stepSize);
      if (ret.find(histogramPosition)==ret.end())
        ret[histogramPosition] = 1;
      else
        ++ret[histogramPosition];
    }
  }
  return ret;
}
