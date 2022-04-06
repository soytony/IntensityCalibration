// This software is licenced under the Creative Commons
// (Attribution-NonCommercial-ShareAlike):
// http://creativecommons.org/licenses/by-nc-sa/4.0

// Author: Bastian Steder <steder@informatik.uni-freiburg.de>


#include "g2oGraphElements.h"

//bool numericalJacobian = true;
//bool numericalJacobian = false;

RemissionLearnerVertex::RemissionLearnerVertex() :
  BaseClass()
{
}

bool RemissionLearnerVertex::read(std::istream& is)
{
  std::cerr << __PRETTY_FUNCTION__<<": Not implemented. :-(\n";
  return true;
}

bool RemissionLearnerVertex::write(std::ostream& os) const
{
  std::cerr << __PRETTY_FUNCTION__<<": Not implemented. :-(\n";
  return os.good();
}

RemissionLearnerUnaryEdge::RemissionLearnerUnaryEdge() : BaseClass() {
  _measurement = 0.0;
  information() = InformationType::Identity();
}

bool RemissionLearnerUnaryEdge::read(std::istream& is) {
  std::cerr << __PRETTY_FUNCTION__<<": Not implemented. :-(\n";
  return true;
}

bool RemissionLearnerUnaryEdge::write(std::ostream& os) const {
  std::cerr << __PRETTY_FUNCTION__<<": Not implemented. :-(\n";
  return os.good();
}
 
void RemissionLearnerUnaryEdge::linearizeOplus() {
  //if (numericalJacobian)
    //linearizeOplusNumeric();
  //else
    linearizeOplusJacobian();
  
  //double& dx = _jacobianOplusXi(0,0);
  
  //std::cout << "--\n";
  //linearizeOplusJacobian();
  //std::cout << "Jacobian: "<<dx;
  
  //linearizeOplusNumeric();
  //std::cout << ", Numeric: "<<dx<<"\n";
}

RemissionLearnerPriorEdge::RemissionLearnerPriorEdge(float value, float weight) : BaseClass() {
  _measurement = value;
  information()(0,0) = weight;
}

bool RemissionLearnerPriorEdge::read(std::istream& is) {
  std::cerr << __PRETTY_FUNCTION__<<": Not implemented. :-(\n";
  return true;
}

bool RemissionLearnerPriorEdge::write(std::ostream& os) const {
  std::cerr << __PRETTY_FUNCTION__<<": Not implemented. :-(\n";
  return os.good();
}

void RemissionLearnerPriorEdge::linearizeOplus() {
  //double x = static_cast<const RemissionLearnerVertex*>(_vertices[0])->estimate();
  double& dx = _jacobianOplusXi(0,0);
  dx =  1.0f;
}

RemissionLearnerMultiEdge_n::RemissionLearnerMultiEdge_n(int dimension) : BaseClass() {
  information() = InformationType::Identity();
  resize(dimension);
}

bool RemissionLearnerMultiEdge_n::read(std::istream& is) {
  std::cerr << __PRETTY_FUNCTION__<<": Not implemented. :-(\n";
  return true;
}

bool RemissionLearnerMultiEdge_n::write(std::ostream& os) const {
  std::cerr << __PRETTY_FUNCTION__<<": Not implemented. :-(\n";
  return os.good();
}
 
void RemissionLearnerMultiEdge_n::linearizeOplus() {
  //if (numericalJacobian)
    linearizeOplusNumeric();
  //else
    //linearizeOplusJacobian();
  
  //double& dx = _jacobianOplusXi(0,0);
  
  //std::cout << "--\n";
  //linearizeOplusJacobian();
  //std::cout << "Jacobian: "<<dx;
  
  //linearizeOplusNumeric();
  //std::cout << ", Numeric: "<<dx<<"\n";
}

RemissionCalibrationGridEdge::RemissionCalibrationGridEdge(RemissionLearnerVertex* v1, RemissionLearnerVertex* v2, double weight) :
    BaseClass()
{
  vertices()[0] = v1;
  vertices()[1] = v2;
  information()(0,0) = weight;
  _measurement = v1->estimate()-v2->estimate();
}

bool RemissionCalibrationGridEdge::read(std::istream& is) {
  std::cerr << __PRETTY_FUNCTION__<<": Not implemented. :-(\n";
  return true;
}

bool RemissionCalibrationGridEdge::write(std::ostream& os) const {
  std::cerr << __PRETTY_FUNCTION__<<": Not implemented. :-(\n";
  return os.good();
}

void RemissionCalibrationGridEdge::linearizeOplus() {
  //double x = static_cast<const RemissionLearnerVertex*>(_vertices[0])->estimate(),
         //y = static_cast<const RemissionLearnerVertex*>(_vertices[1])->estimate();
  //double& dx = _jacobianOplusXi(0,0),
        //& dy = _jacobianOplusXj(0,0);
  BaseClass::linearizeOplus();
  //double numericDx=dx, numericDy=dy;
  //dx =  2.0*y / std::pow(x+y,2);
  //dy = -2.0*x / std::pow(x+y,2);
  //std::cout << numericDx<<"="<<dx<<", "<<numericDy<<"="<<dy<<"\n";
}

RemissionCalibrationEdge_2::RemissionCalibrationEdge_2() :
  BaseClass()
{
  information() = InformationType::Identity();
}

bool RemissionCalibrationEdge_2::read(std::istream& is)
{
  std::cerr << __PRETTY_FUNCTION__<<": Not implemented. :-(\n";
  return true;
}

bool RemissionCalibrationEdge_2::write(std::ostream& os) const
{
  std::cerr << __PRETTY_FUNCTION__<<": Not implemented. :-(\n";
  return os.good();
}
 
void RemissionCalibrationEdge_2::linearizeOplus() {
  //if (numericalJacobian)
    //linearizeOplusNumeric();
  //else
    linearizeOplusJacobian();
  
  //double& dx = _jacobianOplusXi(0,0),
        //& dy = _jacobianOplusXj(0,0);
  
  //std::cout << "--\n";
  //linearizeOplusJacobian();
  //std::cout << "Jacobian: "<<dx<<", "<<dy<<"\n";
  
  //linearizeOplusNumeric();
  //std::cout << "Numeric: "<<dx<<", "<<dy<<"\n";
}

RemissionCalibrationEdgeMode2::RemissionCalibrationEdgeMode2() : BaseClass()
{
  information() = InformationType::Identity();
  resize(4);
}

bool RemissionCalibrationEdgeMode2::read(std::istream& is) {
  std::cerr << __PRETTY_FUNCTION__<<": Not implemented. :-(\n";
  return true;
}

bool RemissionCalibrationEdgeMode2::write(std::ostream& os) const {
  std::cerr << __PRETTY_FUNCTION__<<": Not implemented. :-(\n";
  return os.good();
}
 
void RemissionCalibrationEdgeMode2::linearizeOplus() {
  //if (numericalJacobian)
    //linearizeOplusNumeric();
  //else
    linearizeOplusJacobian();
}

RemissionCalibrationEdgeMode3::RemissionCalibrationEdgeMode3() : BaseClass() {
  information() = InformationType::Identity();
  resize(5);
}

bool RemissionCalibrationEdgeMode3::read(std::istream& is) {
  std::cerr << __PRETTY_FUNCTION__<<": Not implemented. :-(\n";
  return true;
}

bool RemissionCalibrationEdgeMode3::write(std::ostream& os) const {
  std::cerr << __PRETTY_FUNCTION__<<": Not implemented. :-(\n";
  return os.good();
}
 
void RemissionCalibrationEdgeMode3::linearizeOplus() {
  //if (numericalJacobian)
    //linearizeOplusNumeric();
  //else
    linearizeOplusJacobian();
}


g2o::HyperGraphAction* G2oPostIterationAction::operator()(const g2o::HyperGraph* graph, g2o::HyperGraphAction::Parameters* parameters) {
  (void) graph;
  (void) parameters;
  //int iteration = dynamic_cast<g2o::HyperGraphAction::ParametersIteration*>(parameters)->iteration;
  
  if (remissionCalibrationResult==NULL)
    return this;
  
  if (remissionCalibrationResult->mode==1) {
    for (VerticesMode1::const_iterator it=verticesMode1->begin(); it!=verticesMode1->end(); ++it) {
      int laser = it->first;
      const std::vector<RemissionLearnerVertex*>& vertices = it->second;
      std::vector<float>& calibrationValues = remissionCalibrationResult->calibrationValuesMode1[laser];
      for (size_t verticesIdx=0; verticesIdx<vertices.size(); ++verticesIdx)
        calibrationValues[verticesIdx] = vertices[verticesIdx]->estimate();
    }
  }
  else if (remissionCalibrationResult->mode==2) {
    for (VerticesMode2::const_iterator it=verticesMode2->begin(); it!=verticesMode2->end(); ++it) {
      int laser = it->first;
      const CalibIndependentVertices& v = it->second;
      RemissionCalibrationResult::CalibIndependent& c = remissionCalibrationResult->calibrationValuesMode2[laser];
      for (size_t verticesIdx=0; verticesIdx<v.verticesRange.size(); ++verticesIdx)
        c.calibrationValuesRange[verticesIdx] = v.verticesRange[verticesIdx]->estimate();
      for (size_t verticesIdx=0; verticesIdx<v.verticesIncidenceAngle.size(); ++verticesIdx)
        c.calibrationValuesIncidenceAngle[verticesIdx] = v.verticesIncidenceAngle[verticesIdx]->estimate();
      for (size_t verticesIdx=0; verticesIdx<v.verticesIntensity.size(); ++verticesIdx)
        c.calibrationValuesIntensity[verticesIdx] = v.verticesIntensity[verticesIdx]->estimate();
    }
  }
  else {
    for (size_t laserGroupIdx=0; laserGroupIdx<remissionCalibrationResult->calibrationValuesMode3.size(); ++laserGroupIdx) {
      RemissionCalibrationResult::LaserGroup& lg = remissionCalibrationResult->calibrationValuesMode3[laserGroupIdx];
      const LaserGroupVertices& lgv = verticesMode3->at(laserGroupIdx);
      for (size_t cellX=0; cellX<lg.calibrationValuesRange.size(); ++cellX)
        lg.calibrationValuesRange[cellX] = lgv.verticesRange[cellX]->estimate();
      for (size_t cellY=0; cellY<lg.calibrationValuesIncidenceAngle.size(); ++cellY)
        lg.calibrationValuesIncidenceAngle[cellY] = lgv.verticesIncidenceAngle[cellY]->estimate();
      for (size_t cellZ=0; cellZ<lg.calibrationValuesIntensity.size(); ++cellZ)
        lg.calibrationValuesIntensity[cellZ] = lgv.verticesIntensity[cellZ]->estimate();
      for (std::map<int, float>::const_iterator it=lg.calibrationValuesLaser.begin(); it!=lg.calibrationValuesLaser.end(); ++it)
        lg.calibrationValuesLaser[it->first] = lgv.verticesLaser.find(it->first)->second->estimate();
    }
  }
  
  //for (size_t mapCellIdx=0; mapCellIdx<mapCells->size(); ++mapCellIdx) {
    //mapCells->at(mapCellIdx).cellValue = verticesForMapCells->at(mapCellIdx)->estimate();
  //}
  
  //mapCells->setCellValuesToCalibratedMeans(*remissionCalibrationResult);
  
  //for (size_t mapCellIdx=0; mapCellIdx<mapCells->size(); ++mapCellIdx) {
    //verticesForMapCells->at(mapCellIdx)->setEstimate(mapCells->at(mapCellIdx).cellValue);
  //}
  
  //remissionCalibrationResult->normalizeCalibrationFactors(1.0f / mapCells->getAverageCalibratedIntensity(*remissionCalibrationResult));
  //float currentRMSE = mapCells->getRMSE(*remissionCalibrationResult);
  //std::cout << "RMSE of iteration "<<iteration<<" is "<<currentRMSE<<".\n\n";
  std::cout << "."<<std::flush;
  
  
  //if (remissionCalibrationResult->mode==1) {
    //for (std::map<int, std::vector<RemissionLearnerVertex*> >::const_iterator it=verticesMode1->begin(); it!=verticesMode1->end(); ++it) {
      //int laser = it->first;
      //const std::vector<RemissionLearnerVertex*>& vertices = it->second;
      //std::vector<float>& calibrationValues = remissionCalibrationResult->calibrationValuesMode1[laser];
      //for (size_t verticesIdx=0; verticesIdx<vertices.size(); ++verticesIdx)
        //vertices[verticesIdx]->setEstimate(calibrationValues[verticesIdx]);
    //}
  //}
  //else if (remissionCalibrationResult->mode==2) {
    //for (std::map<int, std::vector<RemissionLearnerVertex*> >::const_iterator it=verticesRangePerLaser->begin(); it!=verticesRangePerLaser->end(); ++it) {
      //int laser = it->first;
      //const std::vector<RemissionLearnerVertex*>& verticesRange = it->second;
      //const std::vector<RemissionLearnerVertex*>& verticesIncidenceAngle = verticesIncidenceAnglePerLaser->find(laser)->second;
      //std::vector<float>& calibrationValuesRange = remissionCalibrationResult->calibrationValuesRangePerLaser[laser];
      //std::vector<float>& calibrationValuesIncidenceAngle = remissionCalibrationResult->calibrationValuesIncidenceAnglePerLaser[laser];
      //for (size_t verticesIdx=0; verticesIdx<verticesRange.size(); ++verticesIdx)
        //verticesRange[verticesIdx]->setEstimate(calibrationValuesRange[verticesIdx]);
      //for (size_t verticesIdx=0; verticesIdx<verticesIncidenceAngle.size(); ++verticesIdx)
        //verticesIncidenceAngle[verticesIdx]->setEstimate(calibrationValuesIncidenceAngle[verticesIdx]);
    //}
  //}
  //else {
    //for (size_t laserGroupIdx=0; laserGroupIdx<remissionCalibrationResult->calibrationValuesMode3.size(); ++laserGroupIdx) {
      //RemissionCalibrationResult::LaserGroup& lg = remissionCalibrationResult->calibrationValuesMode3[laserGroupIdx];
      //const LaserGroupVertices& lgv = verticesMode3->at(laserGroupIdx);
      //for (size_t cellX=0; cellX<lg.calibrationValuesRange.size(); ++cellX)
        //lgv.verticesRange[cellX]->setEstimate(lg.calibrationValuesRange[cellX]);
      //for (size_t cellY=0; cellY<lg.calibrationValuesIncidenceAngle.size(); ++cellY)
        //lgv.verticesIncidenceAngle[cellY]->setEstimate(lg.calibrationValuesIncidenceAngle[cellY]);
      //for (std::map<int, float>::const_iterator it=lg.calibrationValuesLaser.begin(); it!=lg.calibrationValuesLaser.end(); ++it)
        //lgv.verticesLaser.find(it->first)->second->setEstimate(lg.calibrationValuesLaser[it->first]);
    //}
  //}
  
  //if (qApplication==NULL)
    //return this;
  
  //if (viewer != NULL)
    //viewer->updateGL();
  //qApplication->processEvents();
  
  return this;
}

