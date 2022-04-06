// This software is licenced under the Creative Commons
// (Attribution-NonCommercial-ShareAlike):
// http://creativecommons.org/licenses/by-nc-sa/4.0

// Author: Bastian Steder <steder@informatik.uni-freiburg.de>


void RemissionLearnerUnaryEdge::computeError() {
  double& e = _error[0];
  double v = static_cast<const RemissionLearnerVertex*>(_vertices[0])->estimate();
  e = 0.0;
  if (v < 1.0)
    e = (1.0/v) + v - 2.0;
}

void RemissionLearnerUnaryEdge::linearizeOplusNumeric() {
  BaseClass::linearizeOplus();
}

void RemissionLearnerUnaryEdge::linearizeOplusJacobian() {
  double& dx = _jacobianOplusXi(0,0);
  double v = static_cast<const RemissionLearnerVertex*>(_vertices[0])->estimate();
  dx = 0.0;
  if (v < 1.0)
    dx = (-1.0/(v*v)) + 1;
}

void RemissionLearnerPriorEdge::computeError() {
  double& e = _error[0];
  double v = static_cast<const RemissionLearnerVertex*>(_vertices[0])->estimate();
  e = v-_measurement;
}

void RemissionLearnerMultiEdge_n::computeError() {
  double averageNormalizationFactor = 0;
  for (size_t i=0; i<_vertices.size(); ++i) {
    averageNormalizationFactor += static_cast<const RemissionLearnerVertex*>(_vertices[i])->estimate();
    //std::cout << static_cast<const RemissionLearnerVertex*>(_vertices[i])->estimate()<<"+";
  }
  //std::cout << " = "<<averageNormalizationFactor<<"\n";
  averageNormalizationFactor /= double(_vertices.size());
  
  _error[0] = averageNormalizationFactor-1.0;
  //std::cout << PVARN(_error[0]);
}

void RemissionLearnerMultiEdge_n::linearizeOplusJacobian() {
  std::cerr << __PRETTY_FUNCTION__<<": Not implemented. :-(\n";
}


void RemissionLearnerMultiEdge_n::linearizeOplusNumeric() {
  BaseClass::linearizeOplus();
}

void RemissionCalibrationGridEdge::computeError() {
  double normalizationFactor1 = static_cast<const RemissionLearnerVertex*>(_vertices[0])->estimate(),
         normalizationFactor2 = static_cast<const RemissionLearnerVertex*>(_vertices[1])->estimate();
  //_error[0] = (normalizationFactor1-normalizationFactor2);
  _error[0] = ((normalizationFactor1-normalizationFactor2)-_measurement)/(normalizationFactor1+normalizationFactor2);
}

void RemissionCalibrationEdge_2::computeError() {
  double normalizationFactor = static_cast<const RemissionLearnerVertex*>(_vertices[0])->estimate(),
         cellIntensity       = static_cast<const RemissionLearnerVertex*>(_vertices[1])->estimate(),
         measuredIntensity   = _measurement;
  _error[0] = normalizationFactor*measuredIntensity - cellIntensity;
}

void RemissionCalibrationEdge_2::linearizeOplusNumeric() {
  BaseClass::linearizeOplus();
}

void RemissionCalibrationEdge_2::linearizeOplusJacobian() {
  //double normalizationFactor = static_cast<const RemissionLearnerVertex*>(_vertices[0])->estimate(),  // x
         //cellIntensity       = static_cast<const RemissionLearnerVertex*>(_vertices[1])->estimate();  // y
  double measuredIntensity = _measurement;
  double& dx = _jacobianOplusXi(0,0),
        & dy = _jacobianOplusXj(0,0);
  dx = measuredIntensity;
  dy = 0;
}

void RemissionCalibrationEdgeMode2::computeError() {
  double normalizationFactorRange     = static_cast<const RemissionLearnerVertex*>(_vertices[0])->estimate(),
         normalizationFactorAngle     = static_cast<const RemissionLearnerVertex*>(_vertices[1])->estimate(),
         normalizationFactorIntensity = static_cast<const RemissionLearnerVertex*>(_vertices[2])->estimate(),
         cellIntensity                = static_cast<const RemissionLearnerVertex*>(_vertices[3])->estimate(),
         measuredIntensity            = _measurement;
  _error[0] = normalizationFactorRange*normalizationFactorAngle*normalizationFactorIntensity*measuredIntensity - cellIntensity;
}

void RemissionCalibrationEdgeMode2::linearizeOplusJacobian() {
  double normalizationFactorRange     = static_cast<const RemissionLearnerVertex*>(_vertices[0])->estimate(), // x1
         normalizationFactorAngle     = static_cast<const RemissionLearnerVertex*>(_vertices[1])->estimate(), // x2
         normalizationFactorIntensity = static_cast<const RemissionLearnerVertex*>(_vertices[2])->estimate(), // x3
         //cellIntensity              = static_cast<const RemissionLearnerVertex*>(_vertices[3])->estimate(), // y
         measuredIntensity            = _measurement;
  double& dx1 = _jacobianOplus[0](0,0),
        & dx2 = _jacobianOplus[1](0,0),
        & dx3 = _jacobianOplus[2](0,0),
        & dy  = _jacobianOplus[3](0,0);
  dx1 = normalizationFactorAngle*normalizationFactorIntensity*measuredIntensity;
  dx2 = normalizationFactorRange*normalizationFactorIntensity*measuredIntensity;
  dx3 = normalizationFactorRange*normalizationFactorRange*measuredIntensity;
  dy  = 0;
}

void RemissionCalibrationEdgeMode2::linearizeOplusNumeric() {
  BaseClass::linearizeOplus();
}

void RemissionCalibrationEdgeMode3::computeError() {
  double normalizationFactorRange     = static_cast<const RemissionLearnerVertex*>(_vertices[0])->estimate(),
         normalizationFactorAngle     = static_cast<const RemissionLearnerVertex*>(_vertices[1])->estimate(),
         normalizationFactorIntensity = static_cast<const RemissionLearnerVertex*>(_vertices[2])->estimate(),
         normalizationFactorLaser     = static_cast<const RemissionLearnerVertex*>(_vertices[3])->estimate(),
         cellIntensity                = static_cast<const RemissionLearnerVertex*>(_vertices[4])->estimate(),
         measuredIntensity            = _measurement;
  _error[0] = normalizationFactorRange*normalizationFactorAngle*normalizationFactorIntensity*normalizationFactorLaser*measuredIntensity -
              cellIntensity;
}

void RemissionCalibrationEdgeMode3::linearizeOplusJacobian() {
  double normalizationFactorRange     = static_cast<const RemissionLearnerVertex*>(_vertices[0])->estimate(), // x1
         normalizationFactorAngle     = static_cast<const RemissionLearnerVertex*>(_vertices[1])->estimate(), // x2
         normalizationFactorIntensity = static_cast<const RemissionLearnerVertex*>(_vertices[2])->estimate(), // x3
         normalizationFactorLaser     = static_cast<const RemissionLearnerVertex*>(_vertices[3])->estimate(), // x4
         //cellIntensity              = static_cast<const RemissionLearnerVertex*>(_vertices[4])->estimate(), // y
         measuredIntensity            = _measurement;
  double& dx1 = _jacobianOplus[0](0,0),
        & dx2 = _jacobianOplus[1](0,0),
        & dx3 = _jacobianOplus[2](0,0),
        & dx4 = _jacobianOplus[3](0,0),
        & dy  = _jacobianOplus[4](0,0);
  dx1 = normalizationFactorAngle*normalizationFactorIntensity*normalizationFactorLaser*measuredIntensity;
  dx2 = normalizationFactorRange*normalizationFactorIntensity*normalizationFactorLaser*measuredIntensity;
  dx3 = normalizationFactorRange*normalizationFactorAngle*normalizationFactorLaser*measuredIntensity;
  dx4 = normalizationFactorRange*normalizationFactorAngle*normalizationFactorIntensity*measuredIntensity;
  dy  = 0;
  
  //std::cout << "-----\n";
  //std::cout << dx1<<", "<<dx2<<", "<<dx3<<", "<<dx4<<", "<<dy<<"\n";
  //linearizeOplusNumeric();
  //std::cout << dx1<<", "<<dx2<<", "<<dx3<<", "<<dx4<<", "<<dy<<"\n";
}

void RemissionCalibrationEdgeMode3::linearizeOplusNumeric() {
  BaseClass::linearizeOplus();
}

