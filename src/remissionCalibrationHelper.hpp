// This software is licenced under the Creative Commons
// (Attribution-NonCommercial-ShareAlike):
// http://creativecommons.org/licenses/by-nc-sa/4.0

// Author: Bastian Steder <steder@informatik.uni-freiburg.de>


float RemissionCalibrationHelper::getCosIncidenceAngle(const PointType& point, const Eigen::Vector3f& normal) {
  //return 1.0f;  // HACK!

  Eigen::Vector3f p = point.getVector3fMap(),
                  sensorPos(point.vp_x, point.vp_y, point.vp_z);
  Eigen::Vector3f viewingDirection = (sensorPos-p).normalized();
  float cosIncidenceAngleToUse = std::abs(viewingDirection.dot(normal));
  //std::cout << PVARN(cosIncidenceAngleToUse);
  cosIncidenceAngleToUse = std::min(1.0f, cosIncidenceAngleToUse);
  
  return cosIncidenceAngleToUse;
}

