#ifndef TRANSFORMATION_REPRESENTATION_H
#define TRANSFORMATION_REPRESENTATION_H

#include <cmath>
#include <Eigen/Geometry>

namespace Ais3dTools {

namespace TransformationRepresentation {

template <typename real>
inline real getRollFromQuaternion(real w, real x, real y, real z) {
  return atan2(2.0 * (y*z + x*w), w*w - x*x - y*y + z*z);
}

template <typename real>
inline real getPitchFromQuaternion(real w, real x, real y, real z) {
  real tmp1 = -2.0 * (x*z - y*w);
  real tmp2 =  2.0 * (y*z + x*w);
  real tmp3 = w*w - x*x - y*y + z*z;
  return std::atan2(tmp1, sqrt(tmp2*tmp2 + tmp3*tmp3));
}

template <typename real>
inline real getYawFromQuaternion(real w, real x, real y, real z) {
  return std::atan2(2.0*(w*z + x*y), (1.0 - 2.0*(y*y + z*z)));
}

template <typename real>
inline void getEulerFromQuaternion(real w, real x, real y, real z, real &roll, real &pitch, real &yaw) {
  roll  = getRollFromQuaternion(w, x, y, z);
  pitch = getPitchFromQuaternion(w, x, y, z);
  yaw   = getYawFromQuaternion(w, x, y, z);
}

template <typename real>
inline void getQuaternionFromEuler(real roll, real pitch, real yaw, real& w, real& x, real& y, real& z) {
  real cosAlphaDiv2=std::cos(0.5f*yaw), sinAlphaDiv2=std::sin(0.5*yaw),  cosBetaDiv2=std::cos(0.5*pitch),
       sinBetaDiv2=std::sin(0.5*pitch), cosGammaDiv2=std::cos(0.5*roll), sinGammaDiv2=std::sin(0.5*roll);
  w = cosAlphaDiv2*cosBetaDiv2*cosGammaDiv2 + sinAlphaDiv2*sinBetaDiv2*sinGammaDiv2;
  x = cosAlphaDiv2*cosBetaDiv2*sinGammaDiv2 - sinAlphaDiv2*sinBetaDiv2*cosGammaDiv2;
  y = cosAlphaDiv2*sinBetaDiv2*cosGammaDiv2 + sinAlphaDiv2*cosBetaDiv2*sinGammaDiv2;
  z = sinAlphaDiv2*cosBetaDiv2*cosGammaDiv2 - cosAlphaDiv2*sinBetaDiv2*sinGammaDiv2;
  
  real normalizationFactor = 1.0 / std::sqrt(w*w + x*x + y*y + z*z);
  
  w *= normalizationFactor;
  x *= normalizationFactor;
  y *= normalizationFactor;
  z *= normalizationFactor;
}

template <typename real>
inline Eigen::Quaternion<real> getQuaternionFromMatrix(const Eigen::Matrix<real,3,3>& m) {
  real tmpQ1 = std::sqrt(fmax(0.0, 1.0 + m(0,0) - m(1,1) - m(2,2)))*0.5,
         tmpQ2 = std::sqrt(fmax(0.0, 1.0 - m(0,0) + m(1,1) - m(2,2)))*0.5,
         tmpQ3 = std::sqrt(fmax(0.0, 1.0 - m(0,0) - m(1,1) + m(2,2)))*0.5;
  
  Eigen::Quaternion<real> ret;
  ret.w() = std::sqrt(fmax(0.0, 1.0 + m(0,0) + m(1,1) + m(2,2)))*0.5,
  ret.x() = (m(2,1) - m(1,2))>=0.0?fabs(tmpQ1):-fabs(tmpQ1),
  ret.y() = (m(0,2) - m(2,0))>=0.0?fabs(tmpQ2):-fabs(tmpQ2),
  ret.z() = (m(1,0) - m(0,1))>=0.0?fabs(tmpQ3):-fabs(tmpQ3);
  return ret;
}

template <typename MatrixType, typename real>
inline void getMatrixFromEuler(real roll, real pitch, real yaw, MatrixType& m) {
  real cosAlpha=std::cos(yaw), sinAlpha=std::sin(yaw), cosBeta=std::cos(pitch),
       sinBeta=std::sin(pitch), cosGamma=std::cos(roll), sinGamma=std::sin(roll);
  m(0,0) = cosAlpha*cosBeta;
  m(0,1) = cosAlpha*sinBeta*sinGamma - sinAlpha*cosGamma;
  m(0,2) = cosAlpha*sinBeta*cosGamma + sinAlpha*sinGamma;
  m(1,0) = sinAlpha*cosBeta;
  m(1,1) = sinAlpha*sinBeta*sinGamma + cosAlpha*cosGamma;
  m(1,2) = sinAlpha*sinBeta*cosGamma - cosAlpha*sinGamma;
  m(2,0) = -sinBeta;
  m(2,1) = cosBeta*sinGamma;
  m(2,2) = cosBeta*cosGamma;
}

template <typename real>
inline Eigen::Matrix<real,3,3> getMatrixFromEuler(real roll, real pitch, real yaw) {
  Eigen::Matrix<real,3,3> ret;
  getMatrixFromEuler<Eigen::Matrix<real,3,3>, real>(roll, pitch, yaw, ret);
  return ret;
}

template <typename MatrixType, typename real>
inline void getMatrixFromTranslationAndEuler(real x, real y, real z, real roll, real pitch, real yaw, MatrixType& m)
{
  getMatrixFromEuler(roll, pitch, yaw, m);
  m(0,3) = x;
  m(1,3) = y;
  m(2,3) = z;
  m(3,0)=m(3,1)=m(3,2)=0.0f, m(3,3)=1.0f;
}


template <typename real>
inline Eigen::Matrix<real,4,4> getMatrixFromTranslationAndEuler(real x, real y, real z, real roll, real pitch, real yaw) {
  Eigen::Matrix<real,4,4> ret;
  getMatrixFromTranslationAndEuler<Eigen::Matrix<real,4,4>, real>(x, y, z, roll, pitch, yaw, ret);
  return ret;
}


template <typename real>
inline Eigen::Matrix<real,3,3> getMatrixFromQuaternion(const Eigen::Quaternion<real> q) {
  real tmp1 = 2.0*q.x()*q.x(), tmp2 = 2.0*q.y()*q.y(), tmp3 = 2.0*q.z()*q.z();
  Eigen::Matrix<real,3,3> ret;
  ret(0,0) = 1.0 - tmp2 - tmp3;
  ret(0,1) = 2.0 * (q.x()*q.y() - q.w()*q.z());
  ret(0,2) = 2.0 * (q.w()*q.y() + q.x()*q.z());
  ret(1,0) = 2.0 * (q.w()*q.z() + q.x()*q.y());
  ret(1,1) = 1.0 - tmp1 - tmp3;
  ret(1,2) = 2.0 * (q.y()*q.z() - q.w()*q.x());
  ret(2,0) = 2.0 * (q.x()*q.z() - q.w()*q.y());
  ret(2,1) = 2.0 * (q.w()*q.x() + q.y()*q.z());
  ret(2,2) = 1.0 - tmp1 - tmp2;
  return ret;
}

template <typename real>
inline void getEulerFromMatrix(const Eigen::Matrix<real,3,3> &m, real &roll, real &pitch, real &yaw) {
   yaw   = std::atan2(m(1,0), m(0,0));
   pitch = std::atan2(-m(2,0), sqrt(m(2,1)*m(2,1) + m(2,2)*m(2,2)));
   roll  = std::atan2(m(2,1), m(2,2));
}

template <typename real>
inline real getYawFromMatrix(const Eigen::Matrix<real,3,3> &m) {
   return std::atan2(m(1,0), m(0,0));
}
template <typename real>
inline real getPitchFromMatrix(const Eigen::Matrix<real,3,3> &m) {
   return std::atan2(-m(2,0), std::sqrt(m(2,1)*m(2,1) + m(2,2)*m(2,2)));
}
template <typename real>
inline real getRollFromMatrix(const Eigen::Matrix<real,3,3> &m) {
   return std::atan2(m(2,1), m(2,2));
}

template <typename MatrixType, typename real>
inline void getEulerAngles(const MatrixType& t, real& roll, real& pitch, real& yaw)
{
  roll  = atan2(t(2,1), t(2,2));
  pitch = asin(-t(2,0));
  yaw   = atan2(t(1,0), t(0,0));
}

template <typename MatrixType, typename real>
inline void getTranslationAndEulerAngles(const MatrixType& t, real& x, real& y, real& z, real& roll, real& pitch, real& yaw)
{
  x = t(0,3);
  y = t(1,3);
  z = t(2,3);
  getEulerAngles(t, roll, pitch, yaw);
}

//inline void getTransFromUnitVectorsZY(const Eigen::Vector3f& z_axis, const Eigen::Vector3f& y_direction, Eigen::Affine3f& transformation)
//{
  //Eigen::Vector3f tmp0 = (y_direction.cross(z_axis)).normalized();
  //Eigen::Vector3f tmp1 = (z_axis.cross(tmp0)).normalized();
  //Eigen::Vector3f tmp2 = z_axis.normalized();
  
  //transformation(0,0)=tmp0[0]; transformation(0,1)=tmp0[1]; transformation(0,2)=tmp0[2]; transformation(0,3)=0.0f;
  //transformation(1,0)=tmp1[0]; transformation(1,1)=tmp1[1]; transformation(1,2)=tmp1[2]; transformation(1,3)=0.0f;
  //transformation(2,0)=tmp2[0]; transformation(2,1)=tmp2[1]; transformation(2,2)=tmp2[2]; transformation(2,3)=0.0f;
  //transformation(3,0)=0.0f;    transformation(3,1)=0.0f;    transformation(3,2)=0.0f;    transformation(3,3)=1.0f;
//}

//! Get the unique 3D rotation that will rotate x_axis into (1,0,0), y_axis into(0,1,0) and z_axis into (0,0,1)
template <typename real1, typename real2>
inline void getRotationFromBasis(const Eigen::Matrix<real1, 3, 1>& x_axis,
                                 const Eigen::Matrix<real1, 3, 1>& y_axis,
                                 const Eigen::Matrix<real1, 3, 1>& z_axis,
                                 Eigen::Matrix<real2, 4, 4>& transformation)
{
  transformation(0,0)=x_axis[0]; transformation(0,1)=x_axis[1]; transformation(0,2)=x_axis[2]; transformation(0,3)=0.0f;
  transformation(1,0)=y_axis[0]; transformation(1,1)=y_axis[1]; transformation(1,2)=y_axis[2]; transformation(1,3)=0.0f;
  transformation(2,0)=z_axis[0]; transformation(2,1)=z_axis[1]; transformation(2,2)=z_axis[2]; transformation(2,3)=0.0f;
  transformation(3,0)=0.0f;      transformation(3,1)=0.0f;      transformation(3,2)=0.0f;      transformation(3,3)=1.0f;
}

//! Get the unique 3D rotation that will rotate x_axis into (1,0,0) and y_axis into(0,1,0)
template <typename real1, typename real2>
inline void getRotationFromBasisXY(const Eigen::Matrix<real1, 3, 1>& x_axis,
                                   const Eigen::Matrix<real1, 3, 1>& y_axis,
                                   Eigen::Matrix<real2, 4, 4>& transformation)
{
  Eigen::Vector3f z_axis = x_axis.cross(y_axis);
  getRotationFromBasis(x_axis, y_axis, z_axis, transformation);
}

//! Get the unique 3D rotation that will rotate x_axis into (1,0,0) and z_axis into (0,0,1)
template <typename real1, typename real2>
inline void getRotationFromBasisXZ(const Eigen::Matrix<real1, 3, 1>& x_axis,
                                   const Eigen::Matrix<real1, 3, 1>& z_axis,
                                   Eigen::Matrix<real2, 4, 4>& transformation)
{
  Eigen::Vector3f y_axis = z_axis.cross(x_axis);
  getRotationFromBasis(x_axis, y_axis, z_axis, transformation);
}

//! Get the unique 3D rotation that will rotate y_axis into(0,1,0) and z_axis into (0,0,1)
template <typename real1, typename real2>
inline void getRotationFromBasisYZ(const Eigen::Matrix<real1, 3, 1>& y_axis,
                                   const Eigen::Matrix<real1, 3, 1>& z_axis,
                                   Eigen::Matrix<real2, 4, 4>& transformation)
{
  Eigen::Vector3f x_axis = y_axis.cross(z_axis);
  getRotationFromBasis(x_axis, y_axis, z_axis, transformation);
}

/** Get the unique 3D rotation that will rotate x_direction into (?,0,0) and y_direction into(?,?,0).
  * The input vectors do not have to be normalized or orthogonal. */
template <typename real1, typename real2>
inline void getRotationFromXY(const Eigen::Matrix<real1, 3, 1>& x_direction,
                              const Eigen::Matrix<real1, 3, 1>& y_direction,
                              Eigen::Matrix<real2, 4, 4>& transformation)
{
  Eigen::Vector3f x_axis = x_direction.normalized(),
                  z_axis = x_direction.cross(y_direction).normalized();
  getRotationFromBasisXZ(x_axis, z_axis, transformation);
}

/** Get the unique 3D rotation that will rotate x_direction into (?,0,0) and z_direction into(?,0,?).
  * The input vectors do not have to be normalized or orthogonal. */
template <typename real1, typename real2>
inline void getRotationFromXZ(const Eigen::Matrix<real1, 3, 1>& x_direction,
                              const Eigen::Matrix<real1, 3, 1>& z_direction,
                              Eigen::Matrix<real2, 4, 4>& transformation)
{
  Eigen::Vector3f x_axis = x_direction.normalized(),
                  y_axis = z_direction.cross(x_direction).normalized();
  getRotationFromBasisXY(x_axis, y_axis, transformation);
}

/** Get the unique 3D rotation that will rotate y_direction into (0,?,0) and z_direction into(0,?,?).
  * The input vectors do not have to be normalized or orthogonal. */
template <typename real1, typename real2>
inline void getRotationFromYZ(const Eigen::Matrix<real1, 3, 1>& y_direction,
                              const Eigen::Matrix<real1, 3, 1>& z_direction,
                              Eigen::Matrix<real2, 4, 4>& transformation)
{
  Eigen::Vector3f y_axis = y_direction.normalized(),
                  x_axis = y_direction.cross(z_direction).normalized();
  getRotationFromBasisXY(x_axis, y_axis, transformation);
}

/** Get transformation that transforms normal1 into normal2. rotation defines the rotation (in radians) around the normal. */
template <typename real>
inline void getRotationBetweenTwoNormals(const Eigen::Matrix<real, 3, 1>& normal1,
                                         const Eigen::Matrix<real, 3, 1>& normal2,
                                         real rotation,
                                         Eigen::Matrix<real, 4, 4>& transformation)
{
  Eigen::Transform<real,3,Eigen::Isometry> rotTrans =
    Eigen::Transform<real,3,Eigen::Isometry>(Eigen::AngleAxis<real>(rotation, Eigen::Matrix<real, 3, 1>(1.0f, 0.0f, 0.0f)));
  if (normal1==normal2) {
    transformation = rotTrans*Eigen::Matrix<real, 4, 4>::Identity();
  }
  else {
    Eigen::Transform<real,3,Eigen::Isometry> transformation1, transformation2;
    getRotationFromXY(normal1, normal2, transformation1.matrix());
    transformation1 = rotTrans*transformation1;
    getRotationFromXY(normal2, normal1, transformation2.matrix());
    transformation = (transformation2.inverse()*transformation1).matrix();
  }
}


} // Namespace end
} // Namespace end

#endif
