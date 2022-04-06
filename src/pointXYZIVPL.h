// This software is licenced under the Creative Commons
// (Attribution-NonCommercial-ShareAlike):
// http://creativecommons.org/licenses/by-nc-sa/4.0

// Author: Bastian Steder <steder@informatik.uni-freiburg.de>


#ifndef POINTXYZIVPL_H
#define POINTXYZIVPL_H

#include <pcl/point_types.h>

namespace pcl {

struct EIGEN_ALIGN16 _PointXYZIVPL {
  PCL_ADD_POINT4D; // This adds the members x,y,z which can also be accessed using the point (which is float[4])
  union
  {
    struct
    {
      float intensity;
      float vp_x;
      float vp_y;
      float vp_z;
    };
    float data2[4];
  };
  union
  {
    struct
    {
      int laser;
      int unused1;
      int unused2;
      int unused3;
    };
    float data3[4];
  };
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};
struct EIGEN_ALIGN16 PointXYZIVPL : public _PointXYZIVPL
{
  inline PointXYZIVPL (const _PointXYZIVPL &p)
  {
    x = p.x; y = p.y; z = p.z; data[3] = 1.0f;
    intensity = p.intensity;
    vp_x = p.vp_x; vp_y = p.vp_y; vp_z = p.vp_z;
    laser = p.laser;
  }
  
  inline PointXYZIVPL (float _x = 0.0f, float _y = 0.0f, float _z = 0.0f,
                      float _intensity = 0.0f, float _vp_x = 0.0f, float _vp_y = 0.0f, float _vp_z = 0.0f, int _laser=0)
  {
    x = _x; y = _y; z = _z;
    data[3] = 1.0f;
    intensity = _intensity;
    vp_x = _vp_x; vp_y = _vp_y; vp_z = _vp_z;
    laser = _laser;
    data3[1]=data3[2]=data3[3] = 1.0f;
  }
};
inline std::ostream& operator << (std::ostream& os, const PointXYZIVPL& p)
{
  os << "("   << p.x    << "," << p.y    << "," << p.z    << " - " << p.intensity
     << " - " << p.vp_x << "," << p.vp_y << "," << p.vp_z << " - " << p.laser << ")";
  return (os);
}

}

POINT_CLOUD_REGISTER_POINT_STRUCT (pcl::_PointXYZIVPL,
    (float, x, x)
    (float, y, y)
    (float, z, z)
    (float, intensity, intensity)
    (float, vp_x, vp_x)
    (float, vp_y, vp_y)
    (float, vp_z, vp_z)
    (int, laser, laser)
)
POINT_CLOUD_REGISTER_POINT_WRAPPER(pcl::PointXYZIVPL, pcl::_PointXYZIVPL)

#endif
