#ifndef AIS3DTOOLS_PCL_DOWNWARDS_COMPATIBILITY_H
#define AIS3DTOOLS_PCL_DOWNWARDS_COMPATIBILITY_H

#include "pcl/pcl_config.h"
#include "pcl/pcl_base.h"

#if PCL_VERSION_COMPARE(<,1,7,0)
  namespace pcl {
    typedef sensor_msgs::PointCloud2 PCLPointCloud2;
  }
  #define fromPCLPointCloud2 fromROSMsg
  #define toPCLPointCloud2 toROSMsg
#endif

#ifdef USE_ROS
  #define USING_PCL_FROM_ROS 1
#else

#endif

#endif
