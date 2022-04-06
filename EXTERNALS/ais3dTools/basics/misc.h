#ifndef STUFF_MISC_H
#define STUFF_MISC_H

#include <cmath>

/** @addtogroup utils **/
// @{

/** \file misc.h
 * \brief some general case utility functions
 *
 *  This file specifies some general case utility functions
 **/

namespace Ais3dTools {

/**
 * return the square value
 */
template <typename T>
inline T square(T x)
{
  return x*x;
}

/**
 * return the hypot of x and y
 */
template <typename T>
inline T hypot(T x, T y)
{
  return (T) (sqrt(x*x + y*y));
}

/**
 * return the squared hypot of x and y
 */
template <typename T>
inline T hypot_sqr(T x, T y)
{
  return x*x + y*y;
}

/**
 * convert from degree to radian
 */
inline double deg2rad(double degree)
{
  return degree * 0.01745329251994329576;
}

/**
 * convert from radian to degree
 */
inline double rad2deg(double rad)
{
  return rad * 57.29577951308232087721;
}

/**
 * normalize the angle
 */
template <typename real>
inline real normalize_theta(real theta)
{
  real multiplier;
  
  if (theta >= -M_PI && theta < M_PI)
    return theta;
  
  multiplier = floor(theta / (2*M_PI));
  theta = theta - multiplier*2*M_PI;
  if (theta >= M_PI)
    theta -= 2*M_PI;
  if (theta < -M_PI)
    theta += 2*M_PI;

  return theta;
}

/**
 * inverse of an angle, i.e., +180 degree
 */
inline double inverse_theta(double th)
{
  return normalize_theta(th + M_PI);
}

/**
 * average two angles
 */
inline double average_angle(double theta1, double theta2)
{
  double x, y;

  x = cos(theta1) + cos(theta2);
  y = sin(theta1) + sin(theta2);
  if(x == 0 && y == 0)
    return 0;
  else
    return std::atan2(y, x);
}

/**
 * sign function.
 * @return the sign of x. +1 for x > 0, -1 for x < 0, 0 for x == 0
 */
template <typename T>
inline int sign(T x)
{
  if (x > 0)
    return 1;
  else if (x < 0)
    return -1;
  else
    return 0;
}

/**
 * clamp x to the interval [l, u]
 */
template <typename T>
inline T clamp(T l, T x, T u) 
{
  if (x < l)
    return l;
  if (x > u)
    return u;
  return x;
}

/**
 * wrap x to be in the interval [l, u]
 */
template <typename T>
inline T wrap(T l, T x, T u) 
{
  T intervalWidth = u - l;
  while (x < l)
    x += intervalWidth;
  while (x > u)
    x -= intervalWidth;
  return x;
}

} // end namespace

// @}

#endif
