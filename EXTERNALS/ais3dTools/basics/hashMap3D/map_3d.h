#ifndef HASH_MAP_3D_H
#define HASH_MAP_3D_H

#include <Eigen/Core>
#include <tr1/unordered_map>

namespace Ais3dTools{

  template <typename CellType>
  //struct HashMap3D : public std::tr1::unordered_map<int, std::tr1::unordered_map<int, std::tr1::unordered_map<int, CellType > > >
  struct HashMap3D : public std::vector< std::tr1::unordered_map<int, CellType > >
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    HashMap3D(float gridSize, Eigen::Vector3i& size): _gridSize(gridSize), _origin(Eigen::Vector3f(0.0f, 0.0f, 0.0f)), _size(size), _center(size/2){
      allocate();

    };

    HashMap3D(float gridSize, Eigen::Vector3f& origin, Eigen::Vector3i& size): _gridSize(gridSize), _origin(origin), _size(size), _center(size/2){
      allocate();
    };

    void allocate(){
      this->resize(_size(0)*_size(1));
    }

    /** has to be called bevore adding first point */
    void setOrigin(Eigen::Vector3f& origin){ _origin = origin;}

    Eigen::Vector3f gridToWorld(Eigen::Vector3i index){
      index -= _center;
      Eigen::Vector3f point;
      point = _gridSize * index.cast<float>();
      point += _origin;
      return point;
    };

    Eigen::Vector3i worldToGrid(Eigen::Vector3f point){
      point = point - _origin;
      Eigen::Vector3i index;
      index(0) = floor(point(0) / _gridSize);
      index(1) = floor(point(1) / _gridSize);
      index(2) = floor(point(2) / _gridSize);
      index += _center;
      return index;
    };

    bool isInside(Eigen::Vector3i& index){
      if(index(0) > 0 && index(0) < _size(0) && index(1) > 0 && index(1) < _size(1))
          return true;
      return false;
    }

    CellType& cell(Eigen::Vector3i& index){
      return (*this)[index(0) + (_size(0) * index(1))][index(2)];
    }

    protected:
      float           _gridSize;
      Eigen::Vector3f _origin;
      Eigen::Vector3i _size;
      Eigen::Vector3i _center;

  };

} //end namespace
#endif