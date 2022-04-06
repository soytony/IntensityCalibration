#ifndef HASH_MAP_3D_H
#define HASH_MAP_3D_H

#include <Eigen/Core>
#include <tr1/unordered_map>

namespace Ais3dTools{

  template <typename CellType>
  //struct HashMap3D : public std::tr1::unordered_map<int, std::tr1::unordered_map<int, std::tr1::unordered_map<int, CellType > > >
  struct HashMap3D : public std::tr1::unordered_map<int, std::tr1::unordered_map<int, std::tr1::unordered_map<int, CellType > > >
  {
    // ===== TYPEDEFS =====
    typedef std::tr1::unordered_map<int, CellType> MapZ;
    typedef std::tr1::unordered_map<int, MapZ>     MapY;
    typedef std::tr1::unordered_map<int, MapY>     MapX;
    typedef typename MapZ::iterator IteratorZ;
    typedef typename MapZ::const_iterator ConstIteratorZ;
    typedef typename MapY::iterator IteratorY;
    typedef typename MapY::const_iterator ConstIteratorY;
    typedef typename MapX::iterator IteratorX;
    typedef typename MapX::const_iterator ConstIteratorX;

    // ===== CONSTRUCTOR =====
    HashMap3D(float gridSize): _gridSize(gridSize), _worldToGrid(1.0f/_gridSize) {}

    // ===== PUBLIC MEMBER FUNCTIONS =====
    inline float getGridSize() const { return _gridSize; }
    
    inline Eigen::Vector3f gridToWorld(const Eigen::Vector3i& index) const {
      Eigen::Vector3f point;
      point = _gridSize * index.cast<float>();
      return point;
    }

    inline Eigen::Vector3i worldToGrid(const Eigen::Vector3f& point) const {
      Eigen::Vector3i index;
      index(0) = floor(_worldToGrid * point(0));
      index(1) = floor(_worldToGrid * point(1));
      index(2) = floor(_worldToGrid * point(2));
      return index;
    };
    
    inline CellType& cell(const Eigen::Vector3i& index) {
      return (*this)[index(0)][index(1)][index(2)];
    }
    
    inline CellType& cellFromPoint(const Eigen::Vector3f& point) {
      Eigen::Vector3i index = worldToGrid(point);
      return (*this)[index(0)][index(1)][index(2)];
    }

    inline bool exists(const Eigen::Vector3i& index) const {
      MapZ& mapZ = (*this)[index(0)][index(1)];
      ConstIteratorZ it = mapZ.find(index(2));
      return it!=mapZ.end();
    }
    
    inline CellType& findOrCreate(const Eigen::Vector3i& index, const CellType& defaultValue) {
      MapZ& mapZ = (*this)[index(0)][index(1)];
      IteratorZ it = mapZ.find(index(2));
      if (it == mapZ.end()) {  // Not existing yet?
        CellType& ret = mapZ[index(2)];
        ret = defaultValue;
        return ret;
      }
      else {
        return it->second;
      }
    }


    protected:
      // ===== PROTECTED MEMBER VARIABLES =====
      float _gridSize;
      float _worldToGrid;

  };

} //end namespace
#endif
