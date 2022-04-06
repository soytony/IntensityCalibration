#ifndef AIS3TOOLS_SPARSE_POINT_CLOUD_H
#define AIS3TOOLS_SPARSE_POINT_CLOUD_H

#include <set>
#include <map>

namespace Ais3dTools {

struct SparsePointCloud {
  //-----STRUCTS&TYPEDEFS-----
  typedef std::set<int>               ContainerZ;
  typedef std::map<int, ContainerZ>   ContainerYZ;
  typedef std::map<int, ContainerYZ>  ContainerXYZ;
  typedef ContainerXYZ::iterator  IteratorX;
  typedef ContainerYZ::iterator   IteratorY;
  typedef ContainerZ::iterator    IteratorZ;
  typedef ContainerXYZ::const_iterator  ConstIteratorX;
  typedef ContainerYZ::const_iterator   ConstIteratorY;
  typedef ContainerZ::const_iterator    ConstIteratorZ;

  struct const_iterator {
    void initialize(const SparsePointCloud* newCloud) {
      cloud = newCloud;
      bool done = false;
      //for (itX=cloud->points.begin(); itX!=cloud.points.end()&&!done; ++itX)
        //for (itY=itX->second.begin(); itY!=itX->second.end()&&!done; ++itY)
          //for (itX=itX->second.begin(); itX!=itY->second.end()&&!done; ++itZ)
            //done = true;
      for (itX=cloud->points.begin(); itX!=cloud->points.end()&&!done;) {
        for (itY=itX->second.begin(); itY!=itX->second.end()&&!done;) {
          itZ = itY->second.begin();
          if (itZ!=itY->second.end()) {
            done = true;
            break;
          }
          ++itY;
        }
        if (!done)
          ++itX;
      }
      //if (!atEnd())
        //std::cout << x()<<";"<<y()<<";"<<z()<<"\n";
    }
    const_iterator& operator++() {
      if (atEnd())
        return *this;
      ++itZ;
      while (itZ==itY->second.end()) {
        ++itY;
        while (itY==itX->second.end()) {
          ++itX;
          if (atEnd())
            return *this;
          itY = itX->second.begin();
        }
        itZ = itY->second.begin();
      }
      return *this;
    }
    inline bool atEnd() const { return itX==cloud->points.end(); }
    inline float x() const { return float(itX->first)*cloud->cellSize; }
    inline float y() const { return float(itY->first)*cloud->cellSize; }
    inline float z() const { return float(*itZ)*cloud->cellSize; }
    //reference operator*() const;
    
    const SparsePointCloud* cloud;
    ConstIteratorX itX;
    ConstIteratorY itY;
    ConstIteratorZ itZ;
  };
  
  //-----CONSTRUCTOR-----
  SparsePointCloud(float cellSize=0.01f) : cellSize(cellSize) {}
  
  //-----MEMBER FUNCTIONS-----
  bool add(float x, float y, float z) {
    return points[lrintf(x/cellSize)][lrintf(y/cellSize)].insert(lrintf(z/cellSize)).second;
  }
  void getIndicesFromCoordinates(float x, float y, float z, int& idxX, int& idxY, int& idxZ) const {
    idxX=lrintf(x/cellSize), idxY=lrintf(y/cellSize), idxZ=lrintf(z/cellSize);
  }
  void getCoordinatesFromIndices(int idxX, int idxY, int idxZ, float& x, float& y, float& z) const {
    x = float(idxX)*cellSize,  y = float(idxY)*cellSize,  z = float(idxZ)*cellSize;
  }
  bool checkExistence(float x, float y, float z) const {
    ConstIteratorX itX = points.find(lrintf(x/cellSize));
    if (itX==points.end())  return false;
    ConstIteratorY itY = itX->second.find(lrintf(y/cellSize));
    if (itY==itX->second.end())  return false;
    ConstIteratorZ itZ = itY->second.find(lrintf(z/cellSize));
    if (itZ==itY->second.end())  return false;
    return true;
  }
  void clear() { points.clear(); }
  bool empty() const { return points.empty(); }
  const_iterator begin() const { const_iterator it; it.initialize(this); return it; }
  
  //-----MEMBER VARIABLES-----
  float cellSize;
  ContainerXYZ points;
};

template <typename PointType_>
struct SparsePointCloudT {
  //-----STRUCTS&TYPEDEFS-----
  typedef PointType_ PointType;
  typedef std::map<int, PointType_>   ContainerZ;
  typedef std::map<int, ContainerZ>   ContainerYZ;
  typedef std::map<int, ContainerYZ>  ContainerXYZ;
  typedef typename ContainerXYZ::iterator  IteratorX;
  typedef typename ContainerYZ::iterator   IteratorY;
  typedef typename ContainerZ::iterator    IteratorZ;
  typedef typename ContainerXYZ::const_iterator  ConstIteratorX;
  typedef typename ContainerYZ::const_iterator   ConstIteratorY;
  typedef typename ContainerZ::const_iterator    ConstIteratorZ;
  
  //-----CONSTRUCTOR-----
  SparsePointCloudT(float cellSize=0.01f) : cellSize(cellSize) {}
  
  //-----MEMBER FUNCTIONS-----
  void add(float x, float y, float z, const PointType_& point) {
    points[lrintf(x/cellSize)][lrintf(y/cellSize)][lrintf(z/cellSize)] = point;
  }
  void getIndicesFromCoordinates(float x, float y, float z, int& idxX, int& idxY, int& idxZ) const {
    idxX=lrintf(x/cellSize), idxY=lrintf(y/cellSize), idxZ=lrintf(z/cellSize);
  }
  void getCoordinatesFromIndices(int idxX, int idxY, int idxZ, float& x, float& y, float& z) const {
    x = float(idxX)*cellSize,  y = float(idxY)*cellSize,  z = float(idxZ)*cellSize;
  }
  void clear() { points.clear(); }
  bool empty() const { return points.empty(); }
  
  
  
  //-----MEMBER VARIABLES-----
  float cellSize;
  ContainerXYZ points;
};

}  // Namespace end

#endif
