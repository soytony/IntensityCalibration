#ifndef HASH_MAP_3D_SEARCH_H
#define HASH_MAP_3D_SEARCH_H

#include <algorithm>    
#include "hash_map_3d.h"
#include "pcl/point_types.h"
#include "pcl/point_cloud.h"


namespace Ais3dTools{

  struct HashMap3DSearchIndexComparator {
    public:
    
    HashMap3DSearchIndexComparator(std::vector<float>* distances) : distanceArray(distances) { ; }
    
    bool operator() (int i, int j) { return (i < j);}

    private:
      std::vector<float>* distanceArray;
  };
  
  
  template <typename PointType>
  struct HashMap3DSearchT : public HashMap3D<std::vector<int> >
  {

    HashMap3DSearchT(float gridSize): HashMap3D<std::vector<int> >(gridSize){
      
    };


    void setInput(typename pcl::PointCloud< PointType >::Ptr cloud)
    {
      pointcloud = cloud;
      for(size_t i = 0; i < cloud->points.size(); ++i){
        Eigen::Vector3i index = worldToGrid(cloud->points[i].getVector3fMap());
        cell(index).push_back(i);
      }
    }    

    void radiusSearch(PointType point, float sqrDistance, std::vector<int>& indices, std::vector<float>& distances){
      indices.clear();
      distances.clear();
      Eigen::Vector3f p = point.getVector3fMap();
      Eigen::Vector3i index = worldToGrid(p);
      int searchRange = lrint(sqrt(sqrDistance) / _gridSize) + 1;
      for(int x = -searchRange; x <= searchRange; ++x)
        for(int y = -searchRange; y <= searchRange; ++y)
          for(int z = -searchRange; z <= searchRange; ++z){
            Eigen::Vector3i cellIndex(x,y,z);
            cellIndex += index;
            
            for(size_t i = 0; i < cell(cellIndex).size(); ++i){
              float sqrDist = (p - (*pointcloud)[cell(cellIndex)[i]].getVector3fMap()).squaredNorm();
              if(sqrDist < sqrDistance){
                indices.push_back(cell(cellIndex)[i]);
                distances.push_back(sqrDist);
              }
            }
          }

              
      ///sort results
      std::vector<int> sorted_indices(indices.size());
      std::vector<float> sorted_distances(indices.size());
      
      std::vector<int> indexWrap(indices.size());
      for(size_t i = 0; i < indexWrap.size(); ++i)
        indexWrap[i] = i;        
      std::sort(indexWrap.begin(), indexWrap.end(), HashMap3DSearchIndexComparator(&distances));
      
      for(size_t i = 0; i < indexWrap.size(); ++i){
        sorted_indices[indexWrap[i]] = indices[i];
        sorted_distances[indexWrap[i]] = distances[i];
      }

      std::swap(sorted_indices, indices);
      std::swap(sorted_distances, distances);      
    }

    void nearestNeighborSearch(PointType& point, float& sqrDistance, std::vector<int>& indices, std::vector<float>& distances){
      indices.clear();
      distances.clear();
      Eigen::Vector3f p = point.getVector3fMap();
      Eigen::Vector3i index = worldToGrid(p);
      int searchRange = lrint(sqrt(sqrDistance) / _gridSize);
      float minDist = sqrDistance + 1e-5;
      int minId = -1;
      for(int x = -searchRange; x <= searchRange; ++x)
        for(int y = -searchRange; y <= searchRange; ++y)
          for(int z = -searchRange; z <= searchRange; ++z){
            Eigen::Vector3i cellIndex(x,y,z);
            cellIndex += index;
            std::vector<int>& cellIndices = cell(cellIndex);
            Eigen::Vector3f cp;
            ///check if min distance to bin is larger that already found best candidate distance
            Eigen::Vector3f cellCenter = gridToWorld(cellIndex);
            float centerdist =  std::max(0.0f, (cellCenter - p).norm() - (1.5f * _gridSize));
            centerdist  *= centerdist;
            
            if(centerdist < minDist)
            for(size_t i = 0; i < cellIndices.size(); ++i){
              cp = pointcloud->points[cellIndices[i]].getVector3fMap();
              float sqrDist = (p(0) - cp(0)) * (p(0) - cp(0)) + (p(1) - cp(1)) * (p(1) - cp(1)) + (p(2) - cp(2)) * (p(2) - cp(2));
              //(p - pointcloud->points[cellIndices[i]].getVector3fMap()).squaredNorm();
              if(sqrDist < minDist){
                minDist = sqrDist;
                minId = cellIndices[i];
              }
            }
          }
      if(minId > 0){
        indices.push_back(minId);
        distances.push_back(minDist);
      }
    }
    
    
    typename pcl::PointCloud< PointType >::Ptr    pointcloud;
    std::tr1::unordered_map<int, std::tr1::unordered_map<int, std::tr1::unordered_map<int, bool > > > valid;
  };

} //end namespace
#endif