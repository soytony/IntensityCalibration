CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

project(remissionCalibration)

# The library prefix
SET(LIB_PREFIX ${PROJECT_NAME}_)

SET(${PROJECT_NAME}_C_FLAGS)
SET(${PROJECT_NAME}_CXX_FLAGS)

# -Wno-error="deprecated-declarations deprecated-copy" 
# disable warnings
set(CMAKE_CXX_FLAGS 
  "${CMAK_CXX_FLAGS} -Wfatal-errors -Wno-deprecated-declarations -Wno-deprecated-copy")

# default built type
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

SET(${PROJECT_NAME}_LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib CACHE PATH "Output path for libraries")
SET(${PROJECT_NAME}_EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin CACHE PATH "Output path for executables")

# Set the output directory for the build executables and libraries
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${${PROJECT_NAME}_LIBRARY_OUTPUT_PATH})
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${${PROJECT_NAME}_EXECUTABLE_OUTPUT_PATH})
IF(WIN32)
  SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${${PROJECT_NAME}_EXECUTABLE_OUTPUT_PATH})
ELSE()
  SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${${PROJECT_NAME}_LIBRARY_OUTPUT_PATH})
ENDIF()

#Build type
if (UNIX)
  SET(BUILD_SHARED_LIBS ON)
ELSEIF (WIN32)
  SET(BUILD_SHARED_LIBS OFF)
ENDIF()


# PCL
find_package(PCL 1.5 REQUIRED)
include_directories(SYSTEM ${PCL_INCLUDE_DIRS})
LINK_DIRECTORIES(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

ADD_SUBDIRECTORY(EXTERNALS)
add_subdirectory(src/BackwardCpp)

INCLUDE_DIRECTORIES(${G2O_INCLUDE_DIRECTORY})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/EXTERNALS)
LINK_DIRECTORIES(${${PROJECT_NAME}_LIBRARY_OUTPUT_PATH})

ADD_LIBRARY(remissionCalibrationLib
  src/g2oGraphElements.h src/g2oGraphElements.hpp src/g2oGraphElements.cpp
  src/pointXYZIVPL.h
  src/remissionCalibrationHelper.h src/remissionCalibrationHelper.hpp src/remissionCalibrationHelper.cpp
  src/remissionCalibrationResult.h src/remissionCalibrationResult.hpp src/remissionCalibrationResult.cpp
)
TARGET_LINK_LIBRARIES(remissionCalibrationLib
  ${PCL_LIBRARIES}
  g2o_core g2o_solver_cholmod g2o_solver_csparse g2o_solver_dense
  #-lg2o_csparse_extension -lg2o_ext_freeglut_minimal -lg2o_opengl_helper -lg2o_solver_csparse -lg2o_solver_dense -lg2o_solver_pcg -lg2o_solver_slam2d_linear -lg2o_solver_structure_only -lg2o_stuff -lg2o_types_data -lg2o_types_icp -lg2o_types_sba -lg2o_types_sclam2d -lg2o_types_sim3 -lg2o_types_slam2d -lg2o_types_slam3d
  ais3dTools_basics
)
# ADD_DEPENDENCIES(remissionCalibrationLib g2o ais3dTools)

ADD_EXECUTABLE(mapCellsExtractor 
  src/mapCellsExtractor.cpp
  ${BACKWARD_ENABLE})
TARGET_LINK_LIBRARIES(mapCellsExtractor 
  remissionCalibrationLib
  dw)

ADD_EXECUTABLE(remissionCalibration 
  src/remissionCalibration.cpp
  ${BACKWARD_ENABLE})
TARGET_LINK_LIBRARIES(remissionCalibration
  remissionCalibrationLib
  dw)

ADD_EXECUTABLE(statisticsGenerator src/statisticsGenerator.cpp)
TARGET_LINK_LIBRARIES(statisticsGenerator remissionCalibrationLib)
