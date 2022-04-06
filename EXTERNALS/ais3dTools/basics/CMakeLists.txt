#INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR} ${QT_QTCORE_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${QT_QTCORE_INCLUDE_DIR}/..)

ADD_LIBRARY(basics
  macros.h
  misc.h
  os_specific.h os_specific.c
  runtime_error.h runtime_error.cpp
  string_tools.h string_tools.cpp
  filesys_tools.h filesys_tools.cpp
  timeutil.h timeutil.cpp
  ais3dTools_stuff_api.h
  getopt.h getopt.c
  parameters_manager.h parameters_manager.cpp
  pngIO.h pngIO.cpp
)
IF (UNIX)
  TARGET_LINK_LIBRARIES(basics rt)
ENDIF(UNIX)

SET_TARGET_PROPERTIES(basics PROPERTIES OUTPUT_NAME ${LIB_PREFIX}basics)

TARGET_LINK_LIBRARIES(basics png)

INSTALL(TARGETS basics
  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
)

#FILE(GLOB headers "${CMAKE_CURRENT_SOURCE_DIR}/*.h" "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")

#INSTALL(FILES ${headers} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/ssa/alu_glwidget)
