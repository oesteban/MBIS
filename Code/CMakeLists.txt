PROJECT(MBIS)

SET (MBIS_VERSION_MAJOR 1)
SET (MBIS_VERSION_MINOR 0)
SET (MBIS_RELEASE "rc1")

configure_file (
  "${PROJECT_SOURCE_DIR}/MBISConfig.h.in"
  "${PROJECT_BINARY_DIR}/MBISConfig.h"
)

cmake_minimum_required(VERSION 2.6)

FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR
          "Please Set ITK_DIR")
ENDIF(ITK_FOUND)

IF (WIN32)
	SET( Boost_NO_SYSTEM_PATHS ON )
	SET ( BOOST_ROOT "C:/Boost" )
	SET ( BOOST_INCLUDEDIR "C:/Boost/include" )
	# SET ( Boost_ADDITIONAL_VERSIONS "1.46" "1.46.1")
	#MESSAGE( "BOOST_ROOT=${BOOST_ROOT}")
ENDIF(WIN32)


#SET( Boost_USE_STATIC_LIBS ON )
FIND_PACKAGE ( Boost COMPONENTS system filesystem program_options REQUIRED)

SET ( MAXFLOW_USE_RESEARCH ON )

IF(MAXFLOW_USE_RESEARCH)
	SET ( MAXFLOW_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/maxflow-3.0 )
	SET ( MAXFLOW_LIBRARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/ThirdParty/maxflow-3.0 )
	ADD_SUBDIRECTORY(ThirdParty/maxflow-3.0)
	MESSAGE( WARNING "Using maxflow-3.0. Be advised that this library is only for research purposes" )
	SET ( MAXFLOW_LIB "maxflow-3.0" )
ELSE(MAXFLOW_USE_RESEARCH) 
	SET ( MAXFLOW_ALG "adjacency_list" )
	SET ( MAXFLOW_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/maxflow-2.21/${MAXFLOW_ALG} )
	SET ( MAXFLOW_LIBRARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/ThirdParty/maxflow-2.21 )
	ADD_SUBDIRECTORY(ThirdParty/maxflow-2.21/${MAXFLOW_ALG})
	SET ( MAXFLOW_LIB "maxflow-2.21" )
ENDIF(MAXFLOW_USE_RESEARCH)

INCLUDE_DIRECTORIES(
	${CMAKE_CURRENT_SOURCE_DIR}
	${MAXFLOW_SOURCE_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}/Numerics/Statistics/include
	${CMAKE_CURRENT_SOURCE_DIR}/Numerics/KdTree/include
	${CMAKE_CURRENT_SOURCE_DIR}/Optimizers/include
	${CMAKE_CURRENT_SOURCE_DIR}/Filters/MRF/include
	${CMAKE_CURRENT_SOURCE_DIR}/Filters/Segmentation/include
	${CMAKE_CURRENT_SOURCE_DIR}/Filters/ImageIntensity/include
	${CMAKE_CURRENT_SOURCE_DIR}/Filters/Statistics/include
	${CMAKE_CURRENT_SOURCE_DIR}/Common/include
	${PROJECT_BINARY_DIR}
)

# ITK4
LINK_DIRECTORIES( ${ITK_LIBRARY_PATH} ${MRF_LIBRARY_DIR} )

ADD_EXECUTABLE(brain_seg brain_seg.cpp)
TARGET_LINK_LIBRARIES( brain_seg
        ${ITK_LIBRARIES}
        ${Boost_LIBRARIES}
        ${MAXFLOW_LIB}
)
 
SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_CRT_SECURE_NO_DEPRECATE -D_SCL_SECURE_NO_WARNINGS")
