# CMake file for creating a library containing MRF algorthims

SET(MAXFLOW_SRC graph.cpp maxflow.cpp)

#include_directories( "${MAXFLOW_ALG_PATH}/" )

ADD_LIBRARY(maxflow-2.21 ${MAXFLOW_SRC})
