# CMake file for creating a library containing MRF algorthims

SET(MAXFLOW_ALG_PATH "./adjacency_list" )

SET(MAXFLOW_SRC
   ${MAXFLOW_ALG_PATH}/graph.cpp
   ${MAXFLOW_ALG_PATH}/maxflow.cpp
)

include_directories( "${MAXFLOW_ALG_PATH}/" )

ADD_LIBRARY(maxflow-2.21 ${MAXFLOW_SRC})
