# CMake file for creating a library containing MRF algorthims

SET(MAXFLOW_SRC
   graph.cpp
   maxflow.cpp
)

ADD_LIBRARY(maxflow-3.0 ${MAXFLOW_SRC})
