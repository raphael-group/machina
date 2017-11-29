cmake_minimum_required( VERSION 2.8 )

project( MACHINA )

set( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR} ${CMAKE_MODULE_PATH} )

set( generatemigrationtrees_src
  src/generatemigrationtreesmain.cpp
  src/utils.cpp
  src/clonetree.cpp
  src/basetree.cpp
  src/migrationtree.cpp
  src/gabowmyers.cpp
)

set( generatemigrationtrees_hdr
  src/utils.h
  src/clonetree.h
  src/basetree.h
  src/migrationtree.h
  src/gabowmyers.h
)

set( pmh_hdr
  src/utils.h
  src/clonetree.h
  src/old_ilps/ilpsolver.h
  src/ilppmhsolver.h
  src/migrationgraph.h
  src/migrationtree.h
  src/gabowmyers.h
)

set( simulation_hdr
  src/simulation/cell.h
  src/simulation/simulation.h
  src/simulation/beta_distribution.hpp
  src/clonetree.h
  src/basetree.h
  src/migrationtree.h
  src/migrationgraph.h
  src/gabowmyers.h
  src/utils.h
)

set( simulation_src
  src/utils.cpp
  src/clonetree.cpp
  src/basetree.cpp
  src/migrationtree.cpp
  src/migrationgraph.cpp
  src/gabowmyers.cpp
  src/simulation/cell.cpp
  src/simulation/simulationmain.cpp
  src/simulation/simulation.cpp
)

set( pmh_src
  src/pmh.cpp
  src/utils.cpp
  src/clonetree.cpp
  src/basetree.cpp
  src/old_ilps/ilpsolver.cpp
  src/ilppmhsolver.cpp
  src/migrationgraph.cpp
  src/migrationtree.cpp
  src/utils.h
  src/clonetree.h
  src/basetree.h
  src/old_ilps/ilpsolver.h
  src/migrationgraph.h
  src/gabowmyers.cpp
)

set( pmh_tr_hdr
  src/utils.h
  src/clonetree.h
  src/basetree.h
  src/old_ilps/ilpsolver.h
  src/old_ilps/ilpbinarizationsolver.h
  src/ilppmhsolver.h
  src/ilppmhtrsolver.h
  src/ilppmhtrsolvercallback.h
  src/migrationgraph.h
  src/migrationtree.h
  src/gabowmyers.h
  src/solutiongraph.h
)

set( pmh_tr_hdr
  src/pmh_tr.cpp
  src/utils.cpp
  src/clonetree.cpp
  src/basetree.cpp
  src/old_ilps/ilpsolver.cpp
  src/old_ilps/ilpsolver.h
  src/old_ilps/ilpbinarizationsolver.cpp
  src/old_ilps/ilpbinarizationsolver.h
  src/ilppmhtrsolvercallback.cpp
  src/ilppmhtrsolver.cpp
  src/ilppmhtrsolver.h
  src/ilppmhsolver.cpp
  src/ilppmhsolver.h
  src/solutiongraph.h
  src/solutiongraph.cpp
  src/migrationgraph.cpp
  src/migrationtree.cpp
  src/utils.h
  src/clonetree.h
  src/basetree.h
  src/migrationgraph.h
  src/gabowmyers.cpp
)

set( generatemutationtrees_hdr
  src/utils.h
  src/basematrix.h
  src/frequencymatrix.h
  src/clonetree.h
  src/enumeratecanonicalclonetrees.h
  src/enumeratemutationtrees.h
  src/spruce/matrix.h
  src/spruce/realmatrix.h
  src/spruce/realtensor.h
  src/spruce/rootedcladisticancestrygraph.h
  src/spruce/rootedcladisticnoisyancestrygraph.h
  src/spruce/rootedcladisticenumeration.h
  src/spruce/rootedcladisticnoisyenumeration.h
  src/spruce/solution.h
  src/spruce/solutionset.h
  src/spruce/statetree.h
  src/spruce/tensor.h
  src/rootedcladisticnoisysparseenumeration.h
)

set( generatemutationtrees_src
  src/generatemutationtreesmain.cpp
  src/utils.cpp
  src/basematrix.cpp
  src/frequencymatrix.cpp
  src/enumeratecanonicalclonetrees.cpp
  src/enumeratemutationtrees.cpp
  src/utils.h
  src/basematrix.h
  src/frequencymatrix.h
  src/basetree.h
  src/basetree.cpp
  src/clonetree.h
  src/clonetree.cpp
  src/rootedcladisticnoisysparseenumeration.h
  src/rootedcladisticnoisysparseenumeration.cpp
  src/spruce/matrix.h
  src/spruce/perfectphylograph.h
  src/spruce/perfectphylomatrix.h
  src/spruce/perfectphylotree.h
  src/spruce/realmatrix.h
  src/spruce/realtensor.h
  src/spruce/rootedcladisticancestrygraph.h
  src/spruce/rootedcladisticnoisyancestrygraph.h
  src/spruce/rootedcladisticenumeration.h
  src/spruce/rootedcladisticnoisyenumeration.h
  src/spruce/statetree.h
  src/spruce/tensor.h
  src/spruce/matrix.cpp
  src/spruce/perfectphylograph.cpp
  src/spruce/perfectphylomatrix.cpp
  src/spruce/perfectphylotree.cpp
  src/spruce/realmatrix.cpp
  src/spruce/realtensor.cpp
  src/spruce/rootedcladisticancestrygraph.cpp
  src/spruce/rootedcladisticnoisyancestrygraph.cpp
  src/spruce/rootedcladisticenumeration.cpp
  src/spruce/rootedcladisticnoisyenumeration.cpp
  src/spruce/solution.cpp
  src/spruce/solutionset.cpp
  src/spruce/statetree.cpp
  src/spruce/tensor.cpp
)

set( ms_hdr
  src/utils.h
  src/basematrix.h
  src/frequencymatrix.h
  src/clonetree.h
  src/spruce/matrix.h
  src/spruce/realmatrix.h
  src/spruce/realtensor.h
  src/spruce/rootedcladisticancestrygraph.h
  src/spruce/rootedcladisticnoisyancestrygraph.h
  src/spruce/rootedcladisticenumeration.h
  src/spruce/rootedcladisticnoisyenumeration.h
  src/spruce/solution.h
  src/spruce/solutionset.h
  src/spruce/statetree.h
  src/spruce/tensor.h
  src/rootedcladisticnoisysparseenumeration.h
)

set( ms_src
  src/ms.cpp
  src/utils.cpp
  src/basematrix.cpp
  src/frequencymatrix.cpp
  src/utils.h
  src/basematrix.h
  src/frequencymatrix.h
  src/migrationtree.h
  src/migrationtree.cpp
  src/basetree.h
  src/basetree.cpp
  src/gabowmyers.h
  src/gabowmyers.cpp
  src/msenumeration.h
  src/msenumeration.cpp
  src/clonetree.h
  src/clonetree.cpp
  src/rootedcladisticnoisysparseenumeration.h
  src/rootedcladisticnoisysparseenumeration.cpp
  src/migrationgraph.cpp
  src/migrationgraph.h
  src/spruce/matrix.h
  src/spruce/perfectphylograph.h
  src/spruce/perfectphylomatrix.h
  src/spruce/perfectphylotree.h
  src/spruce/realmatrix.h
  src/spruce/realtensor.h
  src/spruce/rootedcladisticancestrygraph.h
  src/spruce/rootedcladisticnoisyancestrygraph.h
  src/spruce/rootedcladisticenumeration.h
  src/spruce/rootedcladisticnoisyenumeration.h
  src/spruce/statetree.h
  src/spruce/tensor.h
  src/spruce/matrix.cpp
  src/spruce/perfectphylograph.cpp
  src/spruce/perfectphylomatrix.cpp
  src/spruce/perfectphylotree.cpp
  src/spruce/realmatrix.cpp
  src/spruce/realtensor.cpp
  src/spruce/rootedcladisticancestrygraph.cpp
  src/spruce/rootedcladisticnoisyancestrygraph.cpp
  src/spruce/rootedcladisticenumeration.cpp
  src/spruce/rootedcladisticnoisyenumeration.cpp
  src/spruce/solution.cpp
  src/spruce/solutionset.cpp
  src/spruce/statetree.cpp
  src/spruce/tensor.cpp
)

set( reads_to_frequencies_hdr
  src/utils.h
  src/basematrix.h
  src/frequencymatrix.h
  src/cluster/readmatrix.h
)

set( reads_to_frequencies_src
  src/utils.cpp
  src/basematrix.cpp
  src/frequencymatrix.cpp
  src/cluster/readsmain.cpp
  src/cluster/readmatrix.cpp
)

set( cluster_hdr
  src/utils.h
  src/basematrix.h
  src/frequencymatrix.h
  src/cluster/cluster.h
  src/cluster/readmatrix.h
)

set( cluster_src
  src/utils.cpp
  src/basematrix.cpp
  src/frequencymatrix.cpp
  src/cluster/clustermain.cpp
  src/cluster/cluster.cpp
  src/cluster/readmatrix.cpp
)

set( pmh_ti_hdr
  src/utils.h
  src/clonetree.h
  src/basetree.h
  src/old_ilps/ilpsolver.h
  src/migrationgraph.h
  src/basematrix.h
  src/frequencymatrix.h
  src/gabowmyers.h
  src/enumeratemutationtrees.h
  src/ilppmhtisolver.h
  src/ilppmhtrsolver.h
  src/ilppmhtrsolvercallback.h
  src/ilppmhsolver.h
)

set( pmh_ti_src
  src/pmh_ti.cpp
  src/ilppmhtisolver.cpp
  src/ilppmhtrsolvercallback.cpp
  src/ilppmhtrsolver.cpp
  src/ilppmhsolver.cpp
  src/utils.cpp
  src/clonetree.cpp
  src/basetree.cpp
  src/old_ilps/ilpsolver.cpp
  src/old_ilps/ilpsolverext.cpp
  src/migrationgraph.cpp
  src/basematrix.cpp
  src/frequencymatrix.cpp
  src/gabowmyers.cpp
  src/utils.h
  src/clonetree.h
  src/basetree.h
  src/old_ilps/ilpsolverext.h
  src/old_ilps/ilpsolver.h
  src/migrationgraph.h
  src/migrationtree.h
  src/migrationtree.cpp
  src/gabowmyers.h
  src/rootedcladisticnoisysparseenumeration.h
  src/rootedcladisticnoisysparseenumeration.cpp
  src/enumeratemutationtrees.cpp
  src/spruce/matrix.h
  src/spruce/perfectphylograph.h
  src/spruce/perfectphylomatrix.h
  src/spruce/perfectphylotree.h
  src/spruce/realmatrix.h
  src/spruce/realtensor.h
  src/spruce/rootedcladisticancestrygraph.h
  src/spruce/rootedcladisticnoisyancestrygraph.h
  src/spruce/rootedcladisticenumeration.h
  src/spruce/rootedcladisticnoisyenumeration.h
  src/spruce/statetree.h
  src/spruce/tensor.h
  src/spruce/matrix.cpp
  src/spruce/perfectphylograph.cpp
  src/spruce/perfectphylomatrix.cpp
  src/spruce/perfectphylotree.cpp
  src/spruce/realmatrix.cpp
  src/spruce/realtensor.cpp
  src/spruce/rootedcladisticancestrygraph.cpp
  src/spruce/rootedcladisticnoisyancestrygraph.cpp
  src/spruce/rootedcladisticenumeration.cpp
  src/spruce/rootedcladisticnoisyenumeration.cpp
  src/spruce/solution.cpp
  src/spruce/solutionset.cpp
  src/spruce/statetree.cpp
  src/spruce/tensor.cpp
)

set( stats_hdr
  src/utils.h
  src/clonetree.h
  src/basetree.h
  src/migrationgraph.h
)

set( stats_src
  src/statsmain.cpp
  src/utils.cpp
  src/clonetree.cpp
  src/basetree.cpp
  src/migrationgraph.cpp
)

set( visualizeclonetree_hdr
  src/utils.h
  src/clonetree.h
  src/basetree.h
)

set( visualizeclonetree_src
  src/visualizeclonetreemain.cpp
  src/utils.cpp
  src/clonetree.cpp
  src/basetree.cpp
)

set( robinson_foulds_hdr
  src/utils.h
  src/clonetree.h
  src/basetree.h
)

set( robinson_foulds_src
  src/rfmain.cpp
  src/utils.cpp
  src/clonetree.cpp
  src/basetree.cpp
)

set( visualizemigrationgraph_hdr
  src/migrationgraph.h
  src/utils.h
  src/clonetree.h
  src/basetree.h
  src/migrationtree.h
  src/gabowmyers.h
)

set( visualizemigrationgraph_src
  src/migrationgraph.cpp
  src/utils.cpp
  src/clonetree.cpp
  src/basetree.cpp
  src/migrationtree.cpp
  src/gabowmyers.cpp
  src/visualizemigrationgraphmain.cpp
)

set( pmh_sankoff_hdr
  src/basetree.h
  src/utils.h
  src/mpbase.h
  src/sankoff.h
  src/charactertree.h
  src/clonetree.h
  src/sankofflabeling.h
  src/migrationgraph.h
)

set( pmh_sankoff_src
  src/basetree.cpp
  src/utils.cpp
  src/mpbase.cpp
  src/sankoff.cpp
  src/charactertree.cpp
  src/pmh_sankoff.cpp
  src/clonetree.cpp
  src/sankofflabeling.cpp
  src/migrationgraph.cpp
)

# first we can indicate the documentation build as an option and set it to ON by default
option(BUILD_DOC "Build documentation" ON)

if (CMAKE_BUILD_TYPE MATCHES "^[Rr]elease")
  # check if Doxygen is installed
  find_package(Doxygen)
  if (DOXYGEN_FOUND)
      # set input and output files
      set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
      set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

      # request to configure the file
      configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
      message("Doxygen build started")

      # note the option ALL which allows to build the docs together with the application
      add_custom_target( doc_doxygen ALL
          COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
          COMMENT "Generating API documentation with Doxygen"
          VERBATIM )
  else (DOXYGEN_FOUND) 
    message("Doxygen need to be installed to generate the doxygen documentation")
  endif (DOXYGEN_FOUND)
endif()

set( LIBLEMON_ROOT "$ENV{HOME}/lemon" CACHE PATH "Additional search directory for lemon library" )

set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11" )

if( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-local-typedefs" )
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: Debug Release Profile" FORCE)
endif(NOT CMAKE_BUILD_TYPE)

set( GitCommand1 "symbolic-ref HEAD 2> /dev/null | cut -b 12-" )
execute_process( COMMAND git symbolic-ref HEAD  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
COMMAND tail -c -7 COMMAND tr -d "\\n" OUTPUT_VARIABLE GitBranch )
set( GitCommand2 "log --pretty=format:\\\"%h\\\" -1" )
execute_process( COMMAND git log --pretty=format:%h WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
COMMAND head -n 1  COMMAND tr -d "\\n" OUTPUT_VARIABLE GitHashtag )

set( Boost_USE_STATIC_LIBS ON )
find_package( Boost REQUIRED COMPONENTS thread system filesystem )

include_directories( "${LIBLEMON_ROOT}/include" "src" ${Boost_INCLUDE_DIRS} )
link_directories( "${LIBLEMON_ROOT}/lib" )

set( CommonLibs
  m
  emon
  pthread
  ${Boost_LIBRARIES}
)

find_package( GUROBI REQUIRED )

set( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG" )

add_executable( pmh_sankoff ${pmh_sankoff_src} ${pmh_sankoff_hdr} )
target_link_libraries( pmh_sankoff ${CommonLibs} )

add_executable( pmh ${pmh_src} ${pmh_hdr} )
target_include_directories( pmh PUBLIC "${LIBLEMON_ROOT}/include" "src" ${Boost_INCLUDE_DIRS} ${GUROBI_INCLUDE_DIR} )
target_link_libraries( pmh ${CommonLibs} ${GUROBI_LIBRARIES} )

add_executable( pmh_tr ${pmh_tr_src} ${pmh_tr_hdr} )
target_include_directories( pmh_tr PUBLIC "${LIBLEMON_ROOT}/include" "src" ${Boost_INCLUDE_DIRS} ${GUROBI_INCLUDE_DIR} )
target_link_libraries( pmh_tr ${CommonLibs} ${GUROBI_LIBRARIES} )

add_executable( pmh_ti ${pmh_ti_src} ${pmh_ti_hdr} )
target_include_directories( pmh_ti PUBLIC "${LIBLEMON_ROOT}/include" "src" ${Boost_INCLUDE_DIRS} ${GUROBI_INCLUDE_DIR} )
target_link_libraries( pmh_ti ${CommonLibs} ${GUROBI_LIBRARIES} )

add_executable( ms ${ms_src} ${ms_hdr} )
target_include_directories( ms PUBLIC "${LIBLEMON_ROOT}/include" "src" ${Boost_INCLUDE_DIRS} ${GUROBI_INCLUDE_DIR} )
target_link_libraries( ms ${CommonLibs} ${GUROBI_LIBRARIES})

add_executable( generatemigrationtrees ${generatemigrationtrees_src} ${generatemigrationtrees_hdr} )
target_link_libraries( generatemigrationtrees ${CommonLibs} )

add_executable( generatemutationtrees ${generatemutationtrees_src} ${generatemutationtrees_hdr} )
target_include_directories( generatemutationtrees PUBLIC "${LIBLEMON_ROOT}/include" "src" ${Boost_INCLUDE_DIRS} ${GUROBI_INCLUDE_DIR} )
target_link_libraries( generatemutationtrees ${CommonLibs} ${GUROBI_LIBRARIES})

add_executable( visualizemigrationgraph ${visualizemigrationgraph_src} ${visualizemigrationgraph_hdr} )
target_link_libraries( visualizemigrationgraph ${CommonLibs} )

add_executable( visualizeclonetree ${visualizeclonetree_src} ${visualizeclonetree_hdr} )
target_link_libraries( visualizeclonetree ${CommonLibs} )

add_executable( rf EXCLUDE_FROM_ALL ${robinson_foulds_src} ${robinson_foulds_hdr} )
target_link_libraries( rf ${CommonLibs} )

add_executable( simulate EXCLUDE_FROM_ALL ${simulation_src} ${simulation_hdr} )
target_link_libraries( simulate ${CommonLibs} )

add_executable( cluster ${cluster_src} ${cluster_hdr} )
target_link_libraries( cluster ${CommonLibs} )

add_executable( frequencies EXCLUDE_FROM_ALL ${reads_to_frequencies_src} ${reads_to_frequencies_hdr} )
target_link_libraries( frequencies ${CommonLibs} )

add_executable( stats ${stats_src} ${stats_hdr} )
target_link_libraries( stats ${CommonLibs} )
