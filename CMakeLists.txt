cmake_minimum_required( VERSION 2.8 )

project( MACHINA )

set( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR} ${CMAKE_MODULE_PATH} )

set( generatemigrationtrees_src
  src/generatemigrationtreesmain.cpp
  src/utils.cpp
  src/nonbinaryclonetree.cpp
  src/clonetree.cpp
  src/binarytree.cpp
  src/basetree.cpp
  src/migrationtree.cpp
  src/gabowmyers.cpp
)

set( generatemigrationtrees_hdr
  src/utils.h
  src/nonbinaryclonetree.h
  src/clonetree.h
  src/binarytree.h
  src/basetree.h
  src/migrationtree.h
  src/gabowmyers.h
)

set( pmh_hdr
  src/utils.h
  src/nonbinaryclonetree.h
  src/clonetree.h
  src/binarytree.h
  src/basetree.h
  src/ilpsolver.h
  src/migrationgraph.h
)

set( pmh_src
  src/pmh.cpp
  src/utils.cpp
  src/clonetree.cpp
  src/nonbinaryclonetree.cpp
  src/binarytree.cpp
  src/basetree.cpp
  src/ilpsolver.cpp
  src/migrationgraph.cpp
  src/utils.h
  src/nonbinaryclonetree.h
  src/clonetree.h
  src/binarytree.h
  src/basetree.h
  src/ilpsolver.h
  src/migrationgraph.h
)

set( pmh_pr_hdr
  src/utils.h
  src/nonbinaryclonetree.h
  src/clonetree.h
  src/binarytree.h
  src/basetree.h
  src/ilpsolver.h
  src/ilpbinarizationsolver.h
  src/migrationgraph.h
)

set( pmh_pr_src
  src/pmh_pr.cpp
  src/utils.cpp
  src/clonetree.cpp
  src/nonbinaryclonetree.cpp
  src/binarytree.cpp
  src/basetree.cpp
  src/ilpsolver.cpp
  src/ilpbinarizationsolver.cpp
  src/migrationgraph.cpp
  src/utils.h
  src/nonbinaryclonetree.h
  src/clonetree.h
  src/binarytree.h
  src/basetree.h
  src/ilpsolver.h
  src/ilpbinarizationsolver.h
  src/migrationgraph.h
)

set( pmh_cti_hdr
  src/utils.h
  src/nonbinaryclonetree.h
  src/clonetree.h
  src/binarytree.h
  src/basetree.h
  src/ilpsolver.h
  src/migrationgraph.h
  src/frequencymatrix.h
  src/gabowmyers.h
)

set( pmh_cti_src
  src/pmh_cti.cpp
  src/utils.cpp
  src/clonetree.cpp
  src/nonbinaryclonetree.cpp
  src/binarytree.cpp
  src/basetree.cpp
  src/ilpsolver.cpp
  src/ilpsolverext.cpp
  src/migrationgraph.cpp
  src/frequencymatrix.cpp
  src/gabowmyers.cpp
  src/utils.h
  src/nonbinaryclonetree.h
  src/clonetree.h
  src/binarytree.h
  src/basetree.h
  src/ilpsolverext.h
  src/ilpsolver.h
  src/migrationgraph.h
  src/gabowmyers.h
)

set( visualizeclonetree_hdr
  src/utils.h
  src/nonbinaryclonetree.h
  src/clonetree.h
  src/binarytree.h
  src/basetree.h
)

set( visualizeclonetree_src
  src/visualizeclonetreemain.cpp
  src/utils.cpp
  src/clonetree.cpp
  src/nonbinaryclonetree.cpp
  src/binarytree.cpp
  src/basetree.cpp
)

set( visualizemigrationgraph_hdr
  src/migrationgraph.h
  src/utils.h
  src/nonbinaryclonetree.h
  src/clonetree.h
  src/binarytree.h
  src/basetree.h
  src/migrationtree.h
  src/gabowmyers.h
)

set( visualizemigrationgraph_src
  src/migrationgraph.cpp
  src/utils.cpp
  src/nonbinaryclonetree.cpp
  src/clonetree.cpp
  src/binarytree.cpp
  src/basetree.cpp
  src/migrationtree.cpp
  src/gabowmyers.cpp
  src/visualizemigrationgraphmain.cpp
)

set( pmh_sankoff_hdr
  src/basetree.h
  src/binarytree.h
  src/utils.h
  src/mpbase.h
  src/sankoff.h
  src/charactertree.h
  src/nonbinaryclonetree.h
  src/clonetree.h
  src/sankofflabeling.h
  src/migrationgraph.h
)

set( pmh_sankoff_src
  src/basetree.cpp
  src/binarytree.cpp
  src/utils.cpp
  src/mpbase.cpp
  src/sankoff.cpp
  src/charactertree.cpp
  src/pmh_sankoff.cpp
  src/nonbinaryclonetree.cpp
  src/sankofflabeling.cpp
  src/clonetree.cpp
  src/migrationgraph.cpp
)

# first we can indicate the documentation build as an option and set it to ON by default
option(BUILD_DOC "Build documentation" ON)

#if (CMAKE_BUILD_TYPE MATCHES "^[Rr]elease")
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
#endif()

set( LIBLEMON_ROOT "$ENV{HOME}/lemon" CACHE PATH "Additional search directory for lemon library" )

set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11" )

if( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-local-typedefs" )
endif()

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

add_executable( pmh_pr ${pmh_pr_src} ${pmh_pr_hdr} )
target_include_directories( pmh_pr PUBLIC "${LIBLEMON_ROOT}/include" "src" ${Boost_INCLUDE_DIRS} ${GUROBI_INCLUDE_DIR} )
target_link_libraries( pmh_pr ${CommonLibs} ${GUROBI_LIBRARIES} )

add_executable( pmh_cti ${pmh_cti_src} ${pmh_cti_hdr} )
target_include_directories( pmh_cti PUBLIC "${LIBLEMON_ROOT}/include" "src" ${Boost_INCLUDE_DIRS} ${GUROBI_INCLUDE_DIR} )
target_link_libraries( pmh_cti ${CommonLibs} ${GUROBI_LIBRARIES} )

add_executable( generatemigrationtrees ${generatemigrationtrees_src} ${generatemigrationtrees_hdr} )
target_link_libraries( generatemigrationtrees ${CommonLibs} )

add_executable( visualizemigrationgraph ${visualizemigrationgraph_src} ${visualizemigrationgraph_hdr} )
target_link_libraries( visualizemigrationgraph ${CommonLibs} )

add_executable( visualizeclonetree ${visualizeclonetree_src} ${visualizeclonetree_hdr} )
target_link_libraries( visualizeclonetree ${CommonLibs} )
