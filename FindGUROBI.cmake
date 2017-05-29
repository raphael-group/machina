# copied from https://searchcode.com/codesearch/raw/30585458/
FIND_PATH(GUROBI_INCLUDE_DIR
          NAMES "gurobi_c++.h" "gurobi_c.h"
          PATHS /Library/gurobi651/mac64/include/ /usr/local/gurobi651/linux64/include/
          DOC "Gurobi include directory")

FIND_LIBRARY(GUROBI_CPP_LIB
             NAMES gurobi_c++ 
             PATHS /Library/gurobi651/mac64/lib/ /usr/local/gurobi651/linux64/lib/
             DOC "Gurobi C++ Libraries")

FIND_LIBRARY(GUROBI_LIB
             NAMES gurobi65
             PATHS /Library/gurobi651/mac64/lib/ /usr/local/gurobi651/linux64/lib/
             DOC "Gurobi C Libraries")

set(GUROBI_LIBRARIES ${GUROBI_CPP_LIB} ${GUROBI_LIB})

set(GUROBI_FOUND TRUE)

