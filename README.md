# mFINCH - migration Framework for INferring Cancer Histories

mFINCH is a omputational framework for inferring migration patterns between a primary tumor and metastases using DNA sequencing data.
![Overview of mFINCH](doc/overview.png)

## Compilation instructions

### Dependencies

mFINCH is written in C++11 and thus requires a modern C++ compiler (gcc >= 4, or clang). In addition, mFINCH has the following dependencies.

* [CMake](http://www.cmake.org/) (>= 2.8)
* [Boost](http://www.boost.org) (>= 1.38)
* [LEMON](http://lemon.cs.elte.hu/trac/lemon) graph library (>= 1.3)
* [Gurobi](http://www.gurobi.com) (>= 6.0)

[Graphviz](http://www.graphviz.org) is required to visualize the resulting DOT files, but is not required for compilation.

In case [doxygen](http://www.stack.nl/~dimitri/doxygen/) is available, extended source code documentation will be generated.

### Compilation

To compile SPRUCE, execute the following commands from the root of the repository:

    mkdir build
    cd build
    cmake ..
    make

In case CMake fails to detect LEMON, run the following command with adjusted paths:

    cmake -DLIBLEMON_ROOT=~/lemon ..

The compilation results in the following files in the `build` directory:

EXECUTABLE | DESCRIPTION
-----------|-------------
`pmh_sankoff`  | Enumerates all minimum-migration vertex labelings given a clone tree
`pmh` | Solves the Parsimonious Migration History (PMH) problem under various topological constraints (PS, S, M or R) given a clone tree
`pmh-pr`    | Solves the Parsimonious Migration History with Polytomy Resolution (PMH-PR) under various topological constraints given a clone tree
`pmh-cti`     | Sorts solution trees by the fraction of common edges (solution with rank 0 is the most representative tree)
`generatemigrationtrees` | Generates all migration trees given anatomical site labels

## Usage instructions

### Input formats


### Output formats


