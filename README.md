#mFINCH - migration Framework for INferring Cancer Histories

![overview.png](doc/overview.png)

mFINCH is a omputational framework for inferring migration patterns between a primary tumor and metastases using DNA sequencing data.

##Compilation instructions

###Dependencies

mFINCH is written in C++11 and as such requires a modern C++ compiler
(GCC >= 4, or CLANG). In addition, mFINCH has the following
dependencies.

* [CMake](http://www.cmake.org/) (>= 2.8)
* [Boost](http://www.boost.org) (>= 1.38)
* [LEMON](http://lemon.cs.elte.hu/trac/lemon) graph library (>= 1.3)
* [Gurobi](http://www.gurobi.com) (>= 6.0)

[Graphviz](http://www.graphviz.org) is required to visualize the resulting DOT files, but is not required for compilation.

In case [doxygen](http://www.stack.nl/~dimitri/doxygen/) is available, extended source code documentation will be generated.

###Compilation

To compile SPRUCE, execute the following commands from the root of the repository:

    mkdir build
    cd build
    cmake ..
    make
    
In case CMake fails to detect LEMON, run the following command with adjusted paths:

	cmake -DLIBLEMON_ROOT=~/lemon ..
	
The compilation results in the following files in the `build` directory:

