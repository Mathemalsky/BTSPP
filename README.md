# BBSTP

## Using the whole program

### Prerequisites
- [cmake](https://cmake.org/) version 3.20 or higher is required for compiling
- [dear imgui](https://github.com/ocornut/imgui) is used for graphical output. dear imgui is contained as git submodule, but it needs
- [glfw](https://www.glfw.org/) Can be installed using the following steps:
```
sudo apt-get install xorg-dev
sudo apt-get install libgl1-mesa-dev
git clone https://github.com/glfw/glfw.git
cd glfw && mkdir build && cd build
cmake -G "Unix Makefiles" ..
make
sudo make install
```
- [glew](https://github.com/nigels-com/glew) (can be installed by `sudo apt-get install libglew-dev`)
- [doxygen](https://www.doxygen.nl/) This is optional for generating documentation.
- [Eigen3](https://eigen.tuxfamily.org/) Used as implementation for sparse matrices. (I'm using Eigen 3.4)
- [HiGHS](https://www.maths.ed.ac.uk/hall/HiGHS/#top)
Here you probably need to slightly differ from the standard installation described on their website.
```
git clone https://github.com/ERGO-Code/HiGHS.git
cd HiGhS
mkdir build && cd build
cmake -DFAST_BUILD=ON -DCMAKE_INSTALL_PREFIX=/usr/local/ -DCMAKE_TARGETS=ON -DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -O3 -march=native -fPIC" -DCMAKE_C_FLAGS="${CMAKE_C_FLAGS} -O3 -march=native -fPIC" ..
cmake --build .
ctest
sudo cmake --install .
```
`-fPIC` is necessary, `-O3` and `-march=native` are recommended

Of course `/usr/local/` can be replaced by any other path contained in `CMAKE_SYSTEM_PREFIX_PATH` variable.

### Building the software
```
mkdir build && cd build
cmake ..
make
```

### Running
To launch the program run:
`./<NameOfTheExecutable> <NumberOfNodesInTheGraph>`

Then you can switch modes via the mouse in the settings window or using these shortcuts:

key   | function
------|-------------
`ESC` | close application
`F3`  | toggle visibility of settings window
`R`   | recompute solution(s)
`T`   | switch display of approximation to next mode
`1`   | toggle drawing of BTSP approximation
`2`   | toggle drawing of BTSPP approximation
`3`   | toggle drawing of BTSP exact solution
`4`   | toggle drawing of BTSPP exact solution
`5`   | toggle drawing of TSP exact solution

## Using only the command line program

### Prerequisites
- [cmake](https://cmake.org/) version 3.20 or higher is required for compiling
- [doxygen](https://www.doxygen.nl/) This is optional for generating documentation.
- [Eigen3](https://eigen.tuxfamily.org/) Used as implementation for sparse matrices. (I'm using Eigen 3.4)
- [HiGHS](https://www.maths.ed.ac.uk/hall/HiGHS/#top)

For installation advices see prerequisites for the whole program.

### Building the software
```
mkdir build && cd build
cmake -DVisualization=Off ..
make
```

### Running
To run the application type:
`./<NameOfTheExecutable> <NumberOfNodesInTheGraph> <arg1> <arg2> ...`
Possible arguments are:

argument  | effect
----------|------------------
`-btsp`   | approximates BTSP
`-btspp`  | approximates BTSPP
`-btsp-e` | solves exact BTSP
`-btspp-e`| solves exact BTSPP
`-tsp-e`  | solves exact TSP
