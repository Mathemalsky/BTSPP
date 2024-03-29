# BSTPP

This software contains the implementation of approximation algorithms for 3 variants of the bottleneck travelling salesman problem. Those are the bottleneck travelling salesman problem (BTSP), the bottleneck travelling salesman path problem (BTSPP) and the bottleneck travelling salesman variable path problem (BTSVPP).
In addition, methods to solve TSP, BTSP and BTSPP using a MILP solver are also included.
This program is part of my master thesis.

## Using the whole program

The whole program includes the visualisation for examples on euclidean graphs.

### Prerequisites
- [cmake](https://cmake.org/) version 3.20 or higher is required for compiling
- [dear imgui](https://github.com/ocornut/imgui) is used for graphical output. dear imgui is contained as git submodule, but it needs
- [glfw](https://www.glfw.org/) can be installed using the following steps:
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
- [Eigen3](https://eigen.tuxfamily.org/) Is used as an implementation for sparse matrices. (I'm using Eigen 3.4)
- [HiGHS](https://www.maths.ed.ac.uk/hall/HiGHS/#top)

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
`C`   | toggle collapse mode of settings window
`R`   | recompute solution(s)
`T`   | switch display of approximation to next mode
`1`   | toggle drawing of BTSP approximation
`2`   | toggle drawing of BTSPP approximation
`3`   | toggle drawing of BTSVPP approximation
`4`   | toggle drawing of BTSP exact solution
`5`   | toggle drawing of BTSPP exact solution
`6`   | toggle drawing of TSP exact solution

## Using only the command line program

### Prerequisites
- [cmake](https://cmake.org/) version 3.20 or higher is required for compiling
- [doxygen](https://www.doxygen.nl/) This is optional for generating documentation.
- [Eigen3](https://eigen.tuxfamily.org/) is used as an implementation for sparse matrices. (I'm using Eigen 3.4)
- [HiGHS](https://www.maths.ed.ac.uk/hall/HiGHS/#top)

For installation advice see prerequisites for the whole program.

### Building the software
```
mkdir build && cd build
cmake -DVisualisation=Off ..
make
```

### Running
To run the application type:
`./<NameOfTheExecutable> <NumberOfNodesInTheGraph> <arg1> <arg2> ...`
Possible arguments are:

argument                              | effect
--------------------------------------|------------------
`-btsp`                               | approximates BTSP
`-btspp`                              | approximates BTSPP
`-btsvpp`                             | approximates BTSVPP
`-btsp-e`                             | solves exact BTSP
`-btspp-e`                            | solves exact BTSPP
`-tsp-e`                              | solves exact TSP
`-seed> <int1> ... `                  | set a seed for random generation of graph
`-no-crossing`                        | only if `btsp-e` is set: set extra constraint, that solutions cannot contain crossings
`-logfile:=<filename>`                | specifies a file to write stats to
`-repetitions:=<numberOfRepetitions>` | specifies the number of repetitions
`-no-info`                            | suppress output of calculation information to console
`-no-seed`                            | suppress output of seed to console