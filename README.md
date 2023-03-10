# pathBSTP

## running the software

### prerequisites
- [cmake](https://cmake.org/) is required for compiling
- [dear imgui](https://github.com/ocornut/imgui) is used for graphical output. dear imgui is contained as git submodule but it needs
- [glfw](https://www.glfw.org/) which needs
- [glew](https://github.com/nigels-com/glew) (can be installed by `sudo apt-get install libglew-dev`)
- [doxygen](https://www.doxygen.nl/) This is optinal for generating documentation.
- [Eigen3](https://eigen.tuxfamily.org/) Used as implementation for sparse matrices. (I'm using Eigen 3.4)
- [HiGHS](https://www.maths.ed.ac.uk/hall/HiGHS/#top)
Here you probably need to slightly differ from the standard instalation described on their website.
```
git clone https://github.com/ERGO-Code/HiGHS.git
cd HiGhS
mkdir build && cd build
cmake -DFAST_BUILD=ON -DCMAKE_INSTALL_PREFIX=/usr/local/ -DCMAKE_TARGETS=ON -DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -O3 -march=native -fPIC" -DCMAKE_C_FLAGS="${CMAKE_C_FLAGS} -O3 -march=native -fPIC"..
cmake --build .
ctest
sudo cmake --install .
```
`-fPIC` is necessary, `-O3` and `-march=native` are recommended

Ofcourse `/usr/local/` can be replaced by any other path contained in `CMAKE_SYSTEM_PREFIX_PATH` variable.

### building the software
```
mkdir build && cd build
cmake ..
make
```

### running
key   | function
------|-------------
`ESC` | close application
`F3`  | toggle visibility of settings window
`R`   | recompute solution(s)
`T`   | switch display of aproximation to next mode
`1`   | toggle drawing of BTSP approximation
`3`   | toggle drawing of BTSP exact solution
`4`   | toggle drawing of BTSPP exact solution
`5`   | toggle drawing of TSP exact solution
