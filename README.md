# EMShowerMC
Monte Carlo Electromagnetic Shower Code for the final project of Physics 736.

## Structure
This repo is divided into a few files, a directory with dependencies specific to this program, a python prototyping area, a mean MC model and a more rigorous MC model.


## Running
To compile with ROOT libaries in g++ run
```
g++ script.C -o outName.out `root-config --cflags --glibs
```


## CMake Compiling
This project may be convoluted enough to require the use of CMake
