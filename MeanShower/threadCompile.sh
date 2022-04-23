#!/bin/bash

g++ -std=c++11 -pthread -c ../CalModel/CalGeo.cxx
g++ -std=c++11 -pthread -c ../MathMethods/MatterCalc.cxx
g++ -std=c++11 -pthread meanShower.cxx -c `root-config --cflags --glibs` 
g++ -std=c++11 -pthread CalGeo.o MatterCalc.o meanShower.o -o meanShower.out `root-config --cflags --glibs` 
