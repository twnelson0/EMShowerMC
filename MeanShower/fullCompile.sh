#!/bin/bash

g++ -c ../CalModel/CalGeo.cxx
g++ -c ../MathMethods/MatterCalc.cxx
g++ meanShower.cxx -c `root-config --cflags --glibs`
g++ CalGeo.o MatterCalc.o meanShower.o -o meanShower.out `root-config --cflags --glibs`
