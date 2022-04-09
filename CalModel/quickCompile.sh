#!/bin/bash

#I should technically be using CMake but I'm Lazy so shell scripting it is!
g++ -c CalGeo.cxx 
g++ -c CalGeoTest.cxx
g++ CalGeo.o CalGeoTest.o -o CalTest.out
