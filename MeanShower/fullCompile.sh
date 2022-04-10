#!/bin/bash

g++ -c ../CalModel/CalGeo.cxx 
g++ meanShower.cxx -c `root-config --cflags --glibs`
g++ CalGeo.o meanShower.o -o meanShower.out `root-config --cflags --glibs`