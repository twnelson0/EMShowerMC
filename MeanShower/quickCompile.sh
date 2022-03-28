#!/bin/bash

#I know I should be using CMake but I'm too lazy
g++ meanShower.cxx -o MeanTest.out `root-config --cflags --glibs`
