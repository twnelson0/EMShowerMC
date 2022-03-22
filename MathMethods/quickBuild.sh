#!/bin/bash

g++ -c MatterCalc.cxx
g++ -c testBed.cxx
g++ MatterCalc.o testBed.o -o mathTest.out


