#!/bin/bash

echo Removing .o and .out files

if [[ -e MeanTest.out ]];then
	rm MeanTest.out
fi

if [[ -e meanShower.out ]]; then
	rm CalGeo.o
	rm meanShower.o
	rm meanShower.out
fi
