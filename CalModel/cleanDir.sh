#!/bin/bash

if [[ -e CalTest.out ]]; then
	echo Removing Compiled files
	rm $(pwd)/*.out
	rm $(pwd)/*.o
fi
