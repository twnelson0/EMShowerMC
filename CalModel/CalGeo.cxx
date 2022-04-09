#include "CalGeo.h"
#include <iostream>
#include <cmath>


double *layerTerminus(double startPoint){
	double[2] layarArr;
	double startFloor, layerStart, layerEnd;
	int layerMod;

	//Determine end points of layer
	startFloor = std::floor(startPoint); //Get floor of start point
	layerMod = (int) startFloor % 6; //Get modulus 6 of start point (roughly)

	layerArr[0] = startFloor - (double) layerMod; //Get start point of layer
	layerArr[1] = startFloor + (double) (6 - layerMod); //Get end point of layer

	return layerArr;
}

double layerTrackLen(double edgeVal, bool start = true){
	double trackLen;
	double *layerEnds = layerTerminus(edgeVal); //Get the start and end points of the layer

	//Determine material of current layer (and hence track length)
	if (start){ //If input value is the start point
		if (edgeVal >= *(layerEnds) && edgeVal <= *(layerEnds + 1))trackLen = 4; //Particle in lead
		else trackLen = *(layerEnds + 1) - edgeVal; //How much scintilator the particle will propogate through
	}
	
	//If input value is the end
	if (!start){ 
		if (edgeVal >= *(layerEnds) && edgeVal <= *(layerEnds + 1))trackLen = 0; //Particle in lead
		else trackLen = *(layerEnds + 0) + edgeVal; //How much scintilator the particle will propogate through
	}
	
	
	return trackLen;	
}

double trackLen(double startPoint, double endPoint){
	double *startLayer = layerTerminus(startPoint); //Get first Layer information
	double *endLayer = layerTerminus(endPoint); //Get Ending Layer information
	double trackLen;

	//Get track length through first layer (likely not the full layer)
	trackLen = layerTrackLen(startPoint);

	//Get track length through middle layers
	trackLen += (*(endLayer + 0) - *(startLayer + 1))/6 * 4; //Scintilate through all the middle layers

	//Get track length through final layer (likely not the full layer)
	trackLen += layerTrackLen(endPoint, false);

	return trackLen;

}