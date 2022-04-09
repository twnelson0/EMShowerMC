#include "CalGeo.h"
#include <iostream>
#include <cmath>


double* layerTerminus(double startPoint){
	static double layerArr[2];
	double startFloor, layerStart, layerEnd;
	int layerMod;

	//Determine end points of layer
	startFloor = std::floor(startPoint); //Get floor of start point
	layerMod = (int) startFloor % 6; //Get modulus 6 of start point (roughly)

	layerArr[0] = startFloor - (double) layerMod; //Get start point of layer
	layerArr[1] = startFloor + (double) (6 - layerMod); //Get end point of layer

	return layerArr;
}

double layerTrackLen_scint(double edgeVal, bool start){
	double trackLen;
	double *layerEnds = layerTerminus(edgeVal); //Get the start and end points of the layer

	//Determine material of current layer (and hence track length)
	if (start){ //If input value is the start point
		if (edgeVal >= *(layerEnds) && edgeVal <= *(layerEnds) + 2) trackLen = 4.0; //Particle in lead
		else trackLen = *(layerEnds + 1) - edgeVal; //How much scintilator the particle will propogate through
	}
	
	//If input value is the end
	if (!start){ 
		if (edgeVal >= *(layerEnds) && edgeVal <= *(layerEnds) + 2) trackLen = 0.0; //Particle in lead
		else trackLen = *(layerEnds + 0) + edgeVal; //How much scintilator the particle will propogate through
	}
	
	return trackLen;	
}

double trackLen_scint(double startPoint, double endPoint){
	double *startLayer = layerTerminus(startPoint); //Get first Layer information
	double *endLayer = layerTerminus(endPoint); //Get Ending Layer information
	double trackLen;

	//Get track length through first layer (likely not the full layer)
	trackLen = layerTrackLen_scint(startPoint);

	//Get track length through middle layers
	trackLen += (*(endLayer + 0) - *(startLayer + 1))/6 * 4; //Scintilate through all the middle layers

	//Get track length through final layer (likely not the full layer)
	trackLen += layerTrackLen_scint(endPoint, false);

	return trackLen;
}

double layerTrackLen_pb(double edgeVal, bool start){
	double trackLen;
	double *layerEnds = layerTerminus(edgeVal); //Get the start and end points of the layer

	//Determine material of current layer (and hence track length)
	if (start){ //If input value is the start point
		if (edgeVal >= *(layerEnds) && edgeVal <= *(layerEnds) + 2) trackLen = *(layerEnds) + 2 - edgeVal; //Particle in lead
		else trackLen = 0; //Particle in Scintilator
	}
	
	//If input value is the end
	if (!start){ 
		if (edgeVal >= *(layerEnds) && edgeVal <= *(layerEnds) + 2) trackLen = edgeVal - *(layerEnds); //Particle in lead
		else trackLen = 2; //Particle in Scintilator
	}
	
	
	return trackLen;	
}

double trackLen_pb(double startPoint, double endPoint){
	double *startLayer = layerTerminus(startPoint); //Get first Layer information
	double *endLayer = layerTerminus(endPoint); //Get Ending Layer information
	double trackLen;

	//Get track length through first layer (likely not the full layer)
	trackLen = layerTrackLen_pb(startPoint);

	//Get track length through middle layers
	trackLen += (*(endLayer + 0) - *(startLayer + 1))/6 * 2; //Scintilate through all the middle layers

	//Get track length through final layer (likely not the full layer)
	trackLen += layerTrackLen_pb(endPoint, false);

	return trackLen;
}