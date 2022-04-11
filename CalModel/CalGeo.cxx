#include "CalGeo.h"
#include <iostream>
#include <cmath>

int ind(double x, double lowBound, double upBound){
	if(x >= lowBound && x <= upBound) return 1;
	else return 0;
}

double radLen2Long(double t){return (39.6/25)*t;}

double long2RadLen(double z){return (25/39.6)*z;}

int crntLayer(double crntPoint){return (int) (floor(10*crntPoint) - (int) floor(10*crntPoint) % 4)/6;}

double* layerTerminus(double startPoint){
	static double layerArr[2];
	double startFloor, layerStart, layerEnd;
	int layerMod;

	//Determine end points of layer
	startFloor = std::floor(startPoint*10); //Get floor of start point
	layerMod = (int) startFloor % 6; //Get modulus 6 of start point (roughly)

	layerArr[0] = (startFloor - (double) layerMod)/10; //Get start point of layer
	layerArr[1] = (startFloor + (double) (6 - layerMod))/10; //Get end point of layer

	return layerArr;
}

double layerTrackLen_scint(double edgeVal, bool start){
	double trackLen;
	double *layerEnds = layerTerminus(edgeVal); //Get the start and end points of the layer

	//Determine material of current layer (and hence track length)
	if (start){ //If input value is the start point
		if (edgeVal >= *(layerEnds) && edgeVal <= *(layerEnds) + 0.2) trackLen = 0.4; //Particle in lead
		else trackLen = *(layerEnds + 1) - edgeVal; //How much scintilator the particle will propogate through
	}
	
	//If input value is the end
	if (!start){ 
		if (edgeVal >= *(layerEnds) && edgeVal <= *(layerEnds) + 0.2) trackLen = 0.0; //Particle in lead
		else trackLen = edgeVal - *(layerEnds) - 0.2; //How much scintilator the particle will propogate through
	}
	
	return trackLen;	
}

double trackLen_scint(double startPoint, double endPoint){
	double *startLayer = layerTerminus(startPoint); //Get first Layer information
	double *endLayer = layerTerminus(endPoint); //Get Ending Layer information
	double trackLen = 0;

	//Transversal through part of a single layer
	if (endPoint < *(startLayer + 1)){
		double effStart, effEnd; //Start and end points w.r.t to single layer
		effStart = (startPoint - (floor(startPoint) - *(startLayer)))*10;
		effEnd = (endPoint - (floor(endPoint) - *(startLayer)))*10;
		//Determine how much track is in the scintilator
		trackLen += ((6 - effEnd*ind(effEnd,2,6)) - (effStart*ind(effStart,2,effEnd) - 2))*0.1;

	}else{ //Transveral through 1 or more layers
		//Get track length through first layer (likely not the full layer)
		trackLen += layerTrackLen_scint(startPoint);

		//Get track length through middle layers
		trackLen += (*(endLayer + 0) - *(startLayer + 1))/0.6 * 0.4; //Scintilate track length through all the middle layers

		//Get track length through final layer (likely not the full layer)
		if (endPoint > *(startLayer + 1)) trackLen += layerTrackLen_scint(endPoint, false);

	}

	return trackLen;
}

double layerTrackLen_pb(double edgeVal, bool start){
	double trackLen;
	double *layerEnds = layerTerminus(edgeVal); //Get the start and end points of the layer

	//Determine material of current layer (and hence track length)
	if (start){ //If input value is the start point
		if (edgeVal >= *(layerEnds) && edgeVal <= *(layerEnds) + 0.2) trackLen = *(layerEnds) + 0.2 - edgeVal; //Particle in lead
		else trackLen = 0; //Particle in Scintilator
	}
	
	//If input value is the end
	if (!start){ 
		if (edgeVal >= *(layerEnds) && edgeVal <= *(layerEnds) + 0.2) trackLen = edgeVal - *(layerEnds); //Particle in lead
		else trackLen = 0.2; //Particle in Scintilator
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
	trackLen += (*(endLayer + 0) - *(startLayer + 1))/0.6 * 0.2; //Scintilate through all the middle layers

	//Get track length through final layer (likely not the full layer)
	trackLen += layerTrackLen_pb(endPoint, false);

	return trackLen;
}