#include "CalGeo.h"
#include <iostream>
#include <cmath>

int ind(double x, double lowBound, double upBound){
	if(x >= lowBound && x <= upBound) return 1;
	else return 0;
}

int indSoftUppr(double x, double lowBound, double upBound){
	if(x >= lowBound && x < upBound) return 1;
	else return 0;
}

double radLen2Long(double t){return (39.6/25)*t;}

double long2RadLen(double z){return (25/39.6)*z;}

double radLen2Time(double tInter, double E, double mass){
	double gamma, beta;

	//Obatin the lorentz factor from the energy
	if (mass == 0){ //Photon
		return radLen2Long(tInter)*1e-2/3e8;
	}else{
		//Obtain beta for particle
		gamma = E/mass; 
		beta = sqrt(1 - pow(gamma,-2));
		return radLen2Long(tInter)*1e-2/(beta*3e8);
	}
}

int crntLayer(double crntPoint){return (int) (floor(10*crntPoint) - (int) floor(10*crntPoint) % 6)/6;}

double* layerTerminus(double startPoint){
	double *layerArr = new double[2];
	//double startFloor, layerStart, layerEnd;
	int layerNum = crntLayer(startPoint);
	//int layerMod;

	//Determine end points of layer
	/*startFloor = floor(startPoint*10); //Get floor of start point
	layerMod = ((int) startFloor) % 6; //Get modulus 6 of start point (roughly)

	layerArr[0] = (startFloor - (double) layerMod)/10; //Get start point of layer
	layerArr[1] = (startFloor + (double) (6 - layerMod))/10; //Get end point of layer*/
	*layerArr = 0.6*layerNum;
	*(layerArr + 1) = 0.6*(layerNum + 1);

	return layerArr;
}

double layerTrackLen_scint(double edgeVal, bool start){
	double trackLen;
	double *layerEnds = new double[2];
	layerEnds = layerTerminus(edgeVal); //Get the start and end points of the layer

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
	
	delete[] layerEnds;
	return trackLen;	
}

double trackLen_scint(double startPoint, double endPoint){
	double *startLayer = new double[2];
	double *endLayer = new double[2];
	startLayer = layerTerminus(startPoint); //Get first Layer information
	endLayer = layerTerminus(endPoint); //Get Ending Layer information
	double trackLen = 0;

	//Transversal through part of a single layer
	if (endPoint < *(startLayer + 1)){
		double effStart, effEnd; //Start and end points w.r.t to single layer
		effStart = (startPoint - *(startLayer))*10;
		effEnd = (endPoint - *(startLayer))*10;
		
		//Determine how much track is in the scintilator
		trackLen += (effEnd - (effStart*ind(effStart,2,6) + 2*indSoftUppr(effStart,0,2)))*ind(effEnd,2,6);//(effEnd*ind(effEnd,2,6) - effStart*ind(effStart,2,effEnd))*0.1;

		if (trackLen < 0) {
			std::cout << "!!Track Length is Negative" << std::endl;
			std::cout << effStart << "\n" << effEnd << std::endl;

		}

	}else{ //Transveral through 1 or more layers
		//Get track length through first layer (likely not the full layer)
		trackLen += layerTrackLen_scint(startPoint);

		//Get track length through middle layers
		trackLen += (*(endLayer + 0) - *(startLayer + 1))/0.6 * 0.4; //Scintilate track length through all the middle layers

		//Get track length through final layer (likely not the full layer)
		if (endPoint > *(startLayer + 1)) trackLen += layerTrackLen_scint(endPoint, false);

	}

	delete[] startLayer;
	delete[] endLayer;

	return trackLen;
}

double layerTrackLen_pb(double edgeVal, bool start){
	double trackLen;
	double *layerEnds = new double[2];
	layerEnds = layerTerminus(edgeVal); //Get the start and end points of the layer

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
	
	delete[] layerEnds;

	return trackLen;	
}

double trackLen_pb(double startPoint, double endPoint){
	double *startLayer = new double[2];
	double *endLayer = new double[2];
	startLayer = layerTerminus(startPoint); //Get first Layer information
	endLayer = layerTerminus(endPoint); //Get Ending Layer information
	double trackLen;

	//Get track length through first layer (likely not the full layer)
	trackLen = layerTrackLen_pb(startPoint);

	//Get track length through middle layers
	trackLen += (*(endLayer + 0) - *(startLayer + 1))/0.6 * 0.2; //Scintilate through all the middle layers

	//Get track length through final layer (likely not the full layer)
	trackLen += layerTrackLen_pb(endPoint, false);

	delete [] startLayer;
	delete [] endLayer;

	return trackLen;
}


void currentLayerMat(double crntPoint){
	int effCrntVal = (int) floor(crntPoint*10);

	if (effCrntVal >= 0 && effCrntVal <= 2) std::cout << "Lead" << std::endl;
	else std::cout << "Scintilator" << std::endl;

}