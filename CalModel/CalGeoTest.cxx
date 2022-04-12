#include <iostream>
#include "CalGeo.h"

int main(){
	std::cout << "Testing Geometry functions" << std::endl;
	double testArr[5] = {.3575,.1168,1.61789,1.96566,0}; //Test points
	//double *layer = new double[2];
	//layer = layerTerminus(0.4);
	//double *bullshit = layerTerminus(0.4);
	//double testArr[5] = {.2,.5,1.2,.1,0};

	//As starting points for a given layer
	/*for (int i = 0; i < 5; i++){
		std::cout << "Start Point = " << testArr[i] << std::endl;
		double scintLen = layerTrackLen_scint(testArr[i]);
		std::cout << "Scint Track Length = " <<  scintLen << std::endl;
		std::cout << "Lead Length = " << layerTrackLen_pb(testArr[i]) << std::endl;
	}
	std::cout << "\n\n" << std::endl;
	//As ending points for a given layer
	for (int i = 0; i < 5; i++){
		std::cout << "End Point = " << testArr[i] << std::endl;
		double scintLen = layerTrackLen_scint(testArr[i],false);
		std::cout << "Scint Track Length = " <<  scintLen << std::endl;
		std::cout << "Lead Length = " << layerTrackLen_pb(testArr[i],false) << std::endl;
	}*/

	std::cout << trackLen_scint(1/500,2/500) << std::endl;
	std::cout << trackLen_scint(0.1,0.4) << std::endl;
	std::cout << trackLen_scint(0.1,0.15) << std::endl;
	std::cout << trackLen_scint(0.4,0.6) << std::endl;
	//std::cout << crntLayer(0.4) << std::endl;
	//std::cout << *layer << " and " << *(layer + 1) << std::endl;

	//delete[] layer;
}
