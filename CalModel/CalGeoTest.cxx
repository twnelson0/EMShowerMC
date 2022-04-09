#include <iostream>
#include "CalGeo.h"

int main(){
	std::cout << "Testing Geometry functions" << std::endl;
	double testArr[5] = {3.575,1.168,16.1789,19.6566,0}; //Test points

	//As starting points for a given layer
	for (int i = 0; i < 5; i++){
		double scintLen = layerTrackLen_scint(testArr[i]);
		std::cout << "Scint Track Length = " <<  scintLen << std::endl;
		std::cout << "Lead Length = " << layerTrackLen_pb(testArr[i]) << std::endl;
	}
	std::cout << "\n\n" << std::endl;
	//As ending points for a given layer
	for (int i = 0; i < 5; i++){
		double scintLen = layerTrackLen_scint(testArr[i],false);
		std::cout << "Scint Track Length = " <<  scintLen << std::endl;
		std::cout << "Lead Length = " << layerTrackLen_pb(testArr[i],false) << std::endl;
	}


}
