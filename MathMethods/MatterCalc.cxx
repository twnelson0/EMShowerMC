#include <math.h>
#include <vector>
#include <iostream>
#include "MatterCalc.h"

using namespace physConstants;

//Radiation length
double radLen(uint Znuc, double Ar, double rho, bool densityBool = false){
	double Z = (double) Znuc; //Convert Z to a usable format
	double X0Inv = (densityBool) ? 4*alpha*pow(r0*1e2,2)*((rho*N_Av)/Ar)*Z*(Z + 1)*log(183./(pow(Z,1/3.))) : 4*alpha*pow(r0*1e2,2)*(N_Av/Ar)*Z*(Z + 1)*log(183./(pow(Z,1/3.)));

	return pow(X0Inv,-1);
}


//Efective Radiation length for double layered material
double effRadLen(std::vector<double> X0Vec, std::vector<double> lenVec){
	if (X0Vec.size() != lenVec.size()) std::cout << "Mismatched Vectors" <<std::endl;

	double X0Inv,totalLayer = 0;
	std::vector<double> fracVec;

	//Set up fractional length vector
	for (int i = 0; i < lenVec.size(); i++){totalLayer+=lenVec.at(i);}
	for (int i = 0; i < lenVec.size(); i++){fracVec.push_back(lenVec.at(i)/totalLayer);}


	//Obtain X0Inv
	for (int i = 0; i < X0Vec.size(); i++){X0Inv+= fracVec.at(i)/X0Vec.at(i);}

	return pow(X0Inv,-1);
}


//Linspace function
std::vector<double> linspace(double startVal, double endVal, uint size){
	std::vector<double> linspaveVec;
	double dx = (endVal - startVal)/(size - 1);

	//Fill vector
	for (uint i = 0; i < size; i++){linspaveVec.push_back(startVal + (double) i*dx);}

	return linspaveVec;
}

//Logspace fucntion
std::vector<double> logspace(double startVal, double endVal, uint size){
	std::vector<double> logspaceVec;
	double scale = exp(log(endVal/startVal)*pow((double) (size - 1), -1));

	for (uint i = 0; i < size; i++){logspaceVec.push_back(startVal*pow(scale,i));}

	return logspaceVec;
}


//Critical Energy of electron/positrons to bremsstralung
/*double critE(){

*/
