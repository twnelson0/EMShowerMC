#include <iostream>
#include <vector>
#include <map>
#include "MatterCalc.h"
//#include <math.h>

using namespace physConstants;
/*class particle{
public:
	particle(double E, double m, double ){}

private:
	//Variables
	double Epart, mPart;
	bool fermi;
}*/

/*
This class attempts to store the information about the particles undergoing decay and EM shower
it assumes a 2d geometry and so the angle theta is with respect to the x-axis

*/

//Generic 2-d particle class with an angle from x-axis and id number determined by PGD
class particleR2{
public:
	particleR2(double E, double theta, int idNum){ //Energy must be in MeV
		double pMag; //Magnitude of 3-momenta
		p0 = E; //Store energy
		partID = idNum; //Store Particle ID

		//Determine type of particle
		if (abs(idNum) == 11) m = m_e;
		if (idNum == 22) m = 0;
		if (idNum == 26) m = 50e3; //LLP Placeholder

		//Obtain the 3-momeneta
		pMag = sqrt(pow(E,2) - (m,2));
		p1 = pMag*cos(theta);
		p2 = pMag*sin(theta);
	};
	//~particleR2(); Destructor not needed?

	//Get particle ID
	int id() {return partID;}

	//Kinematic information
	double E() {return p0;}



private:
	double p0,p1,p2,m;
	int partID;
	
};

//2d Particle Class (with inheritance)
/*class particleR2{
public:
	particleR2(double E, double theta, bool antiPart = false);
	int idNum() {return partID;}

	//particleR2();

private:
	int partID;
	double px, py;
};

//Electron Class
class elecR2 : public particleR2{
public:
	using particleR2::particleR2;
	this->setIDNum();

	//int idNum(){return partID;}

private:
	//Calculate Momenta given the angle
	void setIDNum() {partID = (antiPart) ? -11 : 11;}

};

//Photon Class
class photonR2 : public particleR2{
public:
	//photonR2(){partID = 22;}
	this->setIDNum();
	//int idNum(){return partID;}

private:
	//Calculate momenta given the angle
	void setIDNum() {partID = 22;}
};*/

//LLP Class
/*class LLPR2 : public particleR2{

};*/