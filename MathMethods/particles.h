#include <iostream>
#include <vector>
#include <math.h>

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

//2d Particle Class
class particleR2{
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
};

//LLP Class
/*class LLPR2 : public particleR2{

};*/