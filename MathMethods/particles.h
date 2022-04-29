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
	double p() {return sqrt(pow(p1,2) + pow(p2,2));}
	double theta() {return acos(p1/this->p());}

private:
	double p0,p1,p2,m;
	int partID;
	
};



//2d Shower object 
class showerR2{
public:
	showerR2(particleR2 incPart){
		idVec.push_back(incPart.id());
		EVec.push_back(incPart.E());
		locVec.push_back(0); //Push back particle location
		//thetaVec.push_back(incPart.theta());
		//pVec.push_back(incPart.p());
	}

	int showerSize(){return this->idVec.size();}

	int leptonNumber(){
		int leptonNumber = 0;
		for (int i = 0; i < this->showerSize(); i++){
			if (abs(this->idVec.at(i)) == 11) leptonNumber+=this->idVec.at(i);
			//else continue;
		}

		return leptonNumber/11;
	}

	//Count the number of charged tracks in the shower
	int chargedTracks(){
		int chargedTracks = 0;

		for (int i = 0; i < this->showerSize(); i++){
			if (abs(this->idVec.at(i)) == 11) chargedTracks+=1;
			//else continue;
		}

		return chargedTracks;
	}

	void printPart(int i){
		std::cout << "Particle = " << this->idVec.at(i) << std::endl;
		std::cout << "Energy = " << this->EVec.at(i) << std::endl;
		//std::cout << "Angle = " << this->thetaVec.at(i) << std::endl;
		//std::cout << "Momenta = " << this->pVec.at(i) << std::endl;
	}


	//Remove specific particle from shower structure 
	void clearParticle(int partIndx){
		this->idVec.erase(this->idVec.begin() + partIndx);
		this->EVec.erase(this->EVec.begin() + partIndx);
		this->locVec.erase(this->locVec.begin() + partIndx);
		//this->thetaVec.erase(this->thetaVec.begin() + partIndx);
		//this->pVec.erase(this->pVec.begin() + partIndx);
	}

	//Print all the particles
	void showerDump(){
		for (int i = 0; i < this->showerSize(); i++){
			this->printPart(i);
		}
	}


	//~showerR2();

	//Public varaibles
	std::vector<int> idVec;
	std::vector<double> EVec;
	std::vector<double> locVec;
	//std::vector<double> thetaVec;
	//std::vector<double> pVec;

/*private:
	std::vector<auto>[4] vectorAray; //Array of shower atribute vectors
*/
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