//Standard c++ libaries
#include <iostream>
#include <vector>

//Root and other Libaries
//#include "TRandom3.h"
#include "../MathMethods/MatterCalc.h"
#include "../MathMethods/particles.h"

using namespace physConstants; //Physical constants from MatterCalc.h


/*
Outline
	* In this simplified model I will consider incident e-e+ particles on the calorimeter from the decay volume
	* Since these are the product of the decay of LLPs where in the decay volume they originate will be a exponentially random value with average of LLP lifetime
	* The angle of incidence of the LLP (and therefore the e+e-) will be a uniformly generated value within FASERs 2mrad acceptance
	* After each radiation length a given particle will lose half its energy
	* After 1 rad len a charged lepton will yield a photon and a charged lepton (Bremstralung)
	* After 1 rad len a photon will yield an electron positron pair (Pair production with nuclei really)
	* In principle pair production of photons should only happen every 9/7 Rad Len
	* In principle process for leptons stops at E_c and process stops for photons when E < 2m_e
	* When particles split I'm not sure how to handle direction, in principle if I consider a psudo 1d model all that matters is that energy splits in 2 mometa and location in detector irelevent
	* Negate more detailed interactions
	* In more detailed simulation when interaction occurs and resultent energies and directions should be more accuretly obtained from process cross section

Magnets will pull them apart

*/

//Particle Interaction
/*std::vector<?ParticelObj?> layerEdge(ParticleObj inPart){
	//Use angleular average and RMS info in Tav and SLAC paper 
	if (abs(inPart.idNum) == 11){ //Bremstralung

	}
	else if (inPart.idNum == 22){ //Pair production

	}
}*/

//Get status of particle shower for debugging purposes
void getStatus(std::vector<particleR2> shower){
	std::cout << "!!SHOWER STATUS!!" << std::endl;
	std::cout << "There are " << shower.size() << " particles" << std::endl;

	for (particleR2 & part : shower){
		std::cout << part.id() << std::endl;
		std::cout << part.E() << std::endl;
	}
}

//Check that lepton number is being conserved
int leptonNumber(std::vector<particleR2> shower){
	int leptonNumber = 0;
	for (particleR2 & part : shower){
		if (abs(part.id()) == 11) leptonNumber += part.id();
		else continue;
	}

	return leptonNumber/11;
}

//Paticle interaction in 1d
void showerAction1d(std::vector<particleR2>& crntGen, double E_Crit){ //For a realastic shower of 24 particles this vector will contain 2^24 elements, this does not make c++ happy
	uint incCount = 0;

	//Loop over all incident particles in the shower
	for (particleR2 & part : crntGen){ //I ultimently think that this structure is causing problems
		double incEnergy = part.E(); //Store incident energy of each particle

		//Photon Interactions
		if (part.id() == 22 && incEnergy >= 2*m_e){ //Pair Production
			//std::cout << "!Pair Prdoucing!" << std::endl;
			crntGen.push_back(particleR2(incEnergy/2,0,11)); 
			crntGen.push_back(particleR2(incEnergy/2,0,-11)); 
			
			//crntGen.erase(crntGen.begin());
		}else if (incEnergy < 2*m_e){std::cout << "No Longer pair producing" << std::endl;}
		
		//Electron/Poistron Interactions
		else if(abs(part.id()) == 11 && incEnergy >= E_Crit){ //Bremsstralung
			//std::cout << "!Bremsstralunging!" << std::endl;
			crntGen.push_back(particleR2(incEnergy/2,0,part.id()));  //This is only producing electrons or positrons >:(
			crntGen.push_back(particleR2(incEnergy/2,0,22)); 

			//crntGen.erase(crntGen.begin());
		}else if (incEnergy < E_Crit){std::cout << "No Longer Bremsstralunging" << std::endl;}

		else{
			std::cout << "Not Recognized" << std::endl;
		}

		incCount++;
	}

	crntGen.erase(crntGen.begin(),crntGen.begin() + incCount); //This is not functioning correctly/how I expect it to

	//std::cout << "There are " << incCount << " incident particles" << std::endl;
}


//Second version of the 1 dimeinsional showering function with hopefully less memory problems
void showerAction1d_2(showerR2 inShower, double E_Crit){
	int inCount = inShower.showerSize();
	
	for (int i = 0; i < inCount; i++){
		double partTheta,partP;
		double incEnergy = inShower.EVec.at(i);

		//Photon Interactions
		if (inShower.idVec.at(i) == 22 && incEnergy>= 2*m_e){ //Pair production
			//Electron
			inShower.idVec.push_back(11);
			inShower.EVec.push_back(incEnergy/2);
			inShower.pVec.push_back(sqrt(pow(incEnergy/2,2) - pow(m_e,2)));
			inShower.thetaVec.push_back(0);

			//Positron
			inShower.idVec.push_back(-11);
			inShower.EVec.push_back(incEnergy/2);
			inShower.pVec.push_back(sqrt(pow(incEnergy/2,2) - pow(m_e,2)));
			inShower.thetaVec.push_back(0);

		}else if (incEnergy < 2*m_e){std::cout << "No longer pair producing" << std::endl;}

		//Lepton Interactions
		if (abs(inShower.idVec.at(i)) == 11 && incEnergy >= E_Crit){ //Bremsstralung
			//Lepton
			inShower.idVec.push_back(inShower.idVec.at(i));
			inShower.EVec.push_back(incEnergy/2);
			inShower.pVec.push_back(sqrt(pow(incEnergy/2,2) - pow(m_e,2)));
			inShower.thetaVec.push_back(0);

			//Photon
			inShower.idVec.push_back(22);
			inShower.EVec.push_back(incEnergy/2);
			inShower.pVec.push_back(incEnergy/2);
			inShower.thetaVec.push_back(0);
		}
	}

	//Clear old particles
	inShower.idVec.erase(inShower.idVec.begin(), inShower.idVec.begin() + inCount);
	inShower.EVec.erase(inShower.EVec.begin(), inShower.EVec.begin() + inCount);
	inShower.thetaVec.erase(inShower.thetaVec.begin(), inShower.thetaVec.begin() + inCount);
	inShower.pVec.erase(inShower.pVec.begin(), inShower.pVec.begin() + inCount);
}

/*26/03/2022 22:39, I'm wondering if I need these particle objects in the E/2 splitting model they don't matter, it may matter for a more nuanced model but I'm wondering if there isn't
A more function appraoch to this problem

*/

/*
27/03/2022 16:04, I think I should start with a 1d simple mean showering program everything splits in half as it passes through the detector just to test some of the basics then I should add more nuanced behaviors
there is a lot of complexity present in showering and I don't want to waste my time figuring out the most ideal method of dealing with these sorts of interactions

27/03/2022 20:00, I think in the 1 d model I will have a vector of particles and on every iteration I'll loop over the vector shower that element and erase the first element of the vector leaving me with only
the current generation of the shower
*/


int main(){
	//std::cout << "Test" << std::endl;
	//Test a very basic mean showering model of a 50 GeV positron in 1d going through 24 radiation lengths
	double E0 = 50e3;
	particleR2 initPart = particleR2(E0,0,-11);
	std::vector<particleR2> showerVec;
	showerVec.push_back(initPart);

	int startLeptonNum = leptonNumber(showerVec);

	for (int t = 0; t < 15; t++){
		showerAction1d(showerVec,5); //Shower each bunch of particles
		std::cout << "There are " << showerVec.size() << " particles in the shower after " << t + 1 << " generations" << std::endl;

		//Check Lepton Number conversion
		if (startLeptonNum != leptonNumber(showerVec)) std::cout << "Lepton Number Not Being Conserved" << std::endl;
		double Ecrt = 0;

		for (particleR2 & part : showerVec){
			//std::cout << part.id() << std::endl;
			Ecrt += part.E();
		}

		//Check energy conservation 
		if (Ecrt != E0) {std::cout << "!Energy is Not being conserved!" << std::endl; getStatus(showerVec);}
	}
}


/*
Bug List:
 

Other Stuff:
Running over 24 layers and c++ is not happy a 2^24ish element vector makes c++ seg fault, I may need a better way to store the data on each shower (perhaps an LHE File??)
Shower information needs to be stored in an external file or shower information should not be stored but rather useful information should be calucalted from information about the shower
Need to deal with attenuation of particles otherwise c++ will segfault and die

*/