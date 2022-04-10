//Standard c++ libaries
#include <iostream>
#include <vector>

//Root and other Libaries
#include "TRandom3.h" //*NO Idea if I'm going to use ROOT libaries*
#include "../MathMethods/MatterCalc.h"
#include "../MathMethods/particles.h"
#include "../CalModel/CalGeo.h"

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

//Get status of particle shower for debugging purposes (Depercated)
void getStatus(std::vector<particleR2> shower){
	std::cout << "!!SHOWER STATUS!!" << std::endl;
	std::cout << "There are " << shower.size() << " particles" << std::endl;

	for (particleR2 & part : shower){
		std::cout << part.id() << std::endl;
		std::cout << part.E() << std::endl;
	}
}

//Check that lepton number is being conserved (Depercated)
int leptonNumber(std::vector<particleR2> shower){
	int leptonNumber = 0;
	for (particleR2 & part : shower){
		if (abs(part.id()) == 11) leptonNumber += part.id();
		else continue;
	}

	return leptonNumber/11;
}

//Paticle interaction in 1d WILL CAUSE SEG FAULTS!! (Deperacted)
/*void showerAction1d_Dep(std::vector<particleR2>& crntGen, double E_Crit){ //For a realastic shower of 24 particles this vector will contain 2^24 elements, this does not make c++ happy
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
			std::cout << "Input Not Recognized" << std::endl;
		}

		incCount++;
	}

	crntGen.erase(crntGen.begin(),crntGen.begin() + incCount); //This is not functioning correctly/how I expect it to

	//std::cout << "There are " << incCount << " incident particles" << std::endl;
}*/

//Obtain Energy loss due to ionizaiton assuming MIP over a given range 
double ionizationLoss(double E0, double startVal, double endVal){
	double ELoss;

	//Obatin the initial track loss 
	//Obatin the Initial Scintilating track length

	//Obtain the Initial lead track length

	//Obtain the intervening tracks

	//Obtain the final lead track loss

	//Obtain the final Scintilating track loss


}

//Obatin number of scintilation photons produced over a given length of detector
double scintPhoton(double E0, double startPoint, double endPoint){

	//In MIP Photon Yield = 16000/cm

}

//Second version of the 1 dimeinsional showering function with hopefully less memory problems
void showerAction1d(showerR2 &inShower, double E_Crit, TRandom3 *gen){
	int inCount = inShower.showerSize();
	int showPart = 0; //Count the number of particles undergoing radiative processes and pair production
	
	for (int i = 0; i < inCount; i++){
		double partTheta,partP;
		double incEnergy = inShower.EVec.at(i);
		double splitFrac = gen->Uniform(1);

		//Photon Interactions
		if (inShower.idVec.at(i) == 22 && incEnergy>= 2*m_e){ //Pair production
			//Electron
			inShower.idVec.push_back(11);
			inShower.EVec.push_back(incEnergy*splitFrac);
			inShower.pVec.push_back(sqrt(pow(incEnergy*splitFrac,2) - pow(m_e,2)));
			inShower.thetaVec.push_back(0);

			//Positron
			inShower.idVec.push_back(-11);
			inShower.EVec.push_back(incEnergy*(1-splitFrac));
			inShower.pVec.push_back(sqrt(pow(incEnergy*(1-splitFrac),2) - pow(m_e,2)));
			inShower.thetaVec.push_back(0);

			showPart+=1; //Increment number of incident showering particles by 1
			inShower.clearParticle(i);

		}else if (incEnergy < 2*m_e){std::cout << "No longer pair producing" << std::endl;}

		//Lepton Interactions
		else if (abs(inShower.idVec.at(i)) == 11 && incEnergy >= E_Crit){ //Bremsstralung
			//Lepton
			inShower.idVec.push_back(inShower.idVec.at(i));
			inShower.EVec.push_back(incEnergy*splitFrac);
			inShower.pVec.push_back(sqrt(pow(incEnergy*splitFrac,2) - pow(m_e,2)));
			inShower.thetaVec.push_back(0);

			//Photon
			inShower.idVec.push_back(22);
			inShower.EVec.push_back(incEnergy*(1-splitFrac));
			inShower.pVec.push_back(incEnergy*(1-splitFrac));
			inShower.thetaVec.push_back(0);

			showPart+=1; //Increment number of incident showering particle
			inShower.clearParticle(i);

		}else if (incEnergy < E_Crit){ 
			std::cout << "No longer Bremsstrahlunging" << std::endl;

			//Simulate Ionization loss
			
		}

		else{
			std::cout << "Input Not Recognized" << std::endl;
			inShower.printPart(i);
		}
	}
}


//Determine the physical location of each interaction given 
std::vector<double> getIntrPoints(int NX0, double size){
	std::vector<double> interPoints;
	double X0 = size/(double) NX0;

	for (int n = 1; n < NX0; n++){
		interPoints.push_back(n*X0);
	}

	return interPoints;
}


int main(){
	//std::cout << "Test" << std::endl;
	//Test a very basic mean showering model of a 50 GeV positron in 1d going through 24 radiation lengths
	double E0 = 50e3;
	particleR2 initPart = particleR2(E0,0,-11);
	showerR2 inShower = showerR2(initPart);
	TRandom3 *randGen = new TRandom3();
	//std::vector<particleR2> showerVec;
	//showerVec.push_back(initPart);

	//int startLeptonNum = leptonNumber(inShower);
	int startLeptonNum = inShower.leptonNumber();

	//Propogate over 25 Radiation Lengths
	for (int t = 0; t < 25; t++){
		showerAction1d(inShower,5, randGen);
		std::cout << "There are " << inShower.showerSize() << " particles in the shower" << std::endl;

		//Check Lepton Number conversion
		if (startLeptonNum != inShower.leptonNumber()) std::cout << "Lepton Number Not Being Conserved" << std::endl;
		double Ecrt = 0;

		for (int i = 0; i < inShower.showerSize(); i++){
			Ecrt += inShower.EVec.at(i);
		}

		//Check energy conservation 
		//if (Ecrt != E0){std::cout << "!Energy is Not being conserved!" << Ecrt << " != " << E0 << std::endl;}
	}



	//Memory managment 
	delete randGen;
}


/*
Bug List:
 

Other Stuff:
Running over 24 layers and c++ is not happy a 2^24ish element vector makes c++ seg fault, I may need a better way to store the data on each shower (perhaps an LHE File??)
Shower information needs to be stored in an external file or shower information should not be stored but rather useful information should be calucalted from information about the shower
Need to deal with attenuation of particles otherwise c++ will segfault and die

*/