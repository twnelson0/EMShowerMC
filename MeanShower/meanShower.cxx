//Standard c++ libaries
#include <iostream>
#include <vector>

//ROOT Libaries
#include "TCanvas.h"
#include "TRandom3.h" 
#include "TH1.h"
#include "TGraph.h"
#include "TLegend.h"

//EM Shower Libaries
#include "../MathMethods/MatterCalc.h"
#include "../MathMethods/particles.h" //Shower object
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

//Obtain Energy loss due to ionizaiton assuming MIP over a given range 
double ionizationLoss(double E0, double startVal, double endVal){
	double ELoss;
	double *startLayer, *endLayer;

	//Obatin the start and end points of the first and final layers
	startLayer = layerTerminus(startVal);
	endLayer = layerTerminus(endVal);

	//Obatin the initial track loss 
	ELoss = rho_pb*2*layerTrackLen_pb(startVal); //Energy loss in Lead
	ELoss += rho_scint*2*layerTrackLen_scint(startVal); //Energy Loss due to Scintalting material

	//Obtain the intervening tracks
	for (int i = 0; i <  (int) (*(endLayer + 0) - *(startLayer + 1))/0.6; i++){
		ELoss += rho_pb*2*0.2; //Lead
		ELoss += rho_scint*2*0.4; //Scintilator
	}

	//Obatin the final layer track energy loss
	ELoss += rho_pb*2*layerTrackLen_pb(endVal,false);
	ELoss += rho_scint*2*layerTrackLen_scint(endVal,false);

	return ELoss;

}

//Second version of the 1 dimeinsional showering function with hopefully less memory problems
void showerAction1d(showerR2 &inShower, double E_Crit, double crntRadLen, double dt){
	int inCount = inShower.showerSize();
	int showPart = 0; //Count the number of particles undergoing radiative processes and pair production
	
	for (int i = 0; i < inCount; i++){
		double partTheta,partP;
		double incEnergy = inShower.EVec.at(i);

		//Photon Interactions
		if (inShower.idVec.at(i) == 22 && incEnergy>= 2*m_e && ceil(crntRadLen + dt) == floor(crntRadLen + dt)){ //Pair production
			//Electron
			inShower.idVec.push_back(11);
			inShower.EVec.push_back(incEnergy*0.5);
			inShower.pVec.push_back(sqrt(pow(incEnergy*0.5,2) - pow(m_e,2)));
			inShower.thetaVec.push_back(0);

			//Positron
			inShower.idVec.push_back(-11);
			inShower.EVec.push_back(incEnergy*0.5);
			inShower.pVec.push_back(sqrt(pow(incEnergy*0.5,2) - pow(m_e,2)));
			inShower.thetaVec.push_back(0);

			showPart+=1; //Increment number of incident showering particles by 1
			inShower.clearParticle(i);

		}else if (incEnergy < 2*m_e){continue;}//std::cout << "No longer pair producing" << std::endl;}

		//Lepton Interactions
		else if (abs(inShower.idVec.at(i)) == 11 && incEnergy >= E_Crit && ceil(crntRadLen + dt) == floor(crntRadLen + dt)){ //Bremsstralung
			//Lepton
			inShower.idVec.push_back(inShower.idVec.at(i));
			inShower.EVec.push_back(incEnergy*0.5);
			inShower.pVec.push_back(sqrt(pow(incEnergy*0.5,2) - pow(m_e,2)));
			inShower.thetaVec.push_back(0);

			//Photon
			inShower.idVec.push_back(22);
			inShower.EVec.push_back(incEnergy*0.5);
			inShower.pVec.push_back(incEnergy*0.5);
			inShower.thetaVec.push_back(0);

			showPart+=1; //Increment number of incident showering particle
			inShower.clearParticle(i);

		}else if (incEnergy < E_Crit){  //This may be removed
			//std::cout << "No longer Bremsstrahlunging" << std::endl;

			//Simulate Ionization loss
			inShower.EVec.at(i) = inShower.EVec.at(i) - ionizationLoss(inShower.EVec.at(i),radLen2Long(crntRadLen),radLen2Long(crntRadLen + dt));
			inShower.pVec.at(i) = sqrt(pow(inShower.EVec.at(i),2) - pow(m_e,2));		
		}

		else{
			continue;
			//std::cout << "Input Not Recognized" << std::endl;
			//inShower.printPart(i);
		}
	}
}


//Simulate ionization loss over a "continious" range

//Obatin number of scintilation photons produced over a given length of detector (assume inputs are in units of track length)
double scintPhoton(double E0, double startPoint, double endPoint){return 2*layerTrackLen_scint();}



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
	double E0 = 500e3;
	particleR2 initPart = particleR2(E0,0,11);
	showerR2 inShower = showerR2(initPart);
	TRandom3 *randGen = new TRandom3();
	int Nitr = 500;
	TCanvas *c1 = new TCanvas("c1","c1",500,500);
	int photoArr[26];
	int genArr[26];
	
	//Scintilation photon Histogram
	

	//std::vector<particleR2> showerVec;
	//showerVec.push_back(initPart);

	//int startLeptonNum = leptonNumber(inShower);
	int startLeptonNum = inShower.leptonNumber();
	double dt = 1/ (double) Nitr; //Continous step size in radiation lengths

	//Propogate over 25 Radiation Lengths
	for (int t = 0; t < 25; t++){ //Loop over the generations
		std::cout << "Generation " << t << std::endl;

		//Count the number of charged tracks at the begining of each generation
		int chargedTrackNum = inShower.chargedTracks();
		long nGamma; //Number of Scintilation photons

		for (int i = 0; i < Nitr; i++){ //Simulate behavior between "generations"
			showerAction1d(inShower, 5, (double) t + i*dt, dt);
			//std::cout << "There are " << inShower.showerSize() << " particles in the shower" << std::endl;

			//Check Lepton Number conversion
			if (startLeptonNum != inShower.leptonNumber()) std::cout << "Lepton Number Not Being Conserved" << std::endl;
			//double Ecrt = 0;

			/*for (int i = 0; i < inShower.showerSize(); i++){
				Ecrt += inShower.EVec.at(i);
			}

			//Check energy conservation 
			if (Ecrt != E0){std::cout << "!Energy is Not being conserved!" << Ecrt << " != " << E0 << std::endl;}*/
		} 

		//nGamma = (long) layerTrackLen_scint(radLen2Long((double)t),radLen2Long((double)(t + 1)))*chargedTrackNum*8000; //This is causing integer overflows these numbers are too big
		//std::cout << "Number of scintilation photons = " << nGamma << std::endl;
	}

	/*
		Could treat Bremstrahlung as a discrete process that occurs at the start or end of every loop iteration, and do continous simulations of the middle steps
	*/

	//c1->SaveAs("TestCanvas.pdf");

	//Memory managment 
	delete randGen;
	delete c1;
}

/*
I will wan to count the scintilation photons either during the continous loop or 
after each generation, taking into account the number of photons in each layer

If I use the first option I can write a function that take a given range and determines 
which scintiltaor it falls within (perhaps cutting out any part that isn't in a scintilator)

Still not sure what the best solution is as if I choose the first I have to wory about a range 
enetering into the adjacent lead layer, it almost seems like I need to include a function that cuts
off any part of the track length not contained in the scintlator before getting N-gamma

Scheme:
Input(start,end): function gets intial or final scintilation layer, not valid if |start - end| > 6 mm
From there can determine which scintilation layer I'm in

Want another function that can trim any part of a length outside the scintilator in similar to 
the layerTrackLen_scint function

from there I can map a number of photons produced within a given length of scintilator to a specific layer and store it
This should also allow me to deal with generation lines that fall within a scintilator


*/

