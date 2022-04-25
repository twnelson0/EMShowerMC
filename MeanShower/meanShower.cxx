//Standard c++ libaries
#include <iostream>
#include <vector>
#include <fstream>
#include <thread>
#include <utility>

//ROOT Libaries
#include "TCanvas.h"
#include "TRandom3.h" 
#include "TH1.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TFile.h"
#include "TString.h"

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
void showerAction1d(showerR2 &inShower, double E_Crit, double crntRadLen, bool verbose = false){
	int inCount = inShower.showerSize();
	if (verbose) std::cout << "In coming Particles = " << inCount << std::endl;
	int showPart = 0; //Count the number of particles undergoing radiative processes and pair production in a generation
	
	for (int i = 0; i < inCount; i++){
		double partTheta,partP;
		double incEnergy = inShower.EVec.at(i - showPart);

		/*if (verbose){
			std::cout << "Shower Dump #" << i << std::endl;
			inShower.showerDump();
		}*/

		if (verbose) std::cout << "There are " << inShower.showerSize() << " particles in the shower" << std::endl;

		//Photon Interactions
		if(inShower.idVec.at(i - showPart) == 22){
			if (incEnergy>= 2*m_e && ceil(crntRadLen) == floor(crntRadLen)){ //Pair production
				if (verbose) std::cout << "Pair production" << std::endl;

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

				inShower.clearParticle(i - showPart); 
				showPart+=1; //Increment number of incident showering particles by 1

			}else if (incEnergy < 2*m_e){ if (verbose) {std::cout << "No longer pair producing" << std::endl;} else continue;}//inShower.printPart(i - showPart);}}
		}

		//Lepton Interactions
		else if (inShower.idVec.at(i - showPart) == 11 || inShower.idVec.at(i - showPart) == -11){
			if (incEnergy >= E_Crit && ceil(crntRadLen) == floor(crntRadLen)){ //Bremsstralung
				if (verbose) std::cout << "Lepton Interaction" << std::endl;
				//Lepton
				inShower.idVec.push_back(inShower.idVec.at(i - showPart));
				inShower.EVec.push_back(incEnergy*0.5);
				inShower.pVec.push_back(sqrt(pow(incEnergy*0.5,2) - pow(m_e,2)));
				inShower.thetaVec.push_back(0);

				//Photon
				inShower.idVec.push_back(22);
				inShower.EVec.push_back(incEnergy*0.5);
				inShower.pVec.push_back(incEnergy*0.5);
				inShower.thetaVec.push_back(0);

				inShower.clearParticle(i - showPart); 
				showPart+=1; //Increment number of incident showering particles by 1
			}

			else if (incEnergy < E_Crit){  //This may be removed
				if (verbose){
					std::cout << "No longer Bremsstrahlunging" << std::endl;
					//inShower.printPart(i - showPart);
				}
				//Shower Attenuation
				inShower.clearParticle(i - showPart); //Remove non radiative leptons
				showPart+=1; //Increment number of showerin gparticles by one
				
				//continue;
				//if (ceil(crntRadLen + dt) == floor(crntRadLen + dt)) std::cout << "if else logic broken" << std::endl;

				//Simulate Ionization loss
				//inShower.EVec.at(i) = inShower.EVec.at(i) - ionizationLoss(inShower.EVec.at(i),radLen2Long(crntRadLen),radLen2Long(crntRadLen + dt));
				//inShower.pVec.at(i) = sqrt(pow(inShower.EVec.at(i),2) - pow(m_e,2));		
			}
		}

		else{
			if (verbose){std::cout << "Input Not Recognized" << std::endl; inShower.printPart(i);}
			else continue;
		}
	}
}

//Shower for full LLP event consisting of 2 particles
//void showerAction_Full(std::vector<>, double E_Crit, double crntRadLen, bool verbose = false){
	
//}


//Determine the physical location of each interaction given 
std::vector<double> getIntrPoints(int NX0, double size){
	std::vector<double> interPoints;
	double X0 = size/(double) NX0;

	for (int n = 1; n < NX0; n++){
		interPoints.push_back(n*X0);
	}

	return interPoints;
}


//Shower and record total number of scintilation photons produced 
std::vector<int> scintShower(int maxLen, double Einit, int partId, double critVal = 5){
	std::vector<int> showerScint; //Number of scintilation photons
	particleR2 initPart = particleR2(Einit,0,partId); //Set up initial particle
	showerR2 inShower = showerR2(initPart); //Set up shower object
	//int trackArr[2] = 

	//Shower through the calorimter
	for (int i = 0; i < maxLen; i++){
		int nChargeTrack = inShower.chargedTracks(); //Count the number of charged tracks
		showerScint.push_back(16000*layerTrackLen_scint(radLen2Long(i),radLen2Long(i+1))*nChargeTrack*0.12*0.15); //Get scintilation photons detected 
		showerAction1d(inShower,5,(double) i, false); //Shower after 1 radiation length
		if (nChargeTrack == 0) break;
	}

	return showerScint;
}


//Shower and record total number of scintilation photons produced (Multithreaded)
void scintShower_Thread(std::vector<int> &showerScint, int maxLen, double Einit, int partId, double critVal = 5.){
	//std::vector<int> showerScint; //Number of scintilation photons
	std::cout << "Starting Single Shower" << std::endl;
	particleR2 initPart = particleR2(Einit,0,partId); //Set up initial particle
	showerR2 inShower = showerR2(initPart); //Set up shower object
	//int trackArr[2] = 

	//Shower through the calorimter
	for (int i = 0; i < maxLen; i++){
		int nChargeTrack = inShower.chargedTracks(); //Count the number of charged tracks
		if (nChargeTrack == 0) break; //End shower if there are no more charged particles
		//std::cout << "Track Length between " << i << " and " << ++i << " = " << layerTrackLen_scint(radLen2Long(i),radLen2Long(i+1)) << std::endl;
		showerScint.push_back(16000*layerTrackLen_scint(radLen2Long(i),radLen2Long(i+1))*nChargeTrack*0.12*0.15); //Get scintilation photons detected 
		showerAction1d(inShower,5,(double) i, false); //Shower after 1 radiation length
	}

	//return showerScint;
}


int main(){
	//std::cout << "Test" << std::endl;
	//Test a very basic mean showering model of a 50 GeV positron in 1d going through 24 radiation lengths
	double E0 = 500e3; //Starting energy in MeV
	//double E0 = 8;
	particleR2 initPart = particleR2(E0,0,11);
	showerR2 inShower = showerR2(initPart);
	int Nitr = 500;
	int photoArr[25];
	//int layerArr[66];
	//double eLossArr[25];
	int genArr[25];

	//Scintilation photons
	std::vector<double> timeStampVec, leptVec;
	std::vector<double> inputE = linspace(500,5000,10);
	std::vector<double> sumScintPhoto;

	//long testPhoton = scintShower(25,E0,11);
	//std::cout << testPhoton << std::endl;

	//ROOT Objects
	TFile *f1 = new TFile("ScintPhotoOut1.root","RECREATE");
	//TCanvas *c1 = new TCanvas("c1","c1",500,500);
	//TRandom3 *randGen = new TRandom3();

	//Add vector of energies to the ROOT file
	f1->WriteObject(&inputE,"Initial_E");
	int showerNum = 0;

	//Loop over all incident particle energies
	for ( double E : inputE) {
		TString crntName; crntName.Form("Shower_%d",showerNum);
		std::cout << "E = " << E << " GeV" << std::endl;
		//Create 2 threads for each input particle
		//std::thread t1,t2;
		double photoSum = 0;
		std::vector<int> showerVec_elec, showerVec_pos, totalVec; 
		std::vector<int> showerVecArr[2];
		
		std::thread t1(scintShower_Thread,std::ref(showerVec_elec),25,E*1e3,11,5); //Electron Shower
		std::thread t2(scintShower_Thread,std::ref(showerVec_pos),25,E*1e3,-11,5); //Positron Shower

		t1.join();
		t2.join();

		//Single Threading
		//scintShower_Thread(showerVec_elec,25,E,11,5);
		//scintShower_Thread(showerVec_pos,25,E,-11,5);

		//Update the shower vector array with the first index benig the largest vector
		if (showerVec_elec.size() > showerVec_pos.size()){
			showerVecArr[0] = showerVec_elec;
			showerVecArr[1] = showerVec_pos;
		}else if (showerVec_elec.size() < showerVec_pos.size()){
			showerVecArr[0] = showerVec_pos;
			showerVecArr[1] = showerVec_elec;
		}else if (showerVec_elec.size() == showerVec_pos.size()){
			showerVecArr[0] = showerVec_elec;
			showerVecArr[1] = showerVec_pos;
		}
		
		//if (showerVec_pos.size() != showerVec_elec.size()) std::cout <<"size difference" << std::endl;

		//scintShower(25,E*1e3,11); //Old single threaded shower
		//Get total number of scintilation photons for a LLP event
		totalVec = showerVecArr[0];
		for (int i = 0; i < showerVecArr[1].size(); i++){totalVec.at(i) = totalVec.at(i) + showerVecArr[1].at(i);}
		//for (int i = 0; i < showerVec_elec.size(); i++){totalVec.push_back(showerVec_elec.at(i) + showerVec_pos.at(i));}
		f1->WriteObject(&totalVec,crntName); //Write current Vector
		showerNum++;
		//for (int phot : totalVec){std::cout << phot << std::endl;}

		//Combine the scintilations photons into 1 vector
		/*for (int i = 0; i < showerVec.size(); i++){photoSum += (double)showerVec.at(i);}
		sumScintPhoto.push_back(photoSum);*/
	} //Fill Scintilaiton photon vector
	

	f1->Close();
	//delete f1;
	
	//Loop over the total scintilation photon vector
	//for (double v : sumScintPhoto){std::cout << v << std::endl;}

	//Graph SCintilation Photons as Sanity Check
	/*TGraph *g1 = new TGraph(sumScintPhoto.size(),&(inputE[0]),&(sumScintPhoto[0]));
	g1->GetYaxis()->SetTitle("Input Energy (GeV)");
	g1->GetYaxis()->SetTitle("Total Number of Scintilation Photons");
	g1->SetTitle("Total Scintilation Photons");
	g1->Draw("A*");
	c1->SaveAs("Scintilation_PhotonTest_2.pdf");

	delete g1;*/
	//c1->Close();


	//Write Scintilation Photon Count and energies out to a csv file
	/*std::ofstream scintCSV;
	scintCSV.open("TestCSV.csv");
	scintCSV << "Energy (MeV), Scintilation Photons" << "\n";
	for (int i = 0; i < scintPhotoVec.size();i++){
		scintCSV << inputE.at(i) << "," << scintPhoto.at(i) << "\n";
	}

	scintCSV.close();*/



	//Debugging number of leptons
	//int leptonNumber[25];
	//int partArr[25];
	//int genArray[25];

	//Set up photon array
	//for (int i = 0; i < 25; i++){photoArr[i] = 0; genArr[i] = i;}
	
	//Scintilation photon Histogram
	

	//std::vector<particleR2> showerVec;
	//showerVec.push_back(initPart);

	//int startLeptonNum = leptonNumber(inShower);
	//int startLeptonNum = inShower.leptonNumber();
	//double dt = 1/ (double) Nitr; //Continous step size in radiation lengths

	//Propogate over 25 Radiation Lengths
	/*for (int t = 0; t < 25; t++){ //Loop over the generations
		std::cout << "Generation " << t << std::endl;
		//genArray[t] = t;

		//Count the number of charged tracks at the begining of each generation
		int nChargeTrack = inShower.chargedTracks();
		//leptonNumber[t] = nChargeTrack;
		//partArr[t] = inShower.showerSize();
		//inShower.showerDump();
		//if (t == 1) {showerAction1d(inShower, 5, (double) t, true); inShower.showerDump();}
		//else showerAction1d(inShower, 5, (double) t);
		showerAction1d(inShower, 5, (double) t, false);
		if (nChargeTrack == 0) {std::cout << "Ending Shower" << std::endl; break;} //End shower

		//Check energy conservation 
		//if (Ecrt != E0)std::cout << "!Energy is Not being conserved!\n" << Ecrt << " != " << E0 << std::endl;
		

		double nGamma = 0; //Number of Scintilation photons
		nGamma = 16000*layerTrackLen_scint(radLen2Long(t),radLen2Long(t+1))*nChargeTrack;
		//if (t == 0) timeStampVec.push_back(radLen2Time(1,E0/pow(2,t),m_e)); //First track
		//else timeStampVec.push_back(timeStampVec.at(t - 1) + radLen2Time(1,E0/pow(2,t),m_e));
		leptVec.push_back((double) nChargeTrack);
		timeStampVec.push_back((double) t);
		scintPhotoVec.push_back(nGamma*0.15*0.12); //Scale by both PMT quantum efficency and 15% fiber transmission (12% for PMT place holder from LHCB paper)
		//std::cout << nGamma << std::endl;

		std::cout << "There are " << inShower.showerSize() << " particles in the shower" << std::endl;
		//std::cout << "There where " << nChargeTrack << " charged tracks at the start" << std::endl;
		//if (t == 0) inShower.showerDump();
		//nGamma = (long) layerTrackLen_scint(radLen2Long((double)t),radLen2Long((double)(t + 1)))*chargedTrackNum*8000; //This is causing integer overflows these numbers are too big
		//std::cout << "Number of scintilation photons = " << nGamma << std::endl;
	}*/

	//Graph the number of charged tracks and the number of particles
	/*TGraph *gLep = new TGraph(25,genArr,leptonNumber); gLep->SetMarkerStyle(21); gLep->SetTitle("Charged Tracks");
	TGraph *gAll = new TGraph(25,genArr,partArr); gAll->SetMarkerStyle(22); gAll->SetTitle("All Particles");
	TLegend *l1 = new TLegend();
	TMultiGraph *partGraph = new TMultiGraph();
	partGraph->Add(gLep);
	partGraph->Add(gAll);
	partGraph->GetXaxis()->SetTitle("Depth X_{0}");
	partGraph->GetYaxis()->SetTitle("Number of Particles");
	partGraph->SetTitle("Charged Tracks Plot");
	partGraph->Draw("AP");
	c1->BuildLegend();
	c1->SaveAs("Shower_Particle_Plot.pdf");
	delete gLep;
	delete gAll;
	c1->Clear();*/


	//Graph the number of Scintilation Photons
	/*TGraph *gPhoto = new TGraph(scintPhotoVec.size(),&(timeStampVec[0]),&(scintPhotoVec[0]));
	gPhoto->GetXaxis()->SetTitle("Generation");
	gPhoto->GetYaxis()->SetTitle("Number of Scintilation Photons");
	gPhoto->Draw("AL*");
	c1->SaveAs("ScintPhotoPlot_1.pdf");
	c1->Clear();
	delete gPhoto;

	//Graph the number of charged tracks per generation
	TGraph *gCharge = new TGraph(scintPhotoVec.size(),&(timeStampVec[0]),&(leptVec[0]));
	gCharge->GetXaxis()->SetTitle("Depth X_{0}");
	gCharge->GetYaxis()->SetTitle("Number of ChargedTracks");
	gCharge->SetTitle("Nubmer of Charged Tracks");
	gCharge->Draw("AL*");
	c1->SaveAs("ChargedTrackPlot.pdf");*/

	//Graph the number of scintilation photons
	/*TGraph *g2 = new TGraph(25,genArr,photoArr);
	g2->GetXaxis()->SetTitle("Depth X_{0}");
	g2->GetYaxis()->SetTitle("Number of Scintilation Photons");
	g2->SetTitle("");*/

	//c1->SaveAs("Electron_Positron_Number.pdf");

	//Memory managment 
	//delete randGen;
	//delete g1;
	//delete partGraph;
	//delete gCharge;
	//delete c1;
}



/*
I don't understand how to get energy deposited from this simple simulation, I get that the fractional energy deposited will increase but early on since I am halving the particles isn't 
no energy going to be deposited until the charged tracks start scintillating?
I am also a bit confused about when scintilation starts

CUrrently confused on the followng:
*How to simulate scintilation photons
*How to get energy depoosited per distance traveled (the two should be related)
*Do I need to consider energy loss that beyond mean bremstrahulng when simulating early part of the shower
*Attenuation of the shower towards the end No idea where to even start with this
*/
