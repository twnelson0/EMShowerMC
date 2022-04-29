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

//Second version of the 1 dimeinsional showering function with hopefully less memory problems (include radiation length in cm)
void showerAction1d(showerR2 &inShower, double E_Crit, double radLen, TRandom3 *randGen, bool verbose = false){
	int inCount = inShower.showerSize();
	if (verbose) std::cout << "In coming Particles = " << inCount << std::endl;
	int showPart = 0; //Count the number of particles undergoing radiative processes and pair production in a generation
	
	for (int i = 0; i < inCount; i++){
		double partTheta,partP;
		double incEnergy = inShower.EVec.at(i - showPart);
		double crntLoc = inShower.locVec.at(i - showPart);
		double disp, newPos;

		/*if (verbose){
			std::cout << "Shower Dump #" << i << std::endl;
			inShower.showerDump();
		}*/

		if (verbose) std::cout << "There are " << inShower.showerSize() << " particles in the shower" << std::endl;
		if (verbose) std::cout << "i - showPart = " << i - showPart << std::endl;

		//Photon Interactions
		if(inShower.idVec.at(i - showPart) == 22){
			if (incEnergy>= 2*m_e){ //Pair production
				if (verbose) std::cout << "Pair production" << std::endl;
				disp = randGen->Exp(radLen*9/7);//Distance particle travels before undergoing an interaction
				newPos = crntLoc + disp; //Location of decay

				if (newPos >= radLen2Long(25)) {inShower.clearParticle(i - showPart); if (verbose) std::cout << "Outisde Cal" << std::endl;} //Particle outside cal attenuate it
				else{
					//Electron
					inShower.diffLocVec.push_back(disp);
					inShower.locVec.push_back(newPos); 
					inShower.idVec.push_back(11);
					inShower.EVec.push_back(incEnergy*0.5);
					//inShower.pVec.push_back(sqrt(pow(incEnergy*0.5,2) - pow(m_e,2)));
					//inShower.thetaVec.push_back(0);

					//Positron
					inShower.diffLocVec.push_back(disp);
					inShower.locVec.push_back(newPos);
					inShower.idVec.push_back(-11);
					inShower.EVec.push_back(incEnergy*0.5);
					//inShower.pVec.push_back(sqrt(pow(incEnergy*0.5,2) - pow(m_e,2)));
					//inShower.thetaVec.push_back(0);

					inShower.clearParticle(i - showPart); 
					//showPart+=1; //Increment number of incident showering particles by 1
				}
				showPart+=1; //Increment number of incident showering particles by 1

			}else if (incEnergy < 2*m_e){ if (verbose) {std::cout << "No longer pair producing" << std::endl;} else continue;}//inShower.printPart(i - showPart);}}
		}

		//Lepton Interactions
		else if (inShower.idVec.at(i - showPart) == 11 || inShower.idVec.at(i - showPart) == -11){
			if (incEnergy >= E_Crit){ //Bremsstralung
				disp = randGen->Exp(radLen);
				newPos = crntLoc + disp; //Location of decay

				if (newPos >= radLen2Long(25)) {inShower.clearParticle(i - showPart);if (verbose) std::cout << "Outisde Cal" << std::endl;} //Particle outside cal attenuate
				else{
					if (verbose) std::cout << "Lepton Interaction" << std::endl;
					//Lepton
					inShower.diffLocVec.push_back(disp);
					inShower.locVec.push_back(newPos);
					inShower.idVec.push_back(inShower.idVec.at(i - showPart));
					inShower.EVec.push_back(incEnergy*0.5);
					//inShower.pVec.push_back(sqrt(pow(incEnergy*0.5,2) - pow(m_e,2)));
					//inShower.thetaVec.push_back(0);

					//Photon
					inShower.diffLocVec.push_back(disp);
					inShower.locVec.push_back(newPos);
					inShower.idVec.push_back(22);
					inShower.EVec.push_back(incEnergy*0.5);
					//inShower.pVec.push_back(incEnergy*0.5);
					//inShower.thetaVec.push_back(0);

					inShower.clearParticle(i - showPart); 
					//showPart+=1; //Increment number of incident showering particles by 1
				}
				showPart+=1; //Increment number of incident showering particles by 1
			}

			//Range Bug is occuring !!here!!
			else if (incEnergy < E_Crit){  
				if (verbose){
					std::cout << "No longer Bremsstrahlunging" << std::endl;
					//inShower.printPart(i - showPart);
					std::cout << inShower.showerSize() << std::endl;
				}
				//Shower Attenuation
				inShower.clearParticle(i - showPart); //Remove non radiative leptons
				showPart+=1; //Increment number of showering particles by one

				if (verbose) std::cout << "Particles Cleared" << std::endl; 
				
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


//Ionization energy loss due to traversal
void ionELoss(showerR2 &inShower, double startRadLen, double endRadLen){
	int inCount = inShower.showerSize();
	for (int i = 0; i < inCount; i++){
		if (abs(inShower.idVec.at(i)) == 11){ //Ionization loss for leptons only
			inShower.EVec.at(i) -= 2.502*(endRadLen - startRadLen);//Homogenous only
		} 
	}
}

//Shower and record total number of scintilation photons produced 
/*std::vector<int> scintShower(int maxLen, double Einit, int partId, double critVal){
	std::vector<int> showerScint; //Number of scintilation photons
	particleR2 initPart = particleR2(Einit,0,partId); //Set up initial particle
	showerR2 inShower = showerR2(initPart); //Set up shower object
	//int trackArr[2] = 

	//Shower through the calorimter
	for (int i = 0; i < maxLen; i++){
		int nChargeTrack = inShower.chargedTracks(); //Count the number of charged tracks
		showerScint.push_back(16000*layerTrackLen_scint(radLen2Long(i),radLen2Long(i+1))*nChargeTrack*0.12*0.15); //Get scintilation photons detected 
		showerAction1d(inShower,5, X0, gen , false); //Shower after 1 radiation length
		if (nChargeTrack == 0) break;
	}

	return showerScint;
}*/


//Shower and record total number of scintilation photons produced (Multithreaded)
/*void scintShower_Thread_Scint(std::vector<int> &showerScint, int maxLen, double Einit, int partId, double critVal){
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
		
		//showerScint.push_back(gen->Poisson((double) 16000)*layerTrackLen_scint(radLen2Long(i),radLen2Long(i+1))*nChargeTrack*0.12*0.15); //Get scintilation photons detected 
		//ionELoss(inShower,radLen2Long(i),radLen2Long(i + 1)); //Simulate Ionization energy loss
		//showerScint.push_back(gen->Poisson(16000*(39.6/25)*nChargeTrack*0.12*0.15); //Homogenous Cal
		showerAction1d(inShower,93.11,(double) X0, gen, false);
	}

	//return showerScint;
}*/


//Shower and record total number of scintilation photons produced and store information about the number of particles produced (Multithreaded)
void scintShower_Thread_Full(std::vector<int> &showerScint, std::vector<int> &partNum, int maxLen, double Einit, int partId, double critVal, TRandom3 *gen){
	//std::vector<int> showerScint; //Number of scintilation photons
	std::cout << "Starting Single Shower" << std::endl;
	particleR2 initPart = particleR2(Einit,0,partId); //Set up initial particle
	showerR2 inShower = showerR2(initPart); //Set up shower object
	int nChargeTrack = inShower.chargedTracks(); //Get the initial number of charged tracks

	//Shower through the calorimter
	//for (int i = 0; i < maxLen; i++){
	while (nChargeTrack > 0){
		showerAction1d(inShower,93.11,X0, gen, true); //Showering occurs
		nChargeTrack = inShower.chargedTracks(); //Count the number of charged tracks
		std::cout << "There are " << nChargeTrack << std::endl;
		//if (nChargeTrack == 0) break; //End shower if there are no more charged particles
		//std::cout << "Track Length between " << i << " and " << ++i << " = " << layerTrackLen_scint(radLen2Long(i),radLen2Long(i+1)) << std::endl;
		
		//showerScint.push_back(gen->Poisson((double) 16000)*layerTrackLen_scint(radLen2Long(i),radLen2Long(i+1))*nChargeTrack*0.12*0.15); //Get scintilation photons detected 
		//ionELoss(inShower,radLen2Long(i),radLen2Long(i + 1)); //Simulate Ionization energy loss
		//showerScint.push_back(gen->Poisson((double) 16000*(39.6/25)*nChargeTrack)*0.12*0.15); //Homogenous Cal with Poission Fluctuations
		for (int i = 0; i < inShower.showerSize();i++){ //Loop through the shower
			if (abs(inShower.idVec.at(i)) == 11){
				showerScint.push_back(gen->Poisson(16000*trackLen_scint(inShower.locVec.at(i) - inShower.diffLocVec.at(i),inShower.locVec.at(i)))*nChargeTrack*0.12*0.15); //Sampling Cal with Poission Fluctuations
			}else continue;
		}
		partNum.push_back(inShower.showerSize());
	}
}


int main(){
	//std::cout << "Test" << std::endl;
	//Test a very basic mean showering model of a 50 GeV positron in 1d going through 24 radiation lengths
	double E0 = 500; //Starting energy in MeV
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
	//std::vector<double> inputE = linspace(500e-3,1,2);
	std::vector<double> sumScintPhoto;

	//long testPhoton = scintShower(25,E0,11);
	//std::cout << testPhoton << std::endl;

	//Full code Starts here
	//ROOT Objects
	TFile *f1 = new TFile("ScintPhotoOut_Sample_Final_NoPoisson.root","RECREATE");
	//TCanvas *c1 = new TCanvas("c1","c1",500,500);
	TRandom3 *randGen = new TRandom3();

	//Add vector of energies to the ROOT file
	f1->WriteObject(&inputE,"Initial_E");
	int showerNum = 0;

	//Loop over all incident particle energies
	for ( double E : inputE) {
		TString crntNameScintPhoto; crntNameScintPhoto.Form("Shower_%d",showerNum);
		TString crntNameTotalPart; crntNameTotalPart.Form("PartNum_%d",showerNum);
		std::cout << "E = " << E << " GeV" << std::endl;

		//Get stopping length from the analytic model
		double stopX0 = log(E*1e3/93.11)/log(2);
		double stopLen = radLen2Long(stopX0);

		//Create 2 threads for each input particle
		//std::thread t1,t2;
		double photoSum = 0;
		std::vector<int> showerVecScint_elec, showerVecScint_pos,totalVecScint; //Scintilation Photon Vectors
		std::vector<int> showerVecPart_elec, showerVecPart_pos,totalVecPart; //Total number of particles
		std::vector<int> showerVecScintArr[2]; //
		std::vector<int> showerVecPartArr[2]; //

		std::thread t1(scintShower_Thread_Full,std::ref(showerVecScint_elec),std::ref(showerVecPart_elec),25,E*1e3,11,93.11,randGen);
		std::thread t2(scintShower_Thread_Full,std::ref(showerVecScint_pos),std::ref(showerVecPart_pos),25,E*1e3,11,93.11,randGen);

		t1.join();
		t2.join();

		//Single Threading
		//scintShower_Thread_Full(showerVecScint_elec,showerVecPart_elec,25,E*1e3,11,93.11,randGen);
		//scintShower_Thread_Full(showerVecScint_pos,showerVecPart_pos,25,E*1e3,11,93.11,randGen);

		//Update the shower vector array with the first index benig the largest vector
		if (showerVecScint_elec.size() > showerVecScint_pos.size()){
			showerVecScintArr[0] = showerVecScint_elec;
			showerVecScintArr[1] = showerVecScint_pos;
			showerVecPartArr[0] = showerVecPart_elec;
			showerVecPartArr[1] = showerVecPart_pos;
		}else if (showerVecScint_elec.size() < showerVecScint_pos.size()){
			showerVecScintArr[0] = showerVecScint_pos;
			showerVecScintArr[1] = showerVecScint_elec;
			showerVecPartArr[0] = showerVecPart_pos;
			showerVecPartArr[1] = showerVecPart_elec;
		}else if (showerVecScint_elec.size() == showerVecScint_pos.size()){
			showerVecScintArr[0] = showerVecScint_elec;
			showerVecScintArr[1] = showerVecScint_pos;
			showerVecPartArr[0] = showerVecPart_elec;
			showerVecPartArr[1] = showerVecPart_pos;
		}

		//Make scintilation photon adjustments to each shower
		/*double diffArr[2]; //diff0, diff1;
		double scintDiffArr[2]; //scintDiff0, scintDiff1;
		for (int i = 0; i < 2; i++) {
			diffArr[i] = stopX0 - (((double) showerVecScintArr[i].size()) - 1);
			std::cout << diffArr[i] << std::endl;
			//scintDiffArr[i] = showerVecScintArr[i].at(showerVecScintArr[i].size() - 1) + diffArr[i];
		}*/

		//Get the number of photons for each shower
		/*for (int i = 0; i < 2; i++){
			for (int j = 0; j < showerVecScintArr[i].size(); j++){
				showerVecScintArr.at(j)
			}
		}*/

		
		//if (showerVecScint_pos.size() != showerVecScint_elec.size()) std::cout <<"size difference" << std::endl;

		//scintShower(25,E*1e3,11); //Old single threaded shower
		//Get total number of scintilation photons for a LLP event
		totalVecScint = showerVecScintArr[0];
		for (int i = 0; i < showerVecScintArr[1].size(); i++){totalVecScint.at(i) = totalVecScint.at(i) + showerVecScintArr[1].at(i);}

		//Get total number of particles for a LLP event
		totalVecPart = showerVecPartArr[0];
		for (int i = 0; i < showerVecPartArr[1].size(); i++){totalVecPart.at(i) = totalVecPart.at(i) + showerVecPartArr[1].at(i);}
		//for (int i = 0; i < showerVecScint_elec.size(); i++){totalVecScint.push_back(showerVecScint_elec.at(i) + showerVecScint_pos.at(i));}

		//Write current vectors
		f1->WriteObject(&totalVecScint,crntNameScintPhoto); 
		f1->WriteObject(&totalVecPart,crntNameTotalPart); 
		
		showerNum++;
		//for (int phot : totalVec){std::cout << phot << std::endl;}
	} //Fill Scintilaiton photon vector
	
	//Store Active and Passive Laye information
	std::vector<double> activeVec;
	std::vector<double> passiveVec;
	for (int i = 0; i < 25; i++){activeVec.push_back(trackLen_scint(radLen2Long((double) i), radLen2Long((double) i + 1)));passiveVec.push_back(trackLen_pb(radLen2Long((double) i), radLen2Long((double) i + 1)));}

	f1->WriteObject(&activeVec,"ActiveGeo");
	f1->WriteObject(&passiveVec,"PassiveGeo");

	delete randGen;
	f1->Close();
	//delete f1;

	//Informal Debugging
	
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
		showerAction1d(inShower, , (double) t, false);
		if (nChargeTrack == 0) {std::cout << "Ending Shower" << std::endl; break;} //End shower

		//Check energy conservation 
		//if (Ecrt != E0)std::cout << "!Energy is Not being conserved!\n" << Ecrt << " != " << E0 << std::endl;
		

		double nGamma = 0; //Number of Scintilation photons
		//nGamma = 16000*layerTrackLen_scint(radLen2Long(t),radLen2Long(t+1))*nChargeTrack;
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