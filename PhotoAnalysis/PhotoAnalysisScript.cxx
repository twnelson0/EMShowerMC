//C++ libaries
#include <iostream>
#include <vector>
#include <map>

//ROOT Libaries
#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TString.h"
#include "TGraph.h"


//Get intial e+/e- energy from given index value
double indxToEnergy(TFile *f, int indx){
	std::vector<double> *eVec;
	f->GetObject("Initial_E",eVec);

	return eVec->at(indx);
}

//Energy to index map
std::map<double, int> EToIndx(TFile *f1){
	//Energy to index map
	std::map<double, int> outMap;

	//Energy vector
	std::vector<double> *eVec;
	f1->GetObject("Initial_E",eVec);
	std::vector<double> tmp = *eVec;

	int indx = 0;
	for (double E : tmp){
		outMap[E] = indx;
		indx++;
	}

	return outMap;
}

//Energy Estmimation
/*double EnergyReco_1(int totalNGamma){
	double EHat;
	/*
	Use simple algorithm, number of particles = 2^N
	Get a fraction of particles that are charged leptons
	
}*/

//Sum Photons for a given energy value
int photoSum(TFile *f, double EVal){
	int sumVal = 0;

	//Get the Energy to index map
	std::map<double, int> eMap = EToIndx(f);

	//Obtain the histogram title
	TString vecName; vecName.Form("Shower_%d",eMap[EVal]);
	std::vector<int> *tempVec;
	f->GetObject(vecName,tempVec);

	//Sum over all Scintilation Photons
	for (std::vector<int>::iterator it = tempVec->begin(); it != tempVec->end(); it++){
		sumVal+=*(it);
	}

	return sumVal;
} 


//Generate and Save singal shower Graph
void showerPlot(TFile *f, double EVal, bool verbose = false){
	//ROOT Objects
	TCanvas *c1 = new TCanvas("c1","c1",1000,1000);

	//Get the Energy to index map
	std::map<double, int> eMap = EToIndx(f);

	//Obtain the TSTrings and vector name in file
	TString vecName, graphTitle,fileName; vecName.Form("Shower_%d",eMap[EVal]);
	graphTitle.Form("%.0lf GeV Event",EVal*2);
	fileName.Form("%.0lfGeVEventPhotoPlot.png",EVal*2);
	std::vector<int> *EVec;
	f->GetObject(vecName,EVec);
	if (verbose) std::cout << vecName << std::endl;
	
	//Vectors/Arrays for graph
	std::vector<double> radLenVec, photoVec;
	int radLenVal = 0;

	for (std::vector<int>::iterator it = EVec->begin(); it != EVec->end(); it++){
		radLenVec.push_back((double)radLenVal);
		photoVec.push_back((double) *it);
		//if (verbose) {std::cout << *it << std::endl;}
		radLenVal++;
	}
	
	//Graph
	TGraph *g1 = new TGraph(radLenVec.size(),&(radLenVec[0]),&(photoVec[0]));
	g1->SetTitle(graphTitle);
	g1->GetXaxis()->SetTitle("Radiation Lengths Traversed (X_{0})");
	g1->GetYaxis()->SetTitle("Number of Scintillation Photons");
	c1->SetLogy();
	g1->Draw("A*");
	c1->SaveAs(fileName);

	//Memory
	delete g1;
	c1->Close();

}


int main(){
	std::cout << "Test" << std::endl; 
	TFile *f = TFile::Open("ScintPhotoOut_HomogTest_6.root","READ");
	std::cout << indxToEnergy(f, 0) << std::endl;
	std::map<double, int> EMap = EToIndx(f);

	/*for (const auto &p : EMap){
		std::cout << p.first << " --> " << p.second << std::endl;
	}*/

	//showerPlot(f,500);

	std::vector<double> sumPhotVec;
	std::vector<double> energyVec;
	//Get Sum of Photons Vector
	for (const auto &mapPair : EMap){
		sumPhotVec.push_back((double) photoSum(f,mapPair.first));
		//std::cout << sumPhotVec.at(i) << std::endl;
		energyVec.push_back((mapPair.first)*2);
	}

	//Get all the shower plot 
	for (const auto &mapPair : EMap){
		/*if (E == 1000 || E == 3000) showerPlot(f,E/2,true);
		else showerPlot(f,E/2);*/
		showerPlot(f,mapPair.first,true);
	}

	//Make the plot
	TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
	TGraph *g1 = new TGraph(EMap.size(),&(energyVec[0]),&(sumPhotVec[0]));
	g1->SetTitle("LLP Energy Vs. Total Number of Scintillationion Photons");
	g1->GetXaxis()->SetTitle("LLP Energy (GeV)");
	g1->GetYaxis()->SetTitle("Sum of Scintillation Photons");
	g1->Draw("A*");
	c1->SaveAs("ScintPhotoSumPlot_Homog4.png");

	delete g1;
	c1->Close();
	

	/*std::vector<int> *v;
	f->GetObject("Shower_0",v);

	for (std::vector<int>::iterator it = v->begin(); it != v->end();it++){
		std::cout << *it << std::endl;
	}*/

	//Meomry
	delete f;
}

