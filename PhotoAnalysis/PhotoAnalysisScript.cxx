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

//Sum Photons for a given 


int main(){
	std::cout << "Test" << std::endl; 
	TFile *f = TFile::Open("ScintPhotoOut1.root","READ");
	std::cout << indxToEnergy(f, 0) << std::endl;
	std::map<double, int> EMap = EToIndx(f);

	for (const auto &p : EMap){
		std::cout << p.first << " --> " << p.second << std::endl;
	}

	/*std::vector<int> *v;
	f->GetObject("Shower_0",v);

	for (std::vector<int>::iterator it = v->begin(); it != v->end();it++){
		std::cout << *it << std::endl;
	}*/

	//Meomry
	delete f;
}

