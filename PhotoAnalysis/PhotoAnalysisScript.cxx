//C++ libaries
#include <iostream>
#include <vector>

//ROOT Libaries
#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"


int main(){
	std::cout << "Test" << std::endl; 
	TFile *f = TFile::Open("ScintPhotoOut_Test1.root","READ");
	std::vector<int> *v;
	f->GetObject("Shower_0",v);

	for (std::vector<int>::iterator it = v->begin(); it != v->end();it++){
		std::cout << *it << std::endl;
	}

	//Meomry
	delete f;
}

