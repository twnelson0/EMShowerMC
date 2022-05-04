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
#include "TF1.h"
#include "TGraphErrors.h"


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
uint photoSum(TFile *f, double EVal){
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

//Mean of a vector
double vecMean(std::vector<auto> *inVec){
	double meanVal = 0;
	for (int i = 0; i < inVec->size(); i++){meanVal += inVec->at(i);}
	return meanVal/inVec->size();
}

//Mean of square of vecor of data
double vecSqrMean(std::vector<auto> *inVec){
	double meanVal = 0;
	for (int i = 0; i < inVec->size(); i++){meanVal += pow(inVec->at(i),2);}
	return meanVal/inVec->size();
}

//Sample Standard Deviation
double vecStdDev(std::vector<auto> *inVec){
	return vecSqrMean(inVec) - pow(vecMean(inVec),2);
}


//Generate and Save singal shower Graph
void showerScintPlot(TFile *f, double EVal, bool verbose = false){
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
		if (verbose) {std::cout << *it << std::endl;}
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

//Generate and Save singal shower Graph
void showerPartPlot(TFile *f, double EVal, bool verbose = false){
	//ROOT Objects
	TCanvas *c1 = new TCanvas("c1","c1",1000,1000);

	//Get the Energy to index map
	std::map<double, int> eMap = EToIndx(f);

	//Obtain the TSTrings and vector name in file
	TString vecName, graphTitle,fileName; vecName.Form("PartNum_%d",eMap[EVal]);
	graphTitle.Form("%.0lf GeV Event",EVal*2);
	fileName.Form("%.0lfGeVEventPartPlot.png",EVal*2);
	std::vector<int> *EVec;
	f->GetObject(vecName,EVec);
	if (verbose) std::cout << vecName << std::endl;
	
	//Vectors/Arrays for graph
	std::vector<double> radLenVec, partVec;
	int radLenVal = 0;

	for (std::vector<int>::iterator it = EVec->begin(); it != EVec->end(); it++){
		radLenVec.push_back((double)radLenVal);
		partVec.push_back((double) *it);
		if (verbose) {std::cout << *it << std::endl;}
		radLenVal++;
	}
	
	//Graph
	TGraph *g1 = new TGraph(radLenVec.size(),&(radLenVec[0]),&(partVec[0]));
	g1->SetMarkerColor(4);
	g1->SetMarkerStyle(21);

	//Draw Graph
	g1->SetTitle(graphTitle);
	g1->GetXaxis()->SetTitle("Radiation Lengths Traversed (X_{0})");
	g1->GetYaxis()->SetTitle("Number of Particles in Shower");
	c1->SetLogy();
	g1->Draw("A*");

	//Set up analytic function
	TF1 *expf = new TF1("expf","2^(x + 1)",0,radLenVec.at(radLenVec.size() - 1));
	expf->Draw("same");
	c1->SaveAs(fileName);

	//Memory
	delete expf;
	delete g1;
	c1->Close();
}

//Generate and Save singal shower Histogram
void showerPhotoHist(TFile *f, double EVal, bool verbose = false){
	//ROOT Objects
	TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
	TH1D *h1 = new TH1D();

	//Get the Energy to index map
	std::map<double, int> eMap = EToIndx(f);

	//Obtain the TSTrings and vector name in file
	TString histName, histTitle,fileName; histName.Form("PhotoHist_%d",eMap[EVal]);
	histTitle.Form("%.0lf GeV Event",EVal*2);
	fileName.Form("%.0lfGeVEventPhotoHist.png",EVal*2);
	h1 = (TH1D*)f->Get(histName);
	if (verbose) std::cout << histName << std::endl;
	
	//Vectors/Arrays for graph
	std::vector<double> radLenVec, partVec;
	int radLenVal = 0;

	//Draw Graph
	h1->SetStats(0);
	h1->SetTitle(histTitle);
	h1->GetXaxis()->SetTitle("Radiation Lengths (X_{0})");
	//h1->GetYaxis()->SetTitle("Number of Scintillation Photons");
	//c1->SetLogy();
	h1->Draw();

	TF1 *expf = new TF1("expf","(0.15)^2*16000*2^(x + 1)*1.584",0,h1->GetNbinsX());
	//expf->Draw("same");

	/*for (std::vector<int>::iterator it = EVec->begin(); it != EVec->end(); it++){
		radLenVec.push_back((double)radLenVal);
		partVec.push_back((double) *it);
		if (verbose) {std::cout << *it << std::endl;}
		radLenVal++;
	}*/
	
	//Graph
	/*TGraph *g1 = new TGraph(radLenVec.size(),&(radLenVec[0]),&(partVec[0]));
	g1->SetMarkerColor(4);
	g1->SetMarkerStyle(21);

	*/

	//Set up analytic function
	//TF1 *expf = new TF1("expf","2^(x + 1)",0,radLenVec.at(radLenVec.size() - 1));
	//expf->Draw("same");
	c1->SaveAs(fileName);

	//Memory
	//delete expf;
	delete h1;
	delete expf;
	c1->Close();
}


//Generate and Save singal shower Histogram
void showerPartHist(TFile *f, double EVal, bool verbose = false){
	//ROOT Objects
	TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
	TH1D *h1 = new TH1D();

	//Get the Energy to index map
	std::map<double, int> eMap = EToIndx(f);

	//Obtain the TSTrings and vector name in file
	TString histName, graphTitle,fileName; histName.Form("PartHist_%d",eMap[EVal]);
	graphTitle.Form("%.0lf GeV Event",EVal*2);
	fileName.Form("%.0lfGeVEventPartPlot.png",EVal*2);
	std::vector<int> *EVec;
	//f->GetObject(vecName,EVec);
	h1 = (TH1D*)f->Get(histName);
	if (verbose) std::cout << histName << std::endl;
	
	h1->SetStats(0);
	h1->SetTitle(graphTitle);
	h1->GetXaxis()->SetTitle("Radiaton Lengths (X_{0})");
	//h1->GetYaxis()->SetTitle("Number of particles in shower");
	//c1->SetLogy();
	h1->Draw();

	//Set up analytic function
	//TF1 *expf = new TF1("expf","2^(x + 1)",0,h1->GetNbinsX());
	//expf->Draw("same");
	c1->SaveAs(fileName);

	//Memory
	//delete expf;
	delete h1;
	c1->Close();
}


void showerEnergyPlotHomog(TFile *f){
	TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
	std::map<double, int> EMap = EToIndx(f); //Energy map
	std::vector<double> sumPhotVec;
	std::vector<double> energyVec;
	
	//Get Sum of Photons Vector
	for (const auto &mapPair : EMap){
		sumPhotVec.push_back((double) photoSum(f,mapPair.first));
		energyVec.push_back((mapPair.first)*2);
	}

	//Set up the plot
	TGraph *g1 = new TGraph(energyVec.size(),&(energyVec[0]),&(sumPhotVec[0]));
	g1->SetTitle("Intial LLP Energy vs Total Number of Scintillationion Photons");
	g1->GetXaxis()->SetTitle("LLP Energy (GeV)");
	g1->GetYaxis()->SetTitle("Sum of Scintillation Photons");
	g1->Draw("A*");
	c1->SaveAs("ScintPhotoSumPlot_HomogFinal.png");

	//Memory
	delete g1;
	c1->Close();
}

double aHat(std::vector<double> x, std::vector<double> y, std::vector<double> sig){
	double mVal = 0;
	double denom = 0;
	for (int i = 0; i < x.size(); i++){
		mVal += x.at(i)*y.at(i)/pow(sig.at(i),2);
		denom += pow(x.at(i),2)/pow(sig.at(i),2);
	}

	return mVal/denom;
}

double var_a(std::vector<double> x, std::vector<double> sig){
	double varVal = 0;
	for (int i = 0; i < x.size(); i++){
		varVal += pow(x.at(i),2)/pow(sig.at(i),2);
	}

	return pow(varVal,-1);
}

std::vector<double> zeroFill(int N){
	std::vector<double> zeroVec;
	for (int i = 0; i < N; i++) {zeroVec.push_back(0);}

	return zeroVec;
}

double sigE(double N, double a, double aSig, double NSig){
	return 1/a*sqrt(pow(NSig,2) + pow(N,2)/pow(a,2)*pow(aSig,2));
}

void EReco(TFile *f){
	//Read in each energy value
	std::map<double, int> EMap = EToIndx(f);
	std::vector<double> errVec, meanNPhotoVec;
	std::vector<double> inEVec, resVec_Mine;
	double aVal, aSig;
	//std::vector<double> EOverNVec;

	for (const auto &mapPair : EMap){
		inEVec.push_back((mapPair.first)*2);
		//resVec_Tav.push_back(sqrt(93.11e-3/inEVec.at(i)));
	}

	std::vector<double> E_err = zeroFill(EMap.size());

	for (int i = 0; i < EMap.size(); i++){
		//Read out vectors
		TString vecName; vecName.Form("showerScintVec_%d",i);
		std::vector<int> *crntPhotoVec;
		f->GetObject(vecName,crntPhotoVec);

		//Get mean and sigma of data
		meanNPhotoVec.push_back(vecMean(crntPhotoVec));
		errVec.push_back(sqrt(vecStdDev(crntPhotoVec)));
	}

	//Make a plot of the LLP eneryg and number of scintilation photons
	TCanvas *c1 = new TCanvas("c1","c1",1500,1500);
	TGraphErrors *g1 = new TGraphErrors(inEVec.size(),&(inEVec[0]),&(meanNPhotoVec[0]),&(E_err[0]),&(errVec[0]));

	//Graph formatting
	g1->SetTitle("");
	g1->GetXaxis()->SetTitle("LLP Energy (GeV)");
	g1->GetYaxis()->SetTitle("Number of Scintillation Photons");
	g1->SetMarkerStyle(20);
	g1->Draw("ap");
	c1->SaveAs("Final_ELLP_Plot.png");

	delete g1;
	c1->Close();

	//Get estimate for estimator coefficent
	aVal = aHat(inEVec,meanNPhotoVec,errVec);
	aSig = sqrt(var_a(inEVec,errVec));

	std::cout << "aHat = " << aVal << std::endl;
	std::cout << "aSig = " << aSig << std::endl;

	//Plot the sig/E
	for (int i = 0; i < inEVec.size();i++){
		resVec_Mine.push_back(sigE(meanNPhotoVec.at(i),aVal,aSig,errVec.at(i))/inEVec.at(i));
		//std::cout << sigE(meanNPhotoVec.at(i),aVal,aSig,errVec.at(i)) << std::endl;
		//std::cout << sigE(meanNPhotoVec.at(i),aVal,aSig,errVec.at(i)) - << std::endl;
		//std::cout << pow(meanNPhotoVec.at(i),2) << std::endl;
		std::cout << meanNPhotoVec.at(i)/aVal << std::endl;
	}
	
	TCanvas *c2 = new TCanvas("c2","c2",1500,1500);
	TGraph *g2 = new TGraph(inEVec.size(),&(inEVec[0]),&(resVec_Mine[0]));
	c2->SetRightMargin(0.09);
	c2->SetLeftMargin(0.15);
	c2->SetBottomMargin(0.15);
	g2->SetTitle("");
	g2->GetXaxis()->SetTitle("LLP Energy (GeV)");
	g2->GetYaxis()->SetTitle("#\sigma/E");
	g2->SetMarkerStyle(21);
	g2->Draw("ap");

	//TF1 *tavRes = new TF1("tavRes","sqrt(93.11*10^(-3)/x)",1e3,10e3);
	//tavRes->Draw("SAME");
	c2->SaveAs("ResPlot.png");

	
	delete g2;
	//delete tavRes;
	c2->Close();
}

/*void activeGeo(TFile *f){
	std::vector<double> *activeLayers;
	double sumVal = 0;
	f->GetObject("ActiveGeo",activeLayers);

	for (std::vector<double>::iterator it = activeLayers->begin(); it != activeLayers->end(); it++){
		std::cout << *(it) << std::endl;
		sumVal += *(it);
	}

	std::cout << "Sum = " << sumVal << std::endl;
}*/


int main(){
	std::cout << "Test" << std::endl; 
	TFile *f = TFile::Open("ScintPhotoOut_Homogenous_Final_Poisson.root","READ");
	//TFile *f = TFile::Open("ScintPhotoOut_Sample_Final_Poisson.root","READ");
	std::cout << indxToEnergy(f, 0) << std::endl;
	std::map<double, int> EMap = EToIndx(f);

	/*for (const auto &p : EMap){
		std::cout << p.first << " --> " << p.second << std::endl;
	}*/

	//showerScintPlot(f,500);

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
		/*if (E == 1000 || E == 3000) showerScintPlot(f,E/2,true);
		else showerScintPlot(f,E/2);*/
		showerPhotoHist(f,mapPair.first,true);
		showerPartHist(f,mapPair.first);
	}

	//activeGeo(f);

	//Make the plot
	TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
	TGraph *g1 = new TGraph(EMap.size(),&(energyVec[0]),&(sumPhotVec[0]));
	g1->SetTitle("LLP Energy Vs. Total Number of Scintillationion Photons");
	g1->GetXaxis()->SetTitle("LLP Energy (GeV)");
	g1->GetYaxis()->SetTitle("Sum of Scintillation Photons");
	g1->Draw("A*");
	c1->SaveAs("ScintPhotoSumPlot_HomogFinal.png");

	delete g1;
	c1->Close();


	TFile *f2 = TFile::Open("MCOut.root");
	//EReco(f2);
	

	/*std::vector<int> *v;
	f->GetObject("Shower_0",v);

	for (std::vector<int>::iterator it = v->begin(); it != v->end();it++){
		std::cout << *it << std::endl;
	}*/

	//Memory
	delete f;
	delete f2;
}

