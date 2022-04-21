#ifndef MATTERCALC_H
#define MATTERCALC_H

#include <math.h>
#include <vector>
//#include <map>
//#include <iostream>

//Set up Constants that will be useful
namespace physConstants{
	const double alpha = 1/137.;
	const double m_e = 0.511; //MeV
	const double N_Av = 6.022e23;
	const double r0 = 2.818e-15; //Classical Electron radius
	const double X0_pb = 6.37; //Radiation length of lead in g/cm^2
	const double rho_pb =  11.35; //Density in g/cm^3
	const double rho_scint = 1.08; //For Polyesterine based Scintilators
	const double Z_pb = 82;
	const double c = 299792458; //c in m/s (NIST value)
	//const double X0_Sint = ;//Radiation length of scintlator materail
	//const double eCrit_pb = ; //Critical energy of lead
}

double radLen(uint Znuc, double Ar, double rho, bool densityBool);

double effRadLen(std::vector<double> X0Vec, std::vector<double> lenVec); //Effective radiation length for multi layered object (such as a colorimeter)

std::vector<double> linspace(double startVal, double endVal, uint size = 50);

//double critE(); //Critical energy of electron/positron to bremstralung

#endif
