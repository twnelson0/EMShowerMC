#ifndef MATTERCALC_H
#define MATTERCALC_H

#include <math.h>
#include <vector>
//#include <iostream>

//Set up Constants
namespace constants{
	const double alpha = 1/137.;
	const double m_e = 0.511; //MeV
	const double N_Av = 6.022e23;
	const double r0 = 2.818e-15; //Classical Electron radius
	const double X0_pb = 6.37; //Radiation length of lead in g/cm^2
}

double radLen(uint Znuc, double Ar, double rho, bool densityBool);

double effRadLen(std::vector<double> X0Vec, std::vector<double> lenVec); //Effective radiation length for multi layered object (such as a colorimeter)

#endif
