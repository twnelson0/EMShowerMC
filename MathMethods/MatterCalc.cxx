#include <math.h>
#include "MatterCalc.h"

using namespace constants;

//Radiation length
double radLen(uint Znuc, double Ar, double rho, bool densityBool = false){
	double Z = (double) Znuc; //Convert Z to a usable format
	double X0Inv = (densityBool) ? 4*alpha*pow(r0*1e2,2)*((rho*N_Av)/Ar)*Z*(Z + 1)*log(183./(pow(Z,1/3.))) : 4*alpha*pow(r0*1e2,2)*(N_Av/Ar)*Z*(Z + 1)*log(183./(pow(Z,1/3.)));

	return pow(X0Inv,-1);
}
