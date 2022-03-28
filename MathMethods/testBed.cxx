#include <iostream>
#include "MatterCalc.h"
#include "MatterCalc.cxx" //Probably poor form
#include "particles.h"

int main(){
	//std::cout << "The radiation length of Lead is " << radLen(82,207.2,11.350,false) << std::endl;
	//std::cout << "Alpha = " << constants::alpha << std::endl;

	particleR2 testPhoton = particleR2(10,0,22);
	std::cout << testPhoton.id() << std::endl;
	std::cout << testPhoton.E() << std::endl;

	return 0;
}
