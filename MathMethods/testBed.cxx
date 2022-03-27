#include <iostream>
#include "MatterCalc.h"
#include "MatterCalc.cxx" //Probably poor form
#include "particles.h"

int main(){
	//std::cout << "The radiation length of Lead is " << radLen(82,207.2,11.350,false) << std::endl;
	//std::cout << "Alpha = " << constants::alpha << std::endl;

	elecR2 elec = elecR2(1,0,false);
	elecR2 pos = elecR2(1,0,true);
	std::cout << elec.idNum() << std::endl;
	std::cout << pos.idNum() << std::endl;

	return 0;
}
