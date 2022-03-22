#include <iostream>
#include "MatterCalc.h"
#include "MatterCalc.cxx" //Probably poor form

int main(){
	std::cout << "The radiation length of Lead is " << radLen(82,207.2,11.350,false) << std::endl;
	//std::cout << "Alpha = " << constants::alpha << std::endl;

	return 0;
}
