//Standard c++ libaries
#include <iostream>
#include <vector>

//Root and other Libaries
#include "TRandom3.h"
#include "../MathMethods/MatterCalc.h"
#include "../MathMethods/Particles.h"


/*
Outline
	* In this simplified model I will consider incident e-e+ particles on the calorimeter from the decay volume
	* Since these are the product of the decay of LLPs where in the decay volume they originate will be a exponentially random value with average of LLP lifetime
	* The angle of incidence of the LLP (and therefore the e+e-) will be a uniformly generated value within FASERs 2mrad acceptance
	* After each radiation length a given particle will lose half its energy
	* After 1 rad len a charged lepton will yield a photon and a charged lepton (Bremstralung)
	* After 1 rad len a photon will yield an electron positron pair (Pair production with nuclei really)
	* In principle pair production of photons should only happen every 9/7 Rad Len
	* In principle process for leptons stops at E_c and process stops for photons when E < 2m_e
	* When particles split I'm not sure how to handle direction, in principle if I consider a psudo 1d model all that matters is that energy splits in 2 mometa and location in detector irelevent
	* Negate more detailed interactions
	* In more detailed simulation when interaction occurs and resultent energies and directions should be more accuretly obtained from process cross section

Magnets will pull them apart

*/

//Particle Interaction
std::vector<?ParticelObj?> layerEdge(ParticleObj inPart){
	//Use angleular average and RMS info in Tav and SLAC paper 
	if (abs(inPart.idNum) == 11){ //Bremstralung

	}
	else if (inPart.idNum == 22){ //Pair production

	}
}

/*26/03/2022 22:39, I'm wondering if I need these particle objects in the E/2 splitting model they don't matter, it may matter for a more nuanced model but I'm wondering if there isn't
A more function appraoch to this problem

*/

/*
27/03/2022 16:04, I think I should start with a 1d simple mean showering program everything splits in half as it passes through the detector just to test some of the basics then I should add more nuanced behaviors
there is a lot of complexity present in showering and I don't want to waste my time figuring out the most ideal method of dealing with these sorts of interactions
*/


int main(){
	std::cout << "Test" << std::endl;
}
