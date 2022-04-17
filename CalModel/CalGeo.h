#ifndef CALGEO_H
#define CALGEO_H

//Indicator Functions
int ind(double x, double lowBound, double upBound); //Returns 1 if x \in [a,b]

int indSoftUppr(double x, double lowBound, double upBound); // Returns if x \in [a,b)

//All units are implied to be longitudnal distance expcet for radLen2Long or otherwise specified
//All Longintudal distances are assumed to be in cm

//Convert Radiation length to Longitudnal Distance
double radLen2Long(double t);

//Convert longitudnal distance to Radiation Lengths
double long2RadLen(double z);

//Convert radiation length interval to time given a particle's energy
double radLen2Time(double tInter, double E, double mass);

//Obtain layer a given point is located within
int crntLayer(double crntPoint);

//Determine start and end point of a given layer
double *layerTerminus(double startPoint);

//Obtain the Scintilating track length for a given layer (will also work for end points)
double layerTrackLen_scint(double edgeVal, bool start = true);

//For a showering start and end point determine the Scintiliating track length of material
double trackLen_scint(double startPoint, double endPoint);

//Obtain the Lead track length for a given layer (will also work for end points)
double layerTrackLen_pb(double edgeVal, bool start = true);

//For a showering start and end point determine the Lead track length of material
double trackLen_pb(double startPoint, double endPoint);

//Print current layer material
void currentLayerMat(double crntPoint);

//Determine which scintilating layer a small length is contained within (input in rad len)
//int scintNumber(double startVal, double endVal);

#endif
