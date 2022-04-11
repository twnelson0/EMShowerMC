#ifndef CALGEO_H
#define CALGEO_H

//All units are implied to be longitudnal distance expcet for radLen2Long or otherwise specified

//Convert Radiation length to Longitudnal Distance
double radLen2Long(double t);

//Convert longitudnal distance to Radiation Lengths
double long2RadLen(double z);

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

//Determine which scintilating layer a small length is contained within (input in rad len)
//int scintNumber(double startVal, double endVal);

#endif
