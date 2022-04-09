#ifndef CALGEO_H
#define CALGEO_H

//Determine start and end point of a given layer
double *layerTerminus(double startPoint);

//Obtain the Scintilating track length for a given layer
double layerTrackLen(double startPoint);

//For a showering start and end point determine the track length of material
double trackLen(double startPoint, double endPoint);

#endif
