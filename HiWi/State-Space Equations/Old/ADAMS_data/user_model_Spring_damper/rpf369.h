#ifndef _RPF369_H
#define _RPF369_H


const int UNLOADED_RADIUS	= 7;
const int WIDTH				= 8;
const int ASPECT_RATIO		= 10;
const int MASS				= 15;
const int SPRING_CONSTANT_VERTICAL	= 16;
const int DAMPING_COEFF     = 17;


/*Max. string length*/
int MaxStringLen = 128;
int MaxAttributeNameLen = 128;

/*To read [MODEL] block*/
int ModelBlockNameLen = 5;
int UsemodeAttributeNameLen = 12;

/*To read [DIMENSION] block*/
int DimensionBlockNameLen = 9;
int UnloadedradiusAttributeNameLen = 15;
int WidthAttributeNameLen = 5;
int AspectratioAttributeNameLen = 12;

/*To read [PARAMETER] block*/
int ParameterBlockNameLen = 9;
int MassNameLen = 4;
int SpringConstVertNameLen = 24;
int DampingCoefNameLen = 13;

/*To read Table block*/
int ShapeBlockNameLen = 5;
char Format[80];
int Formatlen;

#endif