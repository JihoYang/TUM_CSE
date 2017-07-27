#ifndef _MAIN_H
#define	_MAIN_H

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
enum BOOLEAN{
	FALSE, /*FALSE = 0, TRUE = 1*/
	TRUE
};

#define ONE		1.e0
#define ZERO	0.0e0
#define MAXDIV	10

void act_line_parse(char* Value, int ValueLen);

void fialafx(double* fx, double* kappa, double* mu, 
			 double* fz, double* c_slip);

void fialafy(double* fy,	double* tz,		double* alpha, 
			 double* mu,	double* fz,		double* tirew,		
			 double* c_alpha);

void fialaty(double* ty, double* wy, double* fz, double* rolco);

double sign(double werte);

#endif