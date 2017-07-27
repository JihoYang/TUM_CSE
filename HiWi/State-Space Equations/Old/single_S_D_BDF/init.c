#include "helper.h"
#include "init.h"
#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>

int read_parameters( const char *inFileName,       /* name of the file */
                    double *stiffCoeff,
                    double *dampCoeff,
                    double *mass,
                    double *tEnd,
                    double *tStart,
                    double *stepHeight,
                    double *dt,
                    double *g,
                    double *initialTension
                    ) {

   READ_DOUBLE( inFileName, *stiffCoeff );
   READ_DOUBLE( inFileName, *dampCoeff );
   READ_DOUBLE( inFileName, *mass );
   READ_DOUBLE( inFileName, *tEnd );
   READ_DOUBLE( inFileName, *tStart );
   READ_DOUBLE( inFileName, *stepHeight );
   READ_DOUBLE( inFileName, *dt);
   READ_DOUBLE( inFileName, *g);
   READ_DOUBLE( inFileName, *initialTension);

   return 1;
}


int write_output( char* outFileName, std::vector<double> Stime, std::vector<double> position, std::vector<double> vel,
				   std::vector<double> acc, std::vector<double> force, double massIniPosition) {

	char oFile[100] = "";
	sprintf(oFile,"%s.dat",outFileName);
	// Opening the file
	FILE *fp = fopen (oFile, "w");

	// Writing a small header
	fprintf(fp,"***********************************************************\n");
	fprintf(fp,"******Simulation Results Single Damper System**************\n\n");
	// Writing the column headings
	fprintf(fp,"Time\tDisp\tVel\tAcc\tForce\n");

	// Writing the values in
	for(int i=1; i<Stime.size()-1; i++){
		fprintf(fp,"%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n",Stime.at(i),((massIniPosition - position.at(i))*1000),vel.at(i),acc.at(i),force.at(i));
	}


	// Closing the file
	fclose (fp);

	return 1;
}




