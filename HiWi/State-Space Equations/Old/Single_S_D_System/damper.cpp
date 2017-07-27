/*
 * damper.cpp
 *
 *  Created on: Dec 18, 2012
 *      Author: aditya
 */

#include <stdio.h>
#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include "init.h"

void getGroundProps(double time, double stepHeight, double* displacement, double *velocity);

int main(int argc, char *argv[]) {

	/**
	 * Values are to be read from the input file.
	 */
	double stiffCoeff;		// Stiffness Coefficient.
	double dampCoeff;		// Damping Coefficient.
	double mass; 			// Mass of the object.
	double initialTension; // Initial tension in the spring used to calculate the free length of the spring
	double springFreeLen;	// Free length of the spring :: Should be more than 500 -> mass position

	double tStart;			// Simulation Start time
	double tEnd;			// Simulation End time
	double stepHeight;		// Height of the step input
							// Always maximum at 0.05 seconds after tStart
	double dt;				// time Step
	double g;				// Gravity constant in -y direction

	double massIniPosition = 0.5; // Initial position of the mass
	double springLength = 0;	   // Length of the spring used to calculate spring force

	double simTime = tStart; 	// Current simulation time
	long int count = 0; 		// Iteration count

	std::vector<double> position; 	// Vector holding position of the mass with time
	std::vector<double> acc;  	// Vector holding Acceleration of the mass with time
	std::vector<double> vel;  	// Vector holding Velocity of the mass with time
	std::vector<double> Stime; 	// Vector holding time.
	std::vector<double> force; 	// Vector holding force on the mass.


	/**
	 * Checking if we have both input and output file name.
	 */
	if(argc < 2){
		std::cout<<"Please provide 1. Input file name with path"<<std::endl;
		std::cout<<"               2. Output file name(limited to 90 Characters)"<<std::endl;
		return 1;
	}

	char *inputFile = argv[1]; // Name of the input file including path.
	char *outputFile = argv[2]; // Name of the output file to be written.

	/**
	 *  Initializing the velocity and displacement
	 *  to Zero at the start of Simulation
	 */
	position.push_back(massIniPosition);
	vel.push_back(0.0);
	Stime.push_back(tStart);

	/**
	 * Reading the input file.
	 */
	std::cout<<std::endl;
	read_parameters(inputFile, &stiffCoeff, &dampCoeff, &mass, &tEnd, &tStart, &stepHeight, &dt, &g, &initialTension);
	std::cout<<std::endl;
	double dt_original = dt;
	/**
	 * Calculating the initial length of the spring
	 */
	springFreeLen = massIniPosition + (initialTension/stiffCoeff);

	// Calculating the gravitation force on the mass
	double gravityForce = mass*g;
	/**
	 * Time Loop
	 */
	while(simTime < tEnd){

		/**
		 * Temporary variables holding velocity, displacement and accelerations
		 */
		double xGround;		// Input Displacement
		double vGround;		// Input Velocity
		double x_n1;	// Displacement in next time step
		double v_n1;	// Velocity in next time step
		double acc_n;	// Acceleration in current time step
		/**
		 * Setting the input displacement and velocity
		 * according to the step input
		 */
		getGroundProps(simTime, stepHeight, &xGround, &vGround);

		// Simple Euler time step to calculate velocity in next time step
		acc_n = -(stiffCoeff*(position.at(count) - xGround - springFreeLen)
				+ dampCoeff*(vel.at(count) - vGround)
				+ gravityForce )/mass;
		// Updating the values in the vectors
		acc.push_back(acc_n);

		v_n1 = vel.at(count) + dt*acc.at(count);
		// Updating velocity
		vel.push_back(v_n1);

		// Explicit Euler time step for displacement
		x_n1 = position.at(count) + dt*vel.at(count);
		// Updating Displacement
		position.push_back(x_n1);

		// Updating acceleration
		force.push_back(-acc_n*mass);

		simTime +=  dt;
		Stime.push_back(simTime);
		// Updating the iteration count
		count ++;

	} // End of time Loop


	// Writing the output in the specified file.
	write_output(outputFile,Stime, position, vel, acc, force, massIniPosition);

	std::cout<<"Successfully Completed Simulation !!!"<<std::endl;
	std::cout<<std::endl;

	return 1;
}


/**
 * Method for obtaining the ground displacement and velocity depending on time
 */
void getGroundProps(double time, double stepHeight, double* displacement, double *velocity){

	double stepLength = 0.05;

	if(time < stepLength){
		*displacement = -1.5*(stepHeight/0.05/1000)*time + stepHeight/1000;
		*velocity = -1.5*(stepHeight/0.05/1000);
	} else {
		*displacement  = -stepHeight/1000;
		*velocity = 0;
	}

}

