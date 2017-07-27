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

	double massIniPosition = 0.500; // Initial position of the mass
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
	 * Reading the input file.
	 */
	std::cout<<std::endl;
	read_parameters(inputFile, &stiffCoeff, &dampCoeff, &mass, &tEnd, &tStart, &stepHeight, &dt, &g, &initialTension);
	std::cout<<std::endl;

	/**
	*  Initializing the velocity and displacement
	*  to Zero at the start of Simulation
	*/
	Stime.push_back(tStart);
	position.push_back(massIniPosition);
	vel.push_back(0.0);
	acc.push_back(0.0);
	/**
	 * Calculating the initial length of the spring
	 */
	springFreeLen = massIniPosition + (initialTension/stiffCoeff);
	// Calculating the gravitation force on the mass
	double gravityForce = mass*g;
	simTime = tStart;
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

			// %%%%%%%%%%%%%  FIRST SUB STEP %%%%%%%%%%%%%%
			double k1 = -(stiffCoeff*(position.at(count) - xGround - springFreeLen) + dampCoeff*(vel.at(count) - vGround) + gravityForce)/mass;
			double v1 = vel.at(count) + 0.5*dt*k1;
			double x1 = position.at(count) + 0.5*dt*v1 + 0.5*0.5*dt*dt*k1/2;

			// %%%%%%%%%%%%%  SECOND SUB STEP %%%%%%%%%%%%%%
			double k2 = -(stiffCoeff*(x1 - xGround - springFreeLen) + dampCoeff*(v1 - vGround) + gravityForce)/mass;
			double v2 = vel.at(count) + 0.5*dt*k2;
			double x2 = position.at(count) + 0.5*dt*v2 + 0.5*0.5*dt*dt*k2/2;

			// %%%%%%%%%%%%%  THIRD SUB STEP %%%%%%%%%%%%%%
			double k3 = -(stiffCoeff*(x2 - xGround - springFreeLen) + dampCoeff*(v2 - vGround) + gravityForce)/mass;
			double v3 = vel.at(count) + 1.0*dt*k3;
			double x3 = position.at(count) + 1.0*dt*v3 + + 0.5*0.5*dt*dt*k3/2;

			// %%%%%%%%%%%%%  FOURTH SUB STEP %%%%%%%%%%%%%%
			double k4 = -(stiffCoeff*(x3 - xGround - springFreeLen) + dampCoeff*(v3 - vGround) + gravityForce)/mass;


			// Assembling the RK-4th order method
			 v_n1 = vel.at(count) + dt*(k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;

			 acc_n = (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;
		// Updating the values in the vectors
		acc.push_back(acc_n);
		double diff_acc = (acc.at(count+1) - acc.at(count))/dt;

		// Updating velocity
		vel.push_back(v_n1);
		x_n1 = position.at(count) + dt*(vel.at(count) + 2.0*v1 +2.0*v2 + v_n1)/6;
//		x_n1 = position.at(count) + dt*(v1 + 2.0*v2 +2.0*v3 + v_n1)/6;
		// Explicit Euler time step for displacement
		//x_n1 = position.at(count) + dt*vel.at(count) + dt*dt*acc_n/4 + dt*dt*dt*diff_acc/6;
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
		*displacement = (stepHeight/0.05/1000)*time;
		*velocity = (stepHeight/0.05/1000);
	} else {
		*displacement = stepHeight/1000;
		*velocity = 0;
	}
}

