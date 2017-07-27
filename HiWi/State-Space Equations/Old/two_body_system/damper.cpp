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

#define numBodies 2

void getGroundProps(double time, double stepHeight, double* displacement, double *velocity);

int main(int argc, char *argv[]) {

	double tStart;			// Simulation Start time
	double tEnd;			// Simulation End time
	double stepHeight;		// Height of the step input
							        // Always maximum at 0.05 seconds after tStart
	double dt;				// time Step
	double g;				// Gravity constant in -y direction
	double simTime ; 		// Current simulation time
	long int count = 0; 	// Iteration count

	/**
	 * For Body 1
	 */
	double mass1;			// mass of the body
	double stiffCoeff1;		// Stiffness Coefficient.
	double initialTension1; // Initial tension in the spring used to calculate the free length of the spring
	double springFreeLen1;	// Free length of the spring :: Should be more than 500 -> mass position
	double dampCoeff1;		// Damping Coefficient.
	double massIniPosition1 = 0.500; // Initial position of the mass
	/**
	* For Body 2
	*/
	double mass2;			// mass of the body
	double stiffCoeff2;		// Stiffness Coefficient.
	double initialTension2; // Initial tension in the spring used to calculate the free length of the spring
	double springFreeLen2;	// Free length of the spring :: Should be more than 500 -> mass position
	double dampCoeff2;		// Damping Coefficient.
	double massIniPosition2 = 0.80; // Initial position of the mass

	std::vector<double> Stime; 	// Vector holding time.

	std::vector< std::vector<double> > vel(2); 		// Vector holding Velocity of the mass with time
	std::vector< std::vector<double> > position(2); 	// Vector holding position of the mass with time
	std::vector< std::vector<double> > acc(2);			// Vector holding Acceleration of the mass with time
	std::vector< std::vector<double> > force(2);		// Vector holding force on the mass.
	std::vector<double> stiffCoeff;
	std::vector<double> dampCoeff;
	std::vector<double> initialTension;
	std::vector<double> mass;
	std::vector<double> massIniPosition;
	std::vector<double> springFreeLen;
	std::vector<double> gravityFroce;
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
	read_parameters_1(inputFile, &stiffCoeff1, &dampCoeff1, &mass1, &tEnd, &tStart, &stepHeight, &dt, &g, &initialTension1);
	read_parameters_2(inputFile, &stiffCoeff2, &dampCoeff2, &mass2, &tEnd, &tStart, &stepHeight, &dt, &g, &initialTension2);
	std::cout<<std::endl;

	/**
	*  Initializing the velocity and displacement
	*  to Zero at the start of Simulation
	*/

	massIniPosition.push_back(massIniPosition1);
	massIniPosition.push_back(massIniPosition2);

	mass.push_back(mass1); 	stiffCoeff.push_back(stiffCoeff1);	dampCoeff.push_back(dampCoeff1); 	initialTension.push_back(initialTension1);
	mass.push_back(mass2);	stiffCoeff.push_back(stiffCoeff2);	dampCoeff.push_back(dampCoeff2);	initialTension.push_back(initialTension2);

	vel.at(0).push_back(0.0);	acc.at(0).push_back(0.0);
	vel.at(1).push_back(0.0);	acc.at(1).push_back(0.0);

	position.at(0).push_back(massIniPosition.at(0));
	position.at(1).push_back(massIniPosition.at(1));

	Stime.push_back(tStart);

	/**
	 * Calculating the initial length of the spring
	 */
	for(int i=0; i<numBodies; i++){
		springFreeLen.push_back( massIniPosition.at(i) + (initialTension.at(i)/stiffCoeff.at(i)) );
		gravityFroce.push_back( mass.at(i)*g );
	}

	/**
	 * Acceleration Evaluation Matrix
	 */
	double AEM[7][2] = { { stiffCoeff2/mass1, -stiffCoeff2/mass2 },				// x2
						  { dampCoeff2/mass1, -dampCoeff2/mass1 },				// v2
						  { -(stiffCoeff2+stiffCoeff1)/mass1, stiffCoeff2/mass2 },	//x1
						  { -(dampCoeff1+dampCoeff2)/mass1, dampCoeff2/mass2 },		//v1
						  //{ -(mass2*g - stiffCoeff2*springFreeLen2)/mass2, -(mass1*g + stiffCoeff2*springFreeLen2 - stiffCoeff1*springFreeLen1)/mass1 }, // constants
						  { (-mass1*g - stiffCoeff2*springFreeLen.at(1) + stiffCoeff1*springFreeLen.at(0))/mass1,
								  (-mass2*g + stiffCoeff2*springFreeLen.at(1))/mass2 },
						  { stiffCoeff1/mass1, 0},	// xg
						  { dampCoeff1/mass1, 0} //vg
						};


	/**
	 * Time Loop
	 */
	simTime = tStart;
	while(simTime < tEnd){

		/**
		 * Temporary variables holding velocity, displacement and accelerations
		 */
		double xGround;		// Input Displacement
		double vGround;		// Input Velocity
		std::vector<double> k1(2, 0.0);		std::vector<double> v1(2, 0.0); 	std::vector<double> x1(2, 0.0);
		std::vector<double> k2(2, 0.0);		std::vector<double> v2(2, 0.0);		std::vector<double> x2(2, 0.0);
		std::vector<double> k3(2, 0.0);		std::vector<double> v3(2, 0.0);		std::vector<double> x3(2, 0.0);
		std::vector<double> k4(2, 0.0);
		std::vector<double> v_n1, acc_n, x_n1;
		/**
		 * Setting the input displacement and velocity
		 * according to the step input
		 */
		getGroundProps(simTime, stepHeight, &xGround, &vGround);
		double DOF1[7] = {position.at(1).at(count), vel.at(1).at(count), position.at(0).at(count), vel.at(0).at(count), 1, xGround, vGround };

			// %%%%%%%%%%%%%  FIRST SUB STEP %%%%%%%%%%%%%%
				for(int j=0; j<2; j++){
					for(int i=0; i<7; i++){
						k1.at(j) = k1.at(j) + AEM[i][j]*DOF1[i] ;

					}
					v1.at(j) = vel.at(j).at(count) + 0.5*dt*k1.at(j);
					x1.at(j) = position.at(j).at(count) + 0.5*dt*v1.at(j) + 0.5*0.5*dt*dt*k1.at(j)/2;
				}

		double DOF2[7] = { x1.at(1), v1.at(1), x1.at(0), v1.at(0), 1, xGround, vGround };

			// %%%%%%%%%%%%%  SECOND SUB STEP %%%%%%%%%%%%%%
				for(int j=0; j<2; j++){
					for(int i=0; i<7; i++){
						k2.at(j) = k2.at(j) + AEM[i][j]*DOF2[i] ;
					}

					v2.at(j) = vel.at(j).at(count) + 0.5*dt*k2.at(j);
					x2.at(j) = position.at(j).at(count) + 0.5*dt*v2.at(j) + 0.5*0.5*dt*dt*k2.at(j)/2;
				}

		double DOF3[7] = {x2.at(1), v2.at(1), x2.at(0), v2.at(0), 1, xGround, vGround };

			// %%%%%%%%%%%%%  THIRD SUB STEP %%%%%%%%%%%%%%
				for(int j=0; j<2; j++){
					for(int i=0; i<7; i++){
						k3.at(j) = k3.at(j) + AEM[i][j]*DOF3[i] ;
					}
					v3.at(j) = vel.at(j).at(count) + 0.5*dt*k3.at(j);
					x3.at(j) = position.at(j).at(count) + 0.5*dt*v3.at(j) + 0.5*0.5*dt*dt*k3.at(j)/2;
				}

		double DOF4[7] = {x3.at(1), v3.at(1), x3.at(0), v3.at(0), 1, xGround, vGround };

			// %%%%%%%%%%%%%  FOURTH SUB STEP %%%%%%%%%%%%%%
				for(int j=0; j<2; j++){
					for(int i=0; i<7; i++){
						k4.at(j) = k4.at(j) + AEM[i][j]*DOF4[i] ;
					}
				}

			// Assembling the RK-4th order method
				for(int i=0; i<2; i++){
					 v_n1.push_back( vel.at(i).at(count) + dt*(k1.at(i) + 2.0*k2.at(i) + 2.0*k3.at(i) + k4.at(i))/6.0 );
			 	 	 acc_n.push_back( (k1.at(i) + 2.0*k2.at(i) + 2.0*k3.at(i) + k4.at(i))/6.0 );
				}

				/*for(int i=0; i<2; i++){
					std::cout<<"Calculated Velocity of Body :: "<<i<<"at time ::"<<simTime<<" is :: "<<v_n1.at(i)<<std::endl;
					std::cout<<"Calculated accelera of Body :: "<<i<<"at time ::"<<simTime<<" is :: "<<acc_n.at(i)<<std::endl;
				}*/
				//getchar();


		// Updating the values in the vectors
		std::vector <double>diff_acc;
		for(int i=0; i<2; i++){
			acc.at(i).push_back(acc_n.at(i));		// Updating acceleration Vector
			vel.at(i).push_back(v_n1.at(i));		// Updating Velocity vector
			diff_acc.push_back( (acc.at(i).at(count+1) - acc.at(i).at(count))/dt );
			// Updating Position vector
			position.at(i).push_back( position.at(i).at(count) + dt*(vel.at(i).at(count) + 2.0*v1.at(i) +2.0*v2.at(i) + v3.at(i))/6  );
			// Updating Force
			force.at(i).push_back(-acc_n.at(i)*mass.at(i));

			}

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
void getGroundProps(double time, double stepHeight, double* displacement, double* velocity){

	double stepLength = 0.05;

	if(time < stepLength){
		*displacement = (stepHeight/0.05/1000)*time;
		*velocity = (stepHeight/0.05/1000);
	} else {
		*displacement = stepHeight/1000;
		*velocity = 0;
	}
}

