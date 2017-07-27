/*
 * main.cpp
 *
 *  Created on: Jan 29, 2013
 *      Author: aditya
 */
#include "body.hpp"
#include "joint.hpp"
#include "numerics.hpp"
#include "field.hpp"
#include <stdio.h>
#include <iostream>
#include <string>
#include <cmath>
#include <vector>

#define dataType double
#define pi (3.1416/180)

/**
 * Method for getting the initial conditions for the body 2
 * We only activate body 2 initially.
 */
void getInitialConditions() {

}

// TODO :: 1.	Formulate a way to formulate joint list in each body --> COMPLETED
// TODO :: 2.	Formulate methods to calculate forces and Torques on each body by iterating on the joints
// TODO :: 3.	Formulate methods to integrate Accelerations to get displacements and thetas
// TODO :: 4.	Write a solver to solve the system (diagonal) for obtaining accelerations.

// TODO :: 5.	points 1-4 should be done in numerics class.

template<class T> size_t body<T>::globalNumber = 0;
template<class T> size_t joint<T>::globalNumber = 0;

int main() {

	std::vector<body<dataType> > bodies;
	std::vector<joint<dataType> > joints;
	// Creating bodies
	field<dataType> point(1, 2);
	point[0] = 12.0;
	point[1] = 5.0;
	body<dataType> body1(10, point, 10*pi); // mass, position
	bodies.push_back(body1);

	point[0] = 10.5;
	point[1] = 24.0;
	body<dataType> body2(15, point, 45*pi);
	bodies.push_back(body2);

	point[0] = 31.5;
	point[1] = 18.3;
	body<dataType> body3(8, point, 90*pi);
	bodies.push_back(body3);

	point[0] = 23.5;
	point[1] = 12.3;
	body<dataType> body4(25, point, 0*pi);
	bodies.push_back(body4);

	// Creating Joints
	joint<dataType> joint1(bodies, 1, 1, 2, 1);
	joint<dataType> joint2(bodies, 2, 2, 3, 1);
	joint<dataType> joint3(bodies, 3, 2, 4, 1);
	joint<dataType> joint4(bodies, 4, 2, 1, 2);
	joint<dataType> joint5(bodies, 4, 2, 2, 1);
	joints.push_back(joint1); joints.push_back(joint2);
	joints.push_back(joint3); joints.push_back(joint4);
	joints.push_back(joint5);
	std::cout<<"SHIT"<<std::endl;

	numerics<dataType> Numerics;
	Numerics.formulateForces(bodies,joints);



	std::cout << "Total Number of Bodies :: " << body<dataType>::globalNumber
			<< std::endl;
	std::cout << "Total Number of Joints :: " << joint<dataType>::globalNumber
			<< std::endl;

	std::cout<<"Joint 4 :: x :: "<<joints.at(3).directionVector->operator [](0)<<" y :: "<<joints.at(3).directionVector->operator [](1)<<std::endl;

	return 0;
}

