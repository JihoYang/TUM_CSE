/*
 * body.hpp
 *
 *  Created on: Jan 29, 2013
 *      Author: aditya
 */

#ifndef BODY_HPP_
#define BODY_HPP_

#include <cstdlib>
#include <cassert>
#include <vector>
#include <cmath>
#include <stdio.h>
#include <iostream>
#include "field.hpp"

#define x 0
#define y 1

template<typename T>
class body;

template<typename T>
class body {
private:

	size_t number;
	T length; // Length of the body for calculating I
	T width; // Width of the body for calculating I
	T inertia; // Moment of inertia of the body.
	T mass;
	// The distance of the point from center of mass.
	// One should be negative and one should be positive -> ideally
	// Origin of body co-ordinate system coincides to C.M
	std::vector<T> distPoint;

	void initialize(){
		field<T> vel(1, 2);
		vel[0] = 0.0;
		vel[1] = 0.0;
		velocity.push_back(vel);

		field<T> acc(1, 2);
		acc[0] = 0.0;
		acc[1] = 0.0;
		acceleration.push_back(acc);

	}

public:
	static size_t globalNumber; // Number of the bodies created
	std::vector<int> joints;
	// Initial value is always set to zero;
	std::vector<T> theta; // Contains the theta of the body w.r.t x axis -> including the initial value.

	std::vector< std::vector< field<T> > > points; // First point coordinates where a joint can be attached

	std::vector< std::vector< field<T> > > pointsVelocity; // Velocity of the points

	std::vector<field<T> > velocity; // Velocity of the body. its a vector (u,v)
	std::vector<field<T> > acceleration; // Acceleration of body. Its a vector (ax, ay)
	std::vector<T> alpha; // Angular Acceleration of body about its C.M. Its a Scalar
	std::vector<T> omega; // Angular velocity of the body about its C.M. Its a Scalar
	field<T> *globalForces; // These forces are always in the global Coordinate System
	field<T> *localForces; // These forces are always in the local Coordinate System

	T torque;

	std::vector<field<T> > massCenter; // Contains the position of C.M of the body --> including the initial value.

	// Constructor
	body(T mass, field<T> i_centerMass, T i_theta = 0.0, T i_length = 1.0,
			T i_width = 0.5) {

		number = globalNumber++;
		length = i_length;
		width = i_width;
		inertia = length * pow(width, 3) / 12;
		theta.push_back(i_theta);
		massCenter.push_back(i_centerMass);
		omega.push_back(0.0);

		distPoint.push_back(0.9);
		distPoint.push_back(-0.9);
		initialize(); /// Initializes the velocity and acceleration to zero
		std::vector< field<T> > dummy;
		points.push_back(dummy);
		points.push_back(dummy);
		pointsVelocity.push_back(dummy);
		pointsVelocity.push_back(dummy);
		formulatePoints();
	}


	// Getters and setters
	size_t getNumber() {
		return number;
	}

	T getInertia() {
		return inertia;
	}

	T getPoint1Dist() {
		return distPoint.at(0);
	}

	T getPoint2Dist() {
		return distPoint.at(1);
	}

	// Method to formulate the points of attachment of joints
	void formulatePoints() {

		field<T> point_1(1, 2); // Dummy variable to store the calculated coordinates

		// Formulating Point 1
		point_1[x] = distPoint.at(0) * cos(theta.back()) + massCenter.back()[x];
		point_1[y] = distPoint.at(0) * sin(theta.back()) + massCenter.back()[y];
		points.at(0).push_back(point_1);

		field<T> point_2(1, 2); // Dummy variable to store the calculated coordinates
		// Formulating Point 2
		point_2[x] = distPoint.at(1) * cos(theta.back()) + massCenter.back()[x];
		point_2[y] = distPoint.at(1) * sin(theta.back()) + massCenter.back()[y];
		points.at(1).push_back(point_2);

		 std::cout<<"Body "<<this->number<<" Point 1 :: "<<points.at(0).back()[x]<<"  "<<points.at(0).back()[y]<<std::endl;
		 std::cout<<"Body "<<this->number<<" Point 2 :: "<<points.at(1).back()[x]<<"  "<<points.at(1).back()[y]<<std::endl;

	}

	void formulatePointVelocities() {
		field<T> velocity_1(1, 2);
		// X component
		velocity_1[x] = velocity.back()[x] + distPoint.at(0)*omega.back()*sin(theta.back());
		velocity_1[y] = velocity.back()[y] + distPoint.at(0)*omega.back()*cos(theta.back());
		pointsVelocity.at(0).push_back(velocity_1);

		field<T> velocity_2(1, 2);
		// Y Component
		velocity_2[x] = velocity.back()[x] + distPoint.at(1)*omega.back()*sin(theta.back());
		velocity_2[y] = velocity.back()[y] + distPoint.at(1)*omega.back()*cos(theta.back());
		pointsVelocity.at(1).push_back(velocity_2);
	}
};

#endif /* BODY_HPP_ */
