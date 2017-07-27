/*
 * numerics.hpp
 *
 *  Created on: Jan 29, 2013
 *      Author: aditya
 */

#ifndef NUMERICS_HPP_
#define NUMERICS_HPP_

#include <cstdlib>
#include <cassert>
#include <vector>
#include <cmath>
#include <stdio.h>
#include <iostream>
#include "field.hpp"
#include "body.hpp"
#include "joint.hpp"

template<typename T>
class numerics;

template<typename T>
class numerics {
public:

	// Method to formulate Forces on the body.
	void formulateForces(std::vector< body<T> > &bodyVect, std::vector< joint<T> > &jointVect) {

		// Formulating vectors for joints
		for(int i=0; i<(int)jointVect.size(); i++){

			// Calculating the direction vectors of the joint
			int body1 = jointVect.at(i).connectingBody.at(0) -1;
			int body2 = jointVect.at(i).connectingBody.at(1) -1;

			int body1_point = jointVect.at(i).pointOnBody.at(0) -1;
			int body2_point = jointVect.at(i).pointOnBody.at(1) -1;

			T x1 = bodyVect.at( body1 ).points.at( body1_point ).back()[x];
			T y1 = bodyVect.at( body1 ).points.at( body1_point ).back()[y];

			T x2 = bodyVect.at( body2 ).points.at( body2_point ).back()[x];
			T y2 = bodyVect.at( body2 ).points.at( body2_point ).back()[y];

			T _x = x1-x2;
			T _y = y1-y2;
			T det = sqrt(_x*_x + _y*_y);

			jointVect.at(i).directionVector = new field<T> (1,2);
			jointVect.at(i).directionVector->operator[](0) = _x/1.0;
			jointVect.at(i).directionVector->operator[](1) = _y/1.0;

		}
		// Calculating the Forces
		for(int i=0; i<(int)bodyVect.size(); i++){

			//bodyVect.at(i).formulatePoints();
			bodyVect.at(i).formulatePointVelocities();
			// Loop over all the joints connected to the body to
			// Find the force on the body
			field<T> totalForce(1,2);
			totalForce[x] = 0.0;
			totalForce[y] = 0.0;

			for(int j=0; j<(int)bodyVect.at(i).joints.size(); j++){
				int jointNum = bodyVect.at(i).joints.at(j);
				T jointStiffness = jointVect.at( jointNum ).getStiffness();
				// Spring force
				field<T> springForce(1,2);
				T jX = jointVect.at(jointNum).directionVector->operator[](0);
				T jY = jointVect.at(jointNum).directionVector->operator[](0);
				T jointLength = sqrt( jX*jX + jY*jY );
				T force = jointLength*jointStiffness;

				springForce[x] = force*jX/jointLength;
				springForce[y] = force*jY/jointLength;

				// Damper force
				field<T> dampingForce(1,2);
				T jointDamping =  jointVect.at( jointNum ).getDamping(); // Damping Coefficient of joint
				T jVx = bodyVect.pointsVelocity.back().at(x);	// X component of Velocity
				T jVy = bodyVect.pointsVelocity.back().at(y);	// Y component of Velocity
				dampingForce[x] = jVx*jointDamping;				// X component of Damping Force
				dampingForce[y] = jVy*jointDamping;				// Y component of Damping Force

				// Updating the force
				totalForce[x] += springForce[x] + dampingForce[x];
				totalForce[y] += springForce[y] + dampingForce[y];
			}

			bodyVect.at(i).globalForces = new field(1,2);
		}
	}

};



#endif /* NUMERICS_HPP_ */
