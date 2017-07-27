/*
 * joint.hpp
 *
 *  Created on: Jan 29, 2013
 *      Author: aditya
 */

#ifndef JOINT_HPP_
#define JOINT_HPP_

#include <cstdlib>
#include <cassert>
#include <vector>

template<typename T>
class joint;

template<typename T>
class joint {

private:
	T stiffness; // Stiffness co-efficient
	T damping; // Damping co-efficient
	T freeLength; // Free length of the spring

public:
	static size_t globalNumber; // Number of the joints created
	size_t number;
	std::vector<size_t> connectingBody; // Contains the number of bodies it connects
	std::vector<size_t> pointOnBody; // Contains the number of points on the body it connects simulations.

	field<T>* directionVector;

	// Properties of the gap for this joint


	/**
	 * @param i_body1 	- number of 1st body the joint is attached to
	 * @param i_pOnBody1- number of point on the first body the joint is attached to.
	 * @param i_body2	- number of 2nd body the joint is attached to
	 * @param i_pOnBody2- number of point on the second body the joint is attached to.
	 */
	joint(std::vector< body<T> >& bodies, size_t i_body1, size_t i_pOnBody1, size_t i_body2,
			size_t i_pOnBody2, T i_freeLength = 5.0) {
		number = globalNumber++;

		connectingBody.push_back(i_body1);
		connectingBody.push_back(i_body2);
		pointOnBody.push_back(i_pOnBody1);
		pointOnBody.push_back(i_pOnBody2);

		freeLength = i_freeLength;

		directionVector = new field<T>(1,2);

		// Formulating the joint list in corresponding bodies which this joint connects
		bodies.at(i_body1-1).joints.push_back(this->number);
		bodies.at(i_body2-1).joints.push_back(this->number);
	}

	void setStiffness(T i_stiffness) {
		stiffness = i_stiffness;
	}
	void setDampint(T i_damping) {
		damping = i_damping;
	}

	T getStiffness() {
		return stiffness;
	}
	T getDamping() {
		return damping;
	}

};

#endif /* JOINT_HPP_ */
