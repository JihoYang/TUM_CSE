/*
 * body.hpp
 *
 *  Created on: Jan 4, 2013
 *      Author: aditya
 */

#ifndef BODY_HPP_
#define BODY_HPP_

#include <stdio.h>
#include <iostream>
#include <string>
#include <cmath>
#include <vector>


template <typename T>
class Cbody{

private:


public:

	// Mass of the body
	T mass;

	// Initial position of the body
	T iniPosition;

	// vector for storing the displacement of the body
	std::vector<T> position;

	// vector for storing the velocity of the body
	std::vector<T> vel;

	// Vector for storing the acceleration of the body
	std::vector<T> acc;

	// Vector for storing the force of the body
	std::vector<T> force;


	/**
	 * Constructor(s)
	 */
	Cbody(){
		mass = 0;
		iniPosition = 0;
		position.push_back((T)0.0);
	}

	Cbody(T m, T iniPosit){
		mass = m;
		iniPosition = iniPosit;
		position.push_back(iniPosit);
	}
};


#endif /* BODY_HPP_ */
