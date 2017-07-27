/*
 * Calculator.cpp
 *
 *  Created on: Nov 8, 2015
 *      Author: friedrich
 */

#include "Calculator.h"

Calculator::Calculator(double x, double y): _x(x), _y(y) {
	// TODO Auto-generated constructor stub

}

Calculator::~Calculator() {
	// TODO Auto-generated destructor stub
}

double Calculator::add() const {
	return _x-_y;
}

double Calculator::minus() const {
	return _x-_y;
}

double Calculator::times() const {
	return _x*_y;
}

double Calculator::divide() const{
	return _x/_y;
}
