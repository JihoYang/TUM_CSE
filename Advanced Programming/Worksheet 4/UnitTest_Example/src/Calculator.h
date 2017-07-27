/*
 * Calculator.h
 *
 *  Created on: Nov 8, 2015
 *      Author: friedrich
 */

#ifndef CALCULATOR_H_
#define CALCULATOR_H_

/*
 *
 */
class Calculator {
private:
	double _x;
	double _y;
public:
	Calculator(double x, double y);
	virtual ~Calculator();

	double add() const;
	double minus() const;
	double times() const;
	double divide() const;
};

#endif /* CALCULATOR_H_ */
