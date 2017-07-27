/*
 * Assignment2_main.cpp
 *
 *  Created on: Nov 9, 2015
 *      Author: erik
 */


#include <iostream>
#include <limits>
#include <ostream>

/**
 * function for finding the limit delta where base stays the same after adding
 * and subtracting delta from it
 *
 * @param base : the number in question
 * @return delta: largest number such that base == (base + delta) - delta
 */
float findUpperLimitAdditionFloat(float base){

	/* --------------------INITIALISATION------------------- */

	//Initialise an upper and a lower bound to the delta value. We assume one
	//can add base to itself.
	float lowerBound = base;
	float upperBound = base;

	/* ------------------END-INITIALISATION----------------- */

	//(What will happen if base is really, really, big?)

	/*					*		*		*					 */

	/* -------------FINDING DELTA (APPROXIMATELY)----------- */
	//Run a loop checking with loop condition that before each iteration
	//we can still add and subtract delta to base.

	//Use the upper bound for checking! When the loop finishes, upperBound will
	//be larger than the delta we want to find.
	while (base == ((base + upperBound) - upperBound))
	{

		// If we can still add and subtract upperBound, the delta we want to
		// find is still bigger than upperBound - it's a lower bound for delta.
		// Therefore, set lowerBound to the value of upperBound.
		lowerBound = upperBound;

		// Try in the next iteration with twice as big upper bound!
		upperBound = upperBound * 2.0f;
	}
	//Loop finished: we now have a lower and an upper bound for delta.


	/* -------------END FINDING DELTA (APPROXIMATELY)-------- */


	/*					*		*		*					  */


	/* --------------FINDING DELTA (MORE EXACTLY)------------ */
	//Find a better approximation to the delta by checking at the midpoint
	//between the lower and upper bound multiple times and closing in on the
	//exact value. (bisection)

	//Let's do this multiple times!
	static const int numberOfIterations = 100;
	for (int i = 0; i < numberOfIterations; i++) {

		//Start: finding the middle point between lowerBound and upperBound.
		float m = (lowerBound + upperBound) / 2.0f;

		//Check if the middle point can still be added and subtracted from base.
		if (base == (base + m) - m ) {

			//If the midpoint can be added and subtracted, then the it is a new
			//lower bound for delta.
			lowerBound = m;
		} else {

			//If the midpoint cannot be added and subtracted, then the it is
			//instead a new upper bound for delta.
			upperBound = m;
		}
	}
	/* ----------END FINDING DELTA (MORE EXACTLY)------------ */


	//lowerBound can still be added and subtracted from base, but should be
	//very close to the limit, so we return lowerBound.
	return lowerBound;
}

void printFormattedOutput(double base, double delta)
{
	std::cout	<< "Base:\t\t" 				<< base << "\t\t"
				<< "Biggest delta : \t\t "	<< delta << std::endl;
}

int main(int argc, char **argv) {

	//Calculate the maximum delta that can be added and subtracted from a base
	//number, start 10 ^ (6) and dividing by 10 each iteration, until 10^(-6).
	float base = 1e6f; // = 1 * 10^(6) 											(what happens if we use a power of two?)
	for (int i = 6; i >= -6; --i)
	{
		//Call our findUpperLimitAddition function to find delta
		float delta = findUpperLimitAdditionFloat(base);

		//Print the output
		printFormattedOutput(base,delta);

		//Divide base by 10 for the next iteration
		base /= 10.0;
	}



}



