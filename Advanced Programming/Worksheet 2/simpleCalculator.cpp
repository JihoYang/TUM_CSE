#include <iostream>
#include <cmath>
#include <vector>

/*
 * Implementation of class SimpleCalculator
 */
class SimpleCalculator {
public:
	/*
	 * sum function
	 */
	double sum(double a, double b) {
		return a + b;
	}

	/*
	 * subtract function
	 */
	double subtract(double a, double b) {
		return a - b;
	}

	/*
	 * multiply function
	 */
	double multiply(double a, double b) {
		return a * b;
	}

	/*
	 * divide function
	 */
	double divide(double a, double b, bool& isDenominatorZero) {
		if (b == 0.0) {					// Check if b is zero. If yes, ..
			isDenominatorZero = true;		// .. inform computeResult about this,
			return HUGE_VAL;				// .. return a HUGE Value.
		} else {						// If b is not zero, ..
			return a / b;					// .. return the division result.
		}
	}

	/*
	 * mean function
	 */
	double mean(double a, double b) {
		return 0.5 * (a + b);
	}

	/*
	 * Compute result
	 */
	double computeResult(char oper, double numberA, double numberB) {
		double result;
		if (oper == '+') {
			result = sum(numberA, numberB);
		} else if (oper == '-') {
			result = subtract(numberA, numberB);
		} else if (oper == '*') {
			result = multiply(numberA, numberB);
		} else if (oper == '/') {
			bool isDenominatorZero = false;
			result = divide(numberA, numberB, isDenominatorZero);
			if (isDenominatorZero) {												// If denominator is zero, ..
				std::cout << "ERROR: Division by zero!" << std::endl;					// .. don't append result.
				result = 0;
			}
		} else if (oper == 'm') {
			result = mean(numberA, numberB);
		} else {																	// If operator is invalid, ..
			std::cout << std::endl << "ERROR: Invalid operator!" << std::endl;			// .. don't do anything.
		}

		return result;
	}

	/*
	 * worker function
	 */
	void work() {
		double numberA, numberB;	// Numbers a and b
		char oper;					// Operation to be performed

		std::cout << "Welcome to A Simple Calculator!" << std::endl;			// Politeness.

		std::cout << std::endl << "Enter the operation (+ | - | * | / | s | c | t): ";
		std::cin >> oper;														// .. read operator,

		std::cout << "Enter number A: ";										// .. read a,
		std::cin >> numberA;

		std::cout << "Enter number B: ";										// .. read b,
		std::cin >> numberB;

		double answer = computeResult(oper, numberA, numberB);					// .. compute result,
		std::cout << numberA << " " << oper << " " << numberB << " = " << answer << std::endl;

		std::cout << std::endl << "Ciao!" << std::endl;							// Politeness.
	}
};

int main() {
	SimpleCalculator SC1;			// Create instance of SimpleCalculator.

	SC1.work();						// Ask the caluclator to work.

	return 0;						// Successful execution of program.
}
