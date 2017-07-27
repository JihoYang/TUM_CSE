/*
 * Calculator_test.cpp
 *
 *  Created on: Nov 8, 2015
 *      Author: friedrich
 */

#include "Calculator.h"
#include <cppunit/TestCase.h>
#include <assert.h>

class Calculator_test: public CppUnit::TestFixture{
private:
	Calculator* calculator;
public:
	void setUp(){
		double x = 4;
		double y = 2;
		calculator = new Calculator(x,y);
	}

	void tearDown(){
		delete calculator;
	}

	void testAdd(){
		double result = calculator->add();
		double expectedResult = 6.0;
		CPPUNIT_ASSERT_EQUAL(result,expectedResult);
	}

	void testMinus(){
		double result = calculator->minus();
		double expectedResult = 2.0;
		CPPUNIT_ASSERT_EQUAL(result,expectedResult);
	}

	void testTimes(){
		double result = calculator->times();
		double expectedResult = 8.0;
		CPPUNIT_ASSERT_EQUAL(result,expectedResult);
	}

	void testDivide(){
		double result = calculator->divide();
		double expectedResult = 2.0;
		CPPUNIT_ASSERT_EQUAL(result,expectedResult);
	}
};
