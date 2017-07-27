/*
 * Calulatorsuite.cpp
 *
 *  Created on: Nov 8, 2015
 *      Author: friedrich
 */

#include <cppunit/TestCaller.h>
#include "Calculator_test.cpp"
#include "Calulator_suite.h"

Calulator_suite::Calulator_suite() {
	// TODO Auto-generated constructor stub
}

Calulator_suite::~Calulator_suite() {
	// TODO Auto-generated destructor stub
}

CppUnit::TestSuite* Calulator_suite::suite() {
	CppUnit::TestSuite *calculatorSuite = new CppUnit::TestSuite( "calculatorSuite" );
	calculatorSuite->addTest(new CppUnit::TestCaller<Calculator_test>(
								   "testAdd",
								   &Calculator_test::testAdd ));
	calculatorSuite->addTest(new CppUnit::TestCaller<Calculator_test>(
									   "testMinus",
									   &Calculator_test::testMinus ));
	calculatorSuite->addTest(new CppUnit::TestCaller<Calculator_test>(
									   "testTimes",
									   &Calculator_test::testTimes ));
	calculatorSuite->addTest(new CppUnit::TestCaller<Calculator_test>(
									   "testDivide",
									   &Calculator_test::testDivide ));
	return calculatorSuite;
}
