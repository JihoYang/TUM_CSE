//============================================================================
// Name        : Unittesting.cpp
// Author      : FMenhorn
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <cppunit/ui/text/TestRunner.h>

#include "Calulator_suite.h"
using namespace std;

int main() {
	CppUnit::TextUi::TestRunner runner;
	runner.addTest(Calulator_suite::suite());
	runner.run();

	return 0;
}
