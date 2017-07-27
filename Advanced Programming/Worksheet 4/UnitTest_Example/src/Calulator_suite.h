/*
 * Calulatorsuite.h
 *
 *  Created on: Nov 8, 2015
 *      Author: friedrich
 */

#ifndef CALULATOR_SUITE_H_
#define CALULATOR_SUITE_H_

#include <cppunit/TestSuite.h>
/*
 *
 */
class Calulator_suite{
public:
	Calulator_suite();
	virtual ~Calulator_suite();

	static CppUnit::TestSuite* suite();
};

#endif /* CALULATOR_SUITE_H_ */
