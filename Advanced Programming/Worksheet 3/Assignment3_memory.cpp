/*
 * Assignment3_memory.cpp
 *
 *  Created on: Nov 5, 2015
 *      Author: erik
 */


#include <iostream>
#include <string>



int address_alignment(void * pointer) {

	// interpret pointer as an integer to do arithmetics on the address.
	// store this in a size_t, (an unsigned integer type defined to be large
	// enough to store adresses on our system)
	size_t address = reinterpret_cast<size_t>(pointer);

	//to check which powers of two variables are aligned to, we try to divide
	//it by higher and higher multiples of two until we get a nonzero remainder.
	//The pointer is then aligned with an address that is a multiple of the
	//highest divisor.

	size_t divisor = 1; //initial value: every integer is dividable by 1.
	size_t powerOfTwoOfDivisor = 0; //2^0 = 1.


	//check for nonzero remainder of the next power of two using the modulo
	//operator (% in c++)
	while (address % (2*divisor) == 0) {
		//address was dividable by the next power of two, so
		//try next power of two by multiplying by two
		divisor *= 2;

		// also increment the counter by one!
		powerOfTwoOfDivisor += 1;

	}

	//We have the largest power of two; let's return it.
	return powerOfTwoOfDivisor;
}



//-----------------------------------------------------------------------------


//Simple shorthand for calculating the power of two.
int pow2(int k) {
  int n = 1;
  for (int i = 0; i < k; i++) {
    n *= 2;
  }
  return n;
}


//-----------------------------------------------------------------------------

template <class TYPE>
void printVariableAndInfo(TYPE& variableToPrint, std::string variableTypeName)
{
	//Call functions to get data about our variable.
	void * addressAsVoidPtr = static_cast <void*> (&variableToPrint);
	int alignedToPowerOfTwo = address_alignment(addressAsVoidPtr);

	//Print variable info
	std::cout << "Variable type:\t" << variableTypeName << "\t\t";
	std::cout << "Value:\t\t" << variableToPrint ;
	std::cout << "\t\t(Value as int: \t";
	std::cout << static_cast<int> (variableToPrint) << ")\n";

	//Print other data
	std::cout << "Size: \t" <<sizeof(TYPE) << "\t\t\t";
	std::cout << "Address: \t" <<addressAsVoidPtr << "\t";
	std::cout << "Aligned to: 2^" << alignedToPowerOfTwo << " ( = "
			<< pow2(alignedToPowerOfTwo) << ")\n\n";
}



//-----------------------------------------------------------------------------


void declareAndPrintVariablesStack() {
	char c = 1;
	char c2 = 2;
	double d = 3.0;
	int i = 3;


	printVariableAndInfo(c,	"char");
	printVariableAndInfo(c2,	"char");
	printVariableAndInfo(d,	"double");
	printVariableAndInfo(i,	"int");

}


//-----------------------------------------------------------------------------


void declareAndPrintVariablesHeap() {

	int * i = new int(0);
	char * cs = new char [3] {1,2,3};
	double * d0 = new double(3.142);
	double * d1 = new double(1337);

	printVariableAndInfo(*i,	"int");
	printVariableAndInfo(cs[0],	"char0");
	printVariableAndInfo(cs[1],	"char1");
	printVariableAndInfo(cs[2],	"char2");
	printVariableAndInfo(*d0,	"double0");
	printVariableAndInfo(*d1,	"double1");



	delete i;
	delete [] cs;
	delete d0;
	delete d1;
}


//-----------------------------------------------------------------------------

int main(int argc, char **argv) {

	std::cout << std::endl
			<< "Creating variables on the stack! (inside function)"
			<< std::endl;

	//declare and print some random variables on the stack.
	declareAndPrintVariablesStack();


	//Do the same thing for some variables on the heap.
	std::cout << std::endl
			<< "Creating variables on the heap! (using new) "
			<< std::endl;

	declareAndPrintVariablesHeap();

	//You can also use compiler directives to tell the compiler to pack
	//variables as tightly as possible. #pragma pack(push,n) tells the compiler
	//to start packing as tightly as possible from now on aligning at multiples
	//of n bytes, stopping at #pragma pack(pop).

//	std::cout << "Creating a weird packed structure" << std::endl;
//#pragma pack(push,1)
//	struct {
//		int b;
//		char a;
//		double c;
//	} packedStruct = {1, 2, 2.0};
//#pragma pack(pop)
//
//	printVariableAndInfo(packedStruct.b, "packed int");
//	printVariableAndInfo(packedStruct.a, "packed char");
//	printVariableAndInfo(packedStruct.c, "packed double");


  return 0;
}




