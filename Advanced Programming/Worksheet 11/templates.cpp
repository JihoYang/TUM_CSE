#include <iostream>
#include <cstdlib>
#include <string>

//#define SQUAREMACRO
#ifdef SQUAREMACRO
#define square(a) a*a
#endif

// TODO: ASSIGNMENT 1
/// a function template to compute the square for general numeric types (requrires existence of "operator*")
template<class T>
T square(T a) {
	std::cout << "using generic version to square " << a << std::endl;
	return a*a;
}

/// "int" is so special -- provide a specialized function template
template<>
int square(int a) {
	std::cout << "using specialized generic version to square " << a << std::endl;
	return a*a;
}

/// and there's one more: a non-template square function for int
int square(int a) {
	std::cout << "using explicit version to square " << a << std::endl;
	return a*a;
}

// TODO: ASSIGNMENT 2: implement TArray class
// When you are finished with the implementation activate 
// the following macro to enable the provided test routine.
//#define ASSIGNMENT2_DONE

// test stub implementation
int main(int argc, char **argv) {
    std::cout << "ASSIGNMENT 1:" << std::endl;

	int aInt     = 4;
	float aFloat = 4.2f;

	int aIntSquared = square(aInt + 1);
	std::cout << "int:   " << (aInt + 1) << "^2=" << aIntSquared << std::endl;

    float aFloatSquared = square(aFloat);
	std::cout << "float: " << aFloat << "^2=" << aFloatSquared << std::endl;

    // TODO: How can you call the specialized version of square for integers?    

    // ------------------------------------------------------------ //

#if defined(ASSIGNMENT2_DONE)
    std::cout << "ASSIGNMENT 2:" << std::endl;

    // allocate TArray object
    TArray<double> array_a(10);

    // check implementation of size()
    if (array_a.size() != 10) {
        std::cout << "size of array_a: " << array_a.size()
                  << " but actual size should be 10"
                  << std::endl;
    } else {
        std::cout << "size() works as expected" << std::endl;
    }

    // check implementation of array operator[]
    for (size_t i=0; i < array_a.size(); ++i) {
        array_a[i] = static_cast<double>(i+10);
    }
    bool all_correct = true;
    for (size_t i=0; i < array_a.size(); ++i) {
        if (array_a[i] != static_cast<double>(i+10)) {
            all_correct = false;
            std::cout << "value at position " << i << " is wrong: "
                      << " it should be " << static_cast<double>(i+10)
                      << std::endl;
        }
    }
    
    if (!all_correct) {
        std::cout << "operator[] returns the wrong values" << std::endl;
    } else {
        std::cout << "operator[] works as expected" << std::endl;
    }
    
    std::cout << "tests completed" << std::endl;
    std::cout << std::endl;
#endif

    // ------------------------------------------------------------ //

    return 0;
}
