#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>

// TODO: ASSIGNMENT 1: (part 1)
// implement selection sort for a double array
void mysort(std::vector<double>& array) {

}

// TODO: ASSIGNMENT 1: (part 2)
// convert function to a templated function.

//#define SORT_DONE

// TODO: ASSIGNMENT 2: sort array of pointers in a special way
//#define ASSIGNMENT2_DONE

// test stub implementation
int main(int argc, char **argv) {
#if defined(SORT_DONE)
    std::cout << "ASSIGNMENT 1:" << std::endl;
    
    // seed some random values
    std::vector<double> array_a(10);
    for (size_t i=0; i < array_a.size(); ++i) {
        array_a[i] = static_cast<double>(rand())/RAND_MAX;
    }

    // print unsorted array
    std::cout << "unsorted:" << std::endl;
    for (size_t i=0; i < array_a.size(); ++i) {
        std::cout << " " << array_a[i];
    }
    std::cout << std::endl;

    mysort(array_a);

    // print sorted array
    std::cout << "sorted:" << std::endl;
    for (size_t i=0; i < array_a.size(); ++i) {
        std::cout << " " << array_a[i];
    }
    std::cout << std::endl;

    std::cout << std::endl;
#endif

    // ------------------------------------------------------------ //

#if defined(ASSIGNMENT2_DONE)

    std::cout << "ASSIGNMENT 2:" << std::endl;

    // now we allocate an array of double pointers ....
    std::vector<double*> array_ptrs(10);
    for (size_t i=0; i < array_ptrs.size(); ++i) {
        array_ptrs[i] = new double();
        *(array_ptrs[i]) = static_cast<double>(rand())/RAND_MAX;
    }

    // ... and try to sort them
    std::cout << "unsorted:" << std::endl;
    for (size_t i=0; i < array_ptrs.size(); ++i) {
        std::cout << " " << *(array_ptrs[i]);
    }
    std::cout << std::endl;
  
    mysort(array_ptrs);

    std::cout << "sorted:" << std::endl;
    for (size_t i=0; i < array_ptrs.size(); ++i) {
        std::cout << " " << *(array_ptrs[i]);
    }
    std::cout << std::endl;

    // finally, we cleanup some pointers
    for (size_t i=0; i < array_ptrs.size(); ++i) {
        delete (array_ptrs[i]);
    }

    std::cout << std::endl;

#endif

    // ------------------------------------------------------------ //

    return 0;
}
