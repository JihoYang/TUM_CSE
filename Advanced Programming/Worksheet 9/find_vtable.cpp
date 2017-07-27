#include <iostream>

/*******************************************************************************************
 * Compile without optimization (otherwise the compiler will remove the interesting parts)
 * 
 * g++ -O0 -g -o find_vtable find_vtable.cc 
 */

// ----------------------------------------------- //
// non virtual classes
// ----------------------------------------------- //
class EmptyClass {};

class SimpleClass {
    public:
        SimpleClass() { }
        ~SimpleClass() { }
};

// ----------------------------------------------- //
// virtual classes
// ----------------------------------------------- //
class DummyClass {
    public:
        DummyClass() { }
        virtual ~DummyClass() { }
};

class BaseClass {
    public:
        BaseClass() { }
        virtual void challenge() { 
            std::cout << "base class calling back" << std::endl; 
        }
        virtual ~BaseClass() { }
};

class Subclass : public BaseClass {
    public:
        Subclass() { }
        virtual ~Subclass() { }
};

class OverloadingSubclass : public BaseClass {
    public:
        OverloadingSubclass() { }
        virtual void challenge() {
            std::cout << "overloading subclass calling back" << std::endl;
            BaseClass::challenge();
        }
        virtual ~OverloadingSubclass() { }
};

void print_hidden_table(size_t tableaddr, size_t tablesize) {
    std::cout << "-------- table contents --------" << std::endl;
    // TODO: implement this helper function to print the contents of the hidden table (= fixed array with size_t blocks)
    // tableaddr: specifies the location of the table in memory
    // tablesize: specifies the size of this table in bytes
}

void call_hidden_table_entry(size_t tableaddr, size_t tableentry) {
    std::cout << "[CALL] function from table @ " << tableaddr << " and entry " << tableentry <<  std::endl;
  
    // TODO: implement this helper function to call a method which is given by the hidden table
    // tableaddr: specifies the location of the table in memory
    // tableentry: specifies the entry in the 

    std::cout << "[CALL DONE] function from table" << std::endl;
}

int main(int argc, char **argv) {

	std::cout << "***************************************************************************" << std::endl;
	std::cout << " Task 1:" << std::endl;
	std::cout << " Analysis of class sizes, dep. on whether they are virtual or non-virtual" << std::endl;
	std::cout << "***************************************************************************" << std::endl;

    // TODO: Print size of each predefined class type and discuss the differences.

	std::cout << "***************************************************************************" << std::endl;
	std::cout << " Task 2:" << std::endl;
	std::cout << " Analysis of strides in memory between classes of different type and size," << std::endl;
	std::cout << " also dep. on whether they are virtual or non-virtual" << std::endl;
	std::cout << "***************************************************************************" << std::endl;

    // TODO: Print the values of the hidden members of each class type.
    // HINT: First, use an instance of each predefined class and interpret it as simple chunk of memory.
    //       Then try to interpret this chunk of memory as a datatype which has the same size.
    
    {   // TODO: show hidden member of EmptyClass
        // ... this will probably print junk.
    
    }
 
    {   // TODO: show hidden member of SimpleClass
        // ... this will probably print junk.
    
    }
 
    {   // TODO: show hidden member of DummyClass
        // ... this should NOT print junk.
    
    }

    {   // TODO: show the hidden member of two BaseClass instances
        // ... this should NOT print junk. what is so special about this member? 
    
    }

    {   // TODO: show hidden member of Subclass
        // .... this should NOT print junk.
   
    }

    {   // TODO: show hidden member of OverloadingSubclass
        // ... this should NOT print junk.
   
    }


	std::cout << "***************************************************************************" << std::endl;
	std::cout << " Task 3:" << std::endl;
	std::cout << " Analysis of class contents, for both kinds of classes, virtual and non-virtual" << std::endl;
	std::cout << "***************************************************************************" << std::endl;

    // task 3: - what's behind the address? a function or something more complex? print the offset differences! (how does it change with the number of methods)
    //         - what is the purpose of this special hidden member variable, what's behind it and what does the offset difference indicate?
    // (also think about that these objects are created by the runtime, and the compiler still has to call the right functions)
    // The class types are arranged like a stack in "memory"
    {

        // TODO: - compute the difference of the hidden value of DummyClass and the hidden value of BaseClass.
        //       - implement and use the print_hidden_table function to display the hidden table
 
        // TODO: - compute the difference of the hidden value of Subclass and the hidden value of BaseClass
        //       - implement and use the print_hidden_table function to display the hidden table
 
        // TODO: - compute the difference of the hidden value of OverloadingSubclass and the hidden value of Subclass
        //       - implement and use the print_hidden_table function to display the hidden table
    
    }

    return 0;
}
