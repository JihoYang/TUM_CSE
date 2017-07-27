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
    tablesize /= sizeof(size_t);
    size_t *hiddentable = reinterpret_cast<size_t*>(tableaddr);
    for (size_t element=0; element < tablesize; ++element) {
        std::cout << "hidden[" << std::dec << element << "]: " << std::hex << hiddentable[element] << std::endl;
    }
}

void call_hidden_table_entry(size_t tableaddr, size_t tableentry) {
    std::cout << "[CALL] function from table @ " << tableaddr << " and entry " << tableentry <<  std::endl;
  
    // TODO: implement this helper function to call a method which is given by the hidden table
    // tableaddr: specifies the location of the table in memory
    // tableentry: specifies the entry in the 
    size_t *hiddentable = reinterpret_cast<size_t*>(tableaddr);
    void (*fp)(void) = reinterpret_cast<void (*)(void)>(hiddentable[tableentry]);
    fp();

    std::cout << "[CALL DONE] function from table" << std::endl;
}

int main(int argc, char **argv) {

	std::cout << "***************************************************************************" << std::endl;
	std::cout << " Task 1:" << std::endl;
	std::cout << " Analysis of class sizes, dep. on whether they are virtual or non-virtual" << std::endl;
	std::cout << "***************************************************************************" << std::endl;

    // TODO: Print size of each predefined class type and discuss the differences.
    std::cout << "size of EmptyClass: " << sizeof(EmptyClass) << std::endl;
    std::cout << "size of SimpleClass: " << sizeof(SimpleClass) << std::endl;
    std::cout << "size of DummyClass: " << sizeof(DummyClass) << std::endl;
    std::cout << "size of BaseClass: " << sizeof(BaseClass) << std::endl;


	std::cout << "***************************************************************************" << std::endl;
	std::cout << " Task 2:" << std::endl;
	std::cout << " Analysis of strides in memory between classes of different type and size," << std::endl;
	std::cout << " also dep. on whether they are virtual or non-virtual" << std::endl;
	std::cout << "***************************************************************************" << std::endl;

    // TODO: Print the values of the hidden members of each class type.
    // HINT: First, use an instance of each predefined class and interpret it as a simple chunk of memory.
    //       Then try to interpret this chunk of memory as a datatype which has the same size.
    
    {   // TODO: show hidden member of EmptyClass
        // ... this will probably print junk.
        EmptyClass empty;
        unsigned char value_of_hidden_member = *(reinterpret_cast<unsigned char*>(&empty));
        std::cout << "hidden value of empty (EmptyClass): " << std::hex << (static_cast<unsigned char>(value_of_hidden_member) & 0xFF) << std::endl;
    }
 
    {   // TODO: show hidden member of SimpleClass
        // ... this will probably print junk.
        SimpleClass simple;
        unsigned char value_of_hidden_member = *(reinterpret_cast<unsigned char*>(&simple));
        std::cout << "hidden value of simple (SimpleClass): " << std::hex << (static_cast<unsigned char>(value_of_hidden_member) & 0xFF) << std::endl;
    }
 
    {   // TODO: show hidden member of DummyClass
        // ... this should NOT print junk.
        DummyClass dummy;
        size_t value_of_hidden_member = *(reinterpret_cast<size_t*>(&dummy));
        std::cout << "hidden value of dummy (DummyClass): " << std::hex << static_cast<unsigned int>(value_of_hidden_member) << std::endl;
    }

    {   // TODO: show the hidden member of two BaseClass instances
        // ... this should NOT print junk. what is so special about this member? 
        BaseClass base1;
        size_t value_of_hidden_member = *(reinterpret_cast<size_t*>(&base1));
        std::cout << "hidden value of base 1 (BaseClass): " << std::hex << static_cast<unsigned int>(value_of_hidden_member) << std::endl;
 
        BaseClass base2;
        size_t value_of_hidden_member2 = *(reinterpret_cast<size_t*>(&base2));
        std::cout << "hidden value of base 2 (BaseClass): " << std::hex << static_cast<unsigned int>(value_of_hidden_member2) << std::endl;
    }

    {   // TODO: show hidden member of Subclass
        // .... this should NOT print junk.
        Subclass subclass;
        size_t value_of_hidden_member = *(reinterpret_cast<size_t*>(&subclass));
        std::cout << "hidden value of subclass: " << std::hex << static_cast<unsigned int>(value_of_hidden_member) << std::endl;
    }

    {   // TODO: show hidden member of OverloadingSubclass
        // ... this should NOT print junk.
        OverloadingSubclass overloading;
        size_t value_of_hidden_member = *(reinterpret_cast<size_t*>(&overloading));
        std::cout << "hidden value of overloading subclass: " << std::hex << static_cast<unsigned int>(value_of_hidden_member) << std::endl;
    }


	std::cout << "***************************************************************************" << std::endl;
	std::cout << " Task 3:" << std::endl;
	std::cout << " Analysis of class contents, for both kinds of classes, virtual and non-virtual" << std::endl;
	std::cout << "***************************************************************************" << std::endl;

    // task 3: - what's behind the address? a function or something more complex? print the offset differences! (how does it change with the number of methods)
    //         - what is the purpose of this special hidden member variable, what's behind it and what does the offset difference indiciate?
    // (also think about that these objects are created by the runtime, and the compiler still has to call the right functions)
    // The class types are arranged like a stack in "memory"
    {

        // TODO: - compute the difference of the hidden value of DummyClass and the hidden value of BaseClass.
        //       - implement and use the print_hidden_table function to display the hidden table
        DummyClass dummy;
        size_t dummyclass_hidden_member = *(reinterpret_cast<size_t*>(&dummy));
 
        BaseClass base;
        size_t baseclass_hidden_member = *(reinterpret_cast<size_t*>(&base));
        {
            size_t difference = dummyclass_hidden_member - baseclass_hidden_member;
            std::cout << "difference between base and dummy: " << std::dec << difference << std::endl;
            print_hidden_table(baseclass_hidden_member, difference);
        }
 
        // TODO: - compute the difference of the hidden value of Subclass and the hidden value of BaseClass
        //       - implement and use the print_hidden_table function to display the hidden table
        Subclass Subclass;
        size_t subclass_hidden_member = *(reinterpret_cast<size_t*>(&Subclass));
        {
            size_t difference = baseclass_hidden_member - subclass_hidden_member;
            std::cout << "difference between Subclass and BaseClass: " << std::dec << difference << std::endl;
            print_hidden_table(subclass_hidden_member, difference); // NOW: check with BaseClass and notice similarities!
        }
 
        // TODO: - compute the difference of the hidden value of OverloadingSubclass and the hidden value of Subclass
        //       - implement and use the print_hidden_table function to display the hidden table
        OverloadingSubclass overloading;
        size_t overloadingsubclass_hidden_member = *(reinterpret_cast<size_t*>(&overloading));
        {
            size_t difference = subclass_hidden_member - overloadingsubclass_hidden_member;
            std::cout << "difference between OverloadingSubclass and Subclass: " << std::dec << difference << std::endl;
            print_hidden_table(overloadingsubclass_hidden_member, difference);
        }
    }

	std::cout << "***************************************************************************" << std::endl;
	std::cout << " Task 4: THE CHALLENGE" << std::endl;
	std::cout << " Low-level (and not recommended) way of calling virtual member functions   " << std::endl;
	std::cout << "***************************************************************************" << std::endl;
 
    // --------------------------------------------------------------------------------------------
    // IMPORTANT: - Look for the entries of the hidden tables from task 3 in the output of:
    //              objdump -d <yourbinary>
    //            - Use this information to determine the table entry for the "challenge()" method
    // --------------------------------------------------------------------------------------------

    {   // TODO: execute challenge() method of a BaseClass instance in a low-level way
        BaseClass base;
        size_t baseclass_hidden_member = *(reinterpret_cast<size_t*>(&base));
        call_hidden_table_entry(baseclass_hidden_member, 0); // what function are we calling here?
    }

    {   // TODO: execute challenge() method of a Subclass instance in a low-level way
        Subclass Subclass;
        size_t subclass_hidden_member = *(reinterpret_cast<size_t*>(&Subclass));
        call_hidden_table_entry(subclass_hidden_member, 0); // what function are we calling here?
    }

    {   // TODO: execute challenge() method of a OverloadingSubclass instance in a low-level way
        OverloadingSubclass overloading;
        size_t overloadingsubclass_hidden_member = *(reinterpret_cast<size_t*>(&overloading));
        call_hidden_table_entry(overloadingsubclass_hidden_member, 0); // what function are we calling here?
    }

    return 0;
}
