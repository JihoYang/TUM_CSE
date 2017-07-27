#include <iostream>
#include <iomanip>
#include <typeinfo>

using namespace std;


/**
 * Template class implementing unique pointer semantics for arbitrary types
 */
template<class T>
class UniquePointer {
public:

    // TODO: complete the implementation
    // HINT: do not lose information.
    /// create a unique pointer to hold the instance passed
    UniquePointer(T* obj) {};

    // TODO: complete the implementation, implement the _semantics_ of a unique pointer!
    /// "move" the instance into a new instance of the unique pointer type
    UniquePointer(UniquePointer& uptr)
	{
	};

	/// forbid the default copy constructor explicitly using C++11 features
    UniquePointer(const UniquePointer& uptr) = delete;

    /// destructor that cleans up if necessary, also notifies about actions taken
    ~UniquePointer() {
    	// TODO complete the implementation
    };

    /// assignment operator for unique pointer type 
    // TODO: is it needed? is it clean?
    void
    operator=(UniquePointer& uptr)// = delete;
    {
	    // TODO: implement _clean_ assignment, implementing unique pointer semantics consistently
    };

    // TODO implement function(s) to give convenient access to the instance
 
    /// dereference to work with the instance
    T*
    operator->() { /*TODO:*/ };

    /// TODO: const dereference to work with instance the in const-contexts

private:

    /// the instance all this fuss is about
    T*     _obj;
};


/**
 * Simple test class that does not really matter
 */
struct TestClass
{
    TestClass() {
        std::cout << "Creating instance of TestClass" << std::endl;
    };

    ~TestClass() {
        std::cout << "Deleting instance of TestClass" << std::endl;
    };

    void
    doSomething() const {
        std::cout << "Doing something" << std::endl;
    };
};

int main() {

// TODO remove comment (after completing implementation of UniquePointer class), run the example and interpret
/*
    typedef UniquePointer<TestClass> PtrType;

    std::cout << "test case 1: deleting in inner context" << std::endl;
    {
        PtrType uptr1(new TestClass());
        {
            PtrType uptr2(uptr1);
            const PtrType& uptr2_cref = uptr2;

            //TODO uncomment and make this work by providing the necessary operators.
            //uptr2->doSomething();
            //uptr2_cref->doSomething();
        }
    }

    std::cout << "test case 2: deleting in outer context" << std::endl;
    {
        PtrType uptr3(nullptr);
        {
            PtrType uptr4(new TestClass());
            // TODO forbid the assignment operator and think about how this case can still be made work
            uptr3 = uptr4;
        }
    }
*/

    return 0;
}
