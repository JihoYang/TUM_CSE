{\rtf1\ansi\ansicpg1252\cocoartf1404\cocoasubrtf340
{\fonttbl\f0\fmodern\fcharset0 Courier;}
{\colortbl;\red255\green255\blue255;}
\paperw11900\paperh16840\margl1440\margr1440\vieww10800\viewh8400\viewkind0
\deftab720
\pard\pardeftab720\partightenfactor0

\f0\fs26 \cf0 \expnd0\expndtw0\kerning0
#include <iostream>\
#include <iomanip>\
#include <typeinfo>\
\
using namespace std;\
\
\
/**\
 * Template class implementing unique pointer semantics for arbitrary types\
 */\
template<class T>\
class UniquePointer \{\
public:\
\
    // TODO complete the implementation\
    // HINT: do not lose information.\
    /// create a unique pointer to hold the instance passed\
    UniquePointer(T* obj) : _obj(obj), _index(Index++) \{\};\
\
    // TODO complete the implementation\
    /// "move" the instance into a new instance of the unique pointer type\
    UniquePointer(UniquePointer& uptr) : _obj(uptr._obj), _index(Index++) \{\
        uptr._obj = nullptr;\
    \};\
\
	/// forbid the default copy constructor (and any copy constructor from a const reference)\
    UniquePointer(const UniquePointer& uptr) = delete;\
\
    // TODO complete the implementation\
    /// destructor that cleans up if necessary, also notifies about actions taken\
    ~UniquePointer() \{\
        std::cout << "Deleting index " << _index << " of class " << typeid(_obj).name();\
        if (_obj != nullptr) \{\
            std::cout << ", deleting the object" << std::endl;\
            delete _obj;\
        \} else \{\
            std::cout << ", object already destroyed" << std::endl;\
        \}\
    \};\
\
    // TODO implement _clean_ assignment, and later forbid the operator explicitly (C++11)\
    // assignment operator for unique pointer type -- is it needed? is it clean?\
    void\
    operator=(UniquePointer& uptr) = delete;\
/*    \{\
        if (_obj != nullptr) \{\
            std::cout << "Index " << _index << " of class " << typeid(_obj).name()\
                    << ", deleting the object upon new assignment" << std::endl;\
        \}\
        _obj      = uptr._obj;\
        uptr._obj = nullptr;\
    \};\
*/\
\
    // TODO implement function(s) to give convenient access to the instance\
\
    /// dereference to work with the instance\
    T*\
    operator->() \{ return _obj; \};\
\
    /// TODO: const dereference to work with instance the in const-contexts\
    const T*\
    operator->() const \{ return _obj; \};\
\
private:\
\
    /// the instance all this fuss is about\
    T*     _obj;\
 \
    /////////////////////////////////////////////////////////////// \
    /// the variables below are useful for debugging and \
    /// for better understanding.\
   \
    /// an index to keep track of the destructor calls\
    size_t _index;\
\
    /// counter to identify all instances of this unique pointer type\
    static size_t Index;\
\};\
\
// initilizes the class-wide index counter\
template<class T>\
size_t UniquePointer<T>::Index = 0lu;\
\
\
/**\
 * Simple test class that does not really matter\
 */\
struct TestClass\
\{\
    TestClass() \{\
        std::cout << "Creating instance of TestClass" << std::endl;\
    \};\
\
    ~TestClass() \{\
        std::cout << "Deleting instance of TestClass" << std::endl;\
    \};\
\
    void\
    doSomething() const \{\
        std::cout << "Doing something" << std::endl;\
    \};\
\};\
\
int main() \{\
\
    typedef UniquePointer<TestClass> PtrType;\
\
    std::cout << "test case 1: deleting in inner context" << std::endl;\
    \{\
        PtrType uptr1(new TestClass());\
        \{\
            PtrType uptr2(uptr1);\
            const PtrType& uptr2_cref = uptr2;\
\
            //TODO uncomment and make this work by providing the necessary operators.\
            uptr2->doSomething();\
\
// lambda version\
//            auto doSomethingFunc = [] (const PtrType& uptr) -> void \{\
//                uptr->doSomething();\
//            \};\
//            doSomethingFunc(uptr2);\
\
// non lambda version\
//            uptr2_cref->doSomething();\
        \}\
    \}\
\
    std::cout << "test case 2: deleting in outer context" << std::endl;\
    \{\
        PtrType uptr3(nullptr);\
        \{\
            PtrType uptr4(new TestClass());\
            // TODO forbid the assignment operator and think about how this case can still be made work\
            //uptr3 = uptr4;\
        \}\
    \}\
\
    return 0;\
\}\
\
}