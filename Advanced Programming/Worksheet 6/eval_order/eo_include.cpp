#include <iostream>
#include <cstdlib>

#include "eo_include_1.h"
#include "eo_include_2.h"

void header_guards() {
    A a;
    B b; // TODO: B is not known by the compiler, 
         // though the header file is included? 
         // Why? Fix it!

    std::cout << "we do not want a to be optimized away: " << a.a << std::endl;
    std::cout << "we do not want b to be optimized away: " << b.b << std::endl;
}
 
// TODO: why do we not get the expected value?
// Where do we have to put this declaration instead?
static Depends depends(&base);

int main(int argc, char **argv) {

    header_guards();
    
    std::cout << "expected value: " << base.value 
              << ", value of depends: " << depends.copy << std::endl;
    return 0;
}
