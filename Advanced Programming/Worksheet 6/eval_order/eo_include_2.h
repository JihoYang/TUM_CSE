#ifndef _TEST_INCLUDE_H_
#define _TEST_INCLUDE_H_

class B {
    public:
        int b;
};

class Depends {
    public:
        Depends(Base *base) {
            copy = base->value;
        }

        ~Depends() { }

        int copy;
};
 
#endif
