#ifndef _TEST_INCLUDE_H_
#define _TEST_INCLUDE_H_

class A {
    public:
        int a;
};

class Base {
    public:
        Base() { value = 10; }
        ~Base() { }

        int value;
};

extern Base base;

#endif
