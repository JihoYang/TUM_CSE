// implementation according to worksheet requirements
template <typename DataType> class TArray; 
class Base { public: int a; };
class Derived : public Base { public: int b; };
int main() {
    Base b;
    Derived d;
    b = d; // is working

    TArray<Base> base_array(10);
    TArray<Derived> derived_array(10);
    base_array = derived_array; // is not working
    return 0;
}
