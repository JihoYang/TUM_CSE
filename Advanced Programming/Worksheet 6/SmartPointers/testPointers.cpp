/**
 * Program for testing the pointer types auto, unique and shared.
 * Compile with:
 * 		g++ -std=c++11 testPointers.cpp -o testPointers
 * 
 **/

#include <iostream>
#include <memory>
using namespace std;

class T{
	public:
		T():value(10){};
		int value;
};

void test_autoPtr();
void test_sharedPtr();
void test_uniquePtr();

int main(){
	
	/**
	 * Use this integer to switch between the different pointer type
	 **/ 
	int test = 0;
	switch(test){
		case(0): test_autoPtr();
				 break;
		case(1):
				 test_uniquePtr();
				 break;
		case(2):
				 test_sharedPtr();
				 break;
		default:
				 cout << "test has the wrong value. Set it to 0,1 or 2" << endl;
	}
	
	return 0;
}

void test_autoPtr(){
	cout << "##Test auto_ptr##" << endl;
	
	auto_ptr<T> p1(new T());
	cout << p1.get() << endl;
	cout << p1->value << endl;
	
	auto_ptr<T> p2(p1);
	cout << p1.get() << endl;
	cout << p2.get() << endl;
	cout << p2->value << endl;
}

void test_uniquePtr(){
	cout << "##Test unique_ptr##" << endl;
	
	unique_ptr<T> p1(new T());
	cout << p1.get() << endl;
	cout << p1->value << endl;
	
	unique_ptr<T> p2(std::move(p1));
	cout << p1.get() << endl;
	cout << p2.get() << endl;
	cout << p2->value << endl;
}

void test_sharedPtr(){
	cout << "##Test shared_ptr##" << endl;
	
	shared_ptr<T> p1(new T());
	cout << p1->value << endl;
	
	shared_ptr<T> p2(p1);
	cout << p1.get() << endl;
	cout << p1->value << endl;
	cout << p2.get() << endl;
	cout << p2->value << endl;
}


