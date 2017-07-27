#include <iostream>
int main(int argc, char *argv[]) {
    int xvalues[10];
    for (char i=0; i < 128; i++) {
        std::cin >> xvalues[static_cast<int>(i)];
    }
    double a = static_cast<double>(1/3);
    double b = 5;
    double c = 10;

    for (int i=0; i < 128; i++) {
        int x = xvalues[i];
        double result = a * ((x+1)^2) + b * x + c;
        std::cout << "x=" << xvalues[i]
                  << ",result=" << result 
                  << std::endl;
    }
    return 0;
}
