#include <iostream>
#include <cstdlib>

const int N = 100;

int main(int argc, char **argv) {
    double a[N];
    double b[N];
    double c[N];
    
    for (size_t i=0; i < N; ++i) {
        a[i] = static_cast<double>(rand()) / RAND_MAX;
        b[i] = static_cast<double>(rand()) / RAND_MAX;
        c[i] = 0.0;
    }

    int n=0;
    while (n < N) {
        c[n] = a[n] + b[n];
        n++;
    };

    if (argc == -1) {
        for (size_t i=0; i < N; ++i) {
           std::cout << c[i] << std::endl;
        }
    }
    return 0;
}
