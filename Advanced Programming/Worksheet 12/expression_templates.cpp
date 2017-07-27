#include <iostream>
#include <cstdlib>

// TODO: Assignment 4: Task 1: Implement add / multiply

// TODO: Assignment 4: Task 2: Implement power

// TODO: Assignment 4: Task 3: Implement PolynomialTerm

// TODO: Assignment 4: Task 4: Implement Polynomial
template <typename Data, unsigned int N>
class Polynomial {
    public:
        static Data evaluate(const double *coeffs, const Data& value) {
            return 0;
        }
};

int main(int argc, char **argv) {
    // create our polynomial object with fixed order
    const size_t order = 10;
    const size_t n_coefficients = order+1;
    Polynomial<double,order> polynomial;

    // setup the coefficients with random values
    double coefficients[n_coefficients];
    for (size_t i=0; i < n_coefficients; ++i) {
        coefficients[i] = static_cast<double>(rand()) / RAND_MAX;
    }

    // evaluate our polynomial
    double value = 0.0;
    std::cout << "enter value: ";
    std::cin >> value;
    double result = polynomial.evaluate(coefficients, value);
    std::cout << "f(" << value << ") = " << coefficients[0] << "*(" << value << ")**0";
    for (size_t i=1; i < n_coefficients; ++i) {
        std::cout << " + " << coefficients[i] << "*(" << value << ")**" << i;
    }
    std::cout << " = " << result << std::endl;
    return  0;
}
