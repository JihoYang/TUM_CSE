#include <iostream>
#include <cmath>
#include <cstdlib>

// memcpy
#include <string.h>

// NOTE: this requires C++11 support
#include <chrono>

typedef void (*mm_mult_ptr)(size_t rows, size_t columns, size_t k, double *A, double *B, double *C);
typedef void (*mm_mult_block_ptr)(size_t rows, size_t columns, size_t k, double *A, double *B, double *C, size_t blocksize);

// base version - row - column - k
// C[rows,columns] = A [rows,k] * B [k,columns]
void mm_mult_rck(size_t rows, size_t columns, size_t k, double *A, double *B, double *C) {
    for (int row=0; row < rows; ++row) {
        for (int col = 0; col < columns; ++col) {
            C[row*columns+col] = A[row*k + 0] * B[0*columns + col];
            for (int i=1; i < k; ++i) {
                C[row*columns+col] += A[row*k + i] * B[i*columns + col];
            }
        }
    }
}

// TODO: ASSIGNMENT 1: cache blocking of standard combination - row - column - k
void mm_mult_rck_block(size_t rows, size_t columns, size_t k, double *A, double *B, double *C, size_t blocksize) {
    // reference code without blocking
    for (int row=0; row < rows; ++row) {
        for (int col = 0; col < columns; ++col) {
            C[row*columns+col] = A[row*k + 0] * B[0*columns + col];
            for (int i=1; i < k; ++i) {
                C[row*columns+col] += A[row*k + i] * B[i*columns + col];
            }
        }
    }
}

void testFunction(mm_mult_ptr mm_mult, size_t rows, size_t columns, size_t k, size_t samples, double *A, double *B, double *C) {
    for (int sample=0; sample < samples; ++sample)
    {
        std::chrono::time_point<std::chrono::system_clock> start;
        std::chrono::time_point<std::chrono::system_clock> end;

        start = std::chrono::system_clock::now();

        mm_mult(rows, columns, k, A, B, C);

        end = std::chrono::system_clock::now();

        std::chrono::duration<double> elapsed_seconds = end-start;
        
        const double total_flops = (2*k-1)*rows*columns;
        const double flops_per_second = total_flops / elapsed_seconds.count();
        const double gflops_per_second = flops_per_second / 1000000000.0;

        std::cout << "sample " << sample 
                  << " elapsed time: " << elapsed_seconds.count() << "s"
                  << " GFLOP/s: " << gflops_per_second
                  << std::endl;
    }
}

void testFunctionBlock(mm_mult_block_ptr mm_mult_block, size_t rows, size_t columns, size_t k, size_t samples, size_t blocksize, double *A, double *B, double *C) {
    for (int sample=0; sample < samples; ++sample)
    {
        std::chrono::time_point<std::chrono::system_clock> start;
        std::chrono::time_point<std::chrono::system_clock> end;

        start = std::chrono::system_clock::now();

        mm_mult_block(rows, columns, k, A, B, C, blocksize);

        end = std::chrono::system_clock::now();

        std::chrono::duration<double> elapsed_seconds = end-start;
        
        const double total_flops = (2*k-1)*rows*columns;
        const double flops_per_second = total_flops / elapsed_seconds.count();
        const double gflops_per_second = flops_per_second / 1000000000.0;

        std::cout << "sample " << sample 
                  << " elapsed time: " << elapsed_seconds.count() << "s"
                  << " GFLOP/s: " << gflops_per_second
                  << std::endl;
    }

}

int main(int argc, char **argv) {
    const size_t samples = 4;

    // NOTE: if you run into "out of memory" problems, reduce this number.
    // If you have a lot of memory, feel free to increase this number.
    size_t maximum_number_of_rows_and_columns = 1024;
 
    double *A = new double[maximum_number_of_rows_and_columns*maximum_number_of_rows_and_columns];
    double *B = new double[maximum_number_of_rows_and_columns*maximum_number_of_rows_and_columns];
    double *C = new double[maximum_number_of_rows_and_columns*maximum_number_of_rows_and_columns];

    for (size_t reference_size = 64; reference_size <= maximum_number_of_rows_and_columns; reference_size *=2) {
        const size_t rows = reference_size;
        const size_t columns = reference_size;
        const size_t k = reference_size;
        std::cout << " ================== " << reference_size << " ================ " << std::endl;
 
        std::cout << "== loop order: row - column - k ==" << std::endl;
        testFunction(mm_mult_rck, rows,columns,k, samples,A,B,C);

        for (size_t blocksize = 4; blocksize <= 32; blocksize*=2) {
            std::cout << " ================== " << reference_size << "(" << blocksize << ") ================ " << std::endl;

            const size_t rows = reference_size;
            const size_t columns = reference_size;
            const size_t k = reference_size;

            std::cout << "== blocked loop order: row - column - k ==" << std::endl;
            testFunctionBlock(mm_mult_rck_block,rows,columns,k,samples,blocksize,A,B,C);
        }
 
        std::cout << std::endl << std::endl;
    }

    delete[] A;
    delete[] B;
    delete[] C;

    return 0;
}
