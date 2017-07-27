#include <iostream>
#include <cmath>
#include <cstdlib>

// memcpy
#include <string.h>

// NOTE: this requires C++11 support
#include <chrono>

typedef void (*mm_mult_ptr)(size_t rows, size_t columns, size_t k, double *A, double *B, double *C);
typedef void (*mm_mult_block_ptr)(size_t rows, size_t columns, size_t k, double *A, double *B, double *C, size_t blocksize);

// TODO: ASSIGNMENT 1/2: reference implementation
// base version - row - column - k
// C[rows,columns] = A [rows,k] * B [k,columns]
void mm_mult_rck(size_t rows, size_t columns, size_t k, double *A, double *B, double *C) {
    for (int row=0; row < rows; ++row) {
        for (int col = 0; col < columns; ++col) {
            for (int i=0; i < k; ++i) {
                if (i==0) {
                    C[row*columns+col] = A[row*k + i] * B[i*columns + col];
                } else {
	                C[row*columns+col] += A[row*k + i] * B[i*columns + col];
                }
            }
        }
    }
}

// TODO: ASSIGNMENT 2
// TODO: loop reordering and loop splitting - col - k - row
void mm_mult_ckr(size_t rows, size_t columns, size_t k, double *A, double *B, double *C) {

}

// TODO: ASSIGNMENT 2
// TODO: loop splitting and loop reordering - k - row - col
void mm_mult_krc(size_t rows, size_t columns, size_t k, double *A, double *B, double *C) {

}

// TODO: ASSIGNMENT 2: optional
// TODO: loop reordering - column - row - k
void mm_mult_crk(size_t rows, size_t columns, size_t k, double *A, double *B, double *C) {

}

// TODO: ASSIGNMENT 2: optional
// TODO: loop reordering and loop - k - col - row
void mm_mult_kcr(size_t rows, size_t columns, size_t k, double *A, double *B, double *C) {

}

// TODO: ASSIGNMENT 2: optional
// TODO: loop reordering and loop splitting - row - k - col
void mm_mult_rkc(size_t rows, size_t columns, size_t k, double *A, double *B, double *C) {

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
        const double mflops_per_second = flops_per_second / 1000000.0;

        std::cout << "sample " << sample 
                  << " elapsed time: " << elapsed_seconds.count() << "s"
                  << " MFLOP/s: " << mflops_per_second
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

        std::cout << "== loop order: column - row - k ==" << std::endl;
        testFunction(mm_mult_crk, rows,columns,k, samples,A,B,C);
     
        std::cout << "== loop order: k - row - column ==" << std::endl;
        testFunction(mm_mult_krc, rows,columns,k, samples,A,B,C);
 
        std::cout << "== loop order: k - column - row ==" << std::endl;
        testFunction(mm_mult_kcr, rows,columns,k, samples,A,B,C);

        std::cout << "== loop order: row - k - column ==" << std::endl;
        testFunction(mm_mult_rkc, rows,columns,k, samples,A,B,C);

        std::cout << "== loop order: column - k - row ==" << std::endl;
        testFunction(mm_mult_ckr, rows,columns,k, samples,A,B,C);


    }
 
    /* ----------------------------- */

    delete[] A;
    delete[] B;
    delete[] C;

    return 0;
}
