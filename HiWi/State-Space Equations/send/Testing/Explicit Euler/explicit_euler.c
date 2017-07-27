/*-----------------------------------------------------------------------------*/
/*                                                                             */
/*  Explicit Euler solver for State Space systems                              */
/*                                                                             */  
/*  Developer:  Jiho Yang (MEng)                                               */
/*              M.Sc. candidate, Computational Science & Engineering           */
/*              Technische Universitat Munchen                                 */
/*                                                                             */ 
/*  Work conducted as a student job (HiWi)                                     */
/*  under Seungyong Oh (Dipl.-Ing)                                             */
/*  at Lehrstuhl fur Fordertechnik Materialfluss Logistik                      */
/*  Dept of Mechanical Engineering                                             */
/*  Technische Universitat Munchen, Germany                                    */
/*                                                                             */
/*  Final update date: 03/08/2016                                              */
/*                                                                             */
/*-----------------------------------------------------------------------------*/ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "explicit_euler.h"
#include "helper.h"

void explicit_euler(
                    double **x,
                    double *u,
                    double **A,
                    double **B,
                    double t_end,
                    double dt,
                    int itermax,
                    int dof
                   ){

    /*CPU time variables*/
    clock_t start_t, end_t, total_t;

    /*Loop variables*/
    int i, j;
    int n;
    double t;   

    /*Matrix-vector multiplication variables*/
    double Ax_1;
    double Bu;

    /*Start CPU time measurement*/
    start_t = clock();

    /*Initialise time counters*/
    t = 0;
    n = 0;

    /*Iterate over time*/
    while (t <= t_end){

        for (i = 0; i <= dof * 2 - 1; i++){

            /*Initialise Ax and Bu*/
            Ax_1 = 0;
            Bu   = 0;

            for (j = 0; j <= dof * 2 - 1; j++){

                /*Compute A * x term in numerator*/
                Ax_1 += A[i][j] * x[j][n];

            }

            for (j = 0; j <= dof - 1; j++){

                /*Compute source terms*/
                Bu  += B[i][j] * u[j];

            }

            /*Compute x(t + dt)*/
            x[i][n+1] = x[i][n] + dt * (Ax_1 + Bu);

        }
            
        /*Print out results*/
        printf("Timestep = %f,  x[0][%d] = %f,    x[1][%d] = %f\n", t, n, x[0][n], n, x[1][n]);
        printf("\n");

        t = t + dt;
        n++;

    }

    end_t = clock();
    total_t = (long double)(end_t - start_t) / CLOCKS_PER_SEC;
 
    printf("\n");
    printf("\n");                                                                                                                                                                
    printf("Total time taken by CPU for BDF iteration: %lu\n", total_t);
    printf("\n");

}

