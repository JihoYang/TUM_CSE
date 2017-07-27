/*-----------------------------------------------------------------------------*/
/*                                                                             */
/*  4th order Backward Difference Formula solver for State Space systems       */
/*                                                                             */
/*  Weighted Jacobi method is used for solving the resulting SLE               */
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
/*  Final update date: 31/07/2016                                              */
/*                                                                             */
/*-----------------------------------------------------------------------------*/ 

/*
  The solver generates a SLE (system of linear equations) for computing velocity (or x_2) based on the input matrices A, B and vector x
  To avoid confusion, the SLE for velocity (in form of Ax_1 = b) is denoted as A_vel_1 * x = b_vel_1
  
  Please note that relaxation methods are implemented for solving this SLE, hence the matrix A_vel_1 must be diagonally dominant and/or positive definite
  Weighted Jacobi method is used for both simplicity and performance, but better methods can be implemented with ease (e.g. Gauss-Seidel/Successive over-relaxation)

  Over-relaxation factor (omega) of 1.3 is used
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "bdf.h"
#include "helper.h"

void bdf(
         double **x,
         double *u,
         double **A,
         double **B,
         double t_end,
         double dt,
         int itermax,
         int num_body
        ){

    /*CPU time variables*/
    clock_t start_t, end_t, total_t;

    /*Loop variables*/
    int i, j;
    int n;
    double t;   

    /*Velocity SLE variables*/
    double **A_vel_1, **A_vel_2, **A_vel_3, **A_vel_4;
    double *b_vel_1, *b_vel_2, *b_vel_3, *b_vel_4;
    double *v_1_new, *v_2_new, *v_3_new, *v_4_new;
    double *v_1_old;
    double dt_1, dt_2, dt_3, dt_4;
    int p, k;

    /*Matrix-vector multiplications*/
    double Ax_1, Ax_2, Ax_3, Ax_4;
    double Bu_1, Bu_2, Bu_3, Bu_4;
    double A_velx_1;
    double A_velx_2;
    //A_velx_3, A_velx_4;

    /*Jacobi method iteration condition variables*/
    int it;
    double rloc_1;
    double rloc_2; 
    //rloc_3, rloc_4;
    double res, eps = 0.0000001;
    double omega = 1.0;

    /*Start CPU time measurement*/
    start_t = clock();

    /*Define step sizes for each BDF intermediate steps*/
    dt_1 = dt;
    dt_2 = 2 * dt / 3;
    dt_3 = 6 * dt / 11;
    dt_4 = 12 * dt / 25;

    /*Allocate memory for v(n+1)*/
    v_1_new = malloc((num_body - 1) * sizeof(double));
    v_2_new = malloc((num_body - 1) * sizeof(double)); 
    v_3_new = malloc((num_body - 1) * sizeof(double)); 
    v_4_new = malloc((num_body - 1) * sizeof(double)); 
    v_1_old = malloc((num_body - 1) * sizeof(double));

    /*Allocate memory for velocity SLE*/
    A_vel_1   = matrix(0, num_body - 1, 0, num_body - 1);
    A_vel_2   = matrix(0, num_body - 1, 0, num_body - 1);
    A_vel_3   = matrix(0, num_body - 1, 0, num_body - 1);
    A_vel_4   = matrix(0, num_body - 1, 0, num_body - 1);
    b_vel_1   = malloc((num_body - 1) * sizeof(double));
    b_vel_2   = malloc((num_body - 1) * sizeof(double));
    b_vel_3   = malloc((num_body - 1) * sizeof(double));
    b_vel_4   = malloc((num_body - 1) * sizeof(double));

    /*Initialise time counters*/
    t = 0;
    n = 0;

    /*Create A_vel matrix*/
    for (i = 0; i <= num_body - 1; i++){

        /*Index converter*/
        p = 2 * i + 1;

        /*Diagonal elements*/
        A_vel_1[i][i] = 1 - A[p][p-1] * dt_1 * dt_1 - A[p][p] * dt_1;    
        A_vel_2[i][i] = 1 - A[p][p-1] * dt_2 * dt_2 - A[p][p] * dt_2;
        A_vel_3[i][i] = 1 - A[p][p-1] * dt_3 * dt_3 - A[p][p] * dt_3;
        A_vel_4[i][i] = 1 - A[p][p-1] * dt_4 * dt_4 - A[p][p] * dt_4;

        for (j = 0; j <= num_body - 1; j++){
            
            if (j != i){

                /*Index converter*/
                k = (j + 1) * 2 - 1;

                /*Other elements (upper and lower triangular matrices)*/            
                A_vel_1[i][j] = -A[p][k-1] * dt_1 * dt_1 - A[p][k] * dt_1;
                A_vel_2[i][j] = -A[p][k-1] * dt_2 * dt_2 - A[p][k] * dt_2;  
                A_vel_3[i][j] = -A[p][k-1] * dt_3 * dt_3 - A[p][k] * dt_3;  
                A_vel_4[i][j] = -A[p][k-1] * dt_4 * dt_4 - A[p][k] * dt_4;  

            }

        }

        /*Initialise initial guess for v(k+1)*/
        v_1_old[i] = 0;
        v_2_new[i] = 0;
        v_3_new[i] = 0;
        v_4_new[i] = 0;

    }

    /*Iterate over time*/
    while (t <= t_end){

        /*Create b_vel*/
        for (i = 0; i <= num_body - 1; i++){

            Ax_1 = 0;
            Ax_2 = 0;
            Ax_3 = 0;
            Ax_4 = 0;
            Bu_1 = 0;
            Bu_2 = 0;
            Bu_3 = 0;
            Bu_4 = 0;

            p = 2 * i + 1;

            for (j = 0; j <= num_body - 1; j++){

                k = (j + 1) * 2 - 1;

                /*Compute A*x(k)*/
                Ax_1 += A[p][k-1] * x[2*j][n];
                Ax_2 += A[p][k-1] * x[2*j][n];
                Ax_3 += A[p][k-1] * x[2*j][n];
                Ax_4 += A[p][k-1] * x[2*j][n];

                for (j = num_body; j <= num_body * 2 - 1; j++){

                    /*Compute B*u*/
                    Bu_1 += B[p][j] * u[j];
                    Bu_2 += B[p][j] * u[j];
                    Bu_3 += B[p][j] * u[j];
                    Bu_4 += B[p][j] * u[j];

                }

            }

            b_vel_1[i] = x[p][n] + dt_1 * (Ax_1 + Bu_1);
            b_vel_2[i] = x[p][n] + dt_2 * (Ax_2 + Bu_2);
            b_vel_3[i] = x[p][n] + dt_3 * (Ax_3 + Bu_3);
            b_vel_4[i] = x[p][n] + dt_4 * (Ax_4 + Bu_4);

        }

        /*Initialise Jacobi iteration parameters*/
        it = 0;
        res = 1e6;

        /*Jacobi iteration (1st BDF)*/
        while (it < itermax && res > eps){
            
            for (i = 0; i <= num_body - 1; i++){
                
                p = i * 2 + 1;
                
                /*Initialise A_vel_1 * x*/
                A_velx_1 = 0;
                
                /*Pre-compute A_vel_1 * x*/
                for (j = 0; j <= num_body - 1; j++){
                    
                    if (j != i){
                        
                        A_velx_1 += A_vel_1[i][j] * v_1_old[j]; 
                        
                    }
                    
                }
                
                v_1_new[i] = (1 - omega) * v_1_old[i] +  omega * (b_vel_1[i] - A_velx_1) / A_vel_1[i][i];
                
                v_1_old[i] = v_1_new[i];

                
            }
            
            /*Compute residual*/
            rloc_1 = 0;
            
            for (i = 0; i <= num_body - 1; i++){
                
                for (j = 0; j <= num_body - 1; j++){
                    
                    rloc_1 += pow(b_vel_1[i] - A_vel_1[i][j] * v_1_new[j], 2);
                    
                }
                
            }
            
            res = (rloc_1) / (num_body);
            res = sqrt(res);

            it++;
            
        }

        /*Initialise Jacobi iteration parameters*/
        it = 0;
        res = 1e6; 
//
        /*Jacobi iteration (2nd BDF)*/
        while (it < itermax && res > eps){
            
            for (i = 0; i <= num_body - 1; i++){
                
                p = i * 2 + 1;
                
                /*Initialise A_vel_2 * x*/
                A_velx_2 = 0;
                
                /*Pre-compute A_vel_2 * x*/
                for (j = 0; j <= num_body - 1; j++){
                    
                    if (j != i){
                        
                        A_velx_2 += A[i][j] * x[p][n];
                        
                    }
                    
                }
                
                v_2_new[i] = (1 - omega) * x[p][n] +  omega * (b_vel_2[i] - A_velx_2) / A_vel_2[i][i];
                
                v_2_new[i] = 4 * v_1_new[i] / 3 - x[p][n] / 3 + v_2_new[i];
                
                x[p][n] = v_2_new[i];
                
            }
            
            /*Compute residual*/
            rloc_2 = 0;
            
            for (i = 0; i <= num_body - 1; i++){
                
                for (j = 0; j <= num_body - 1; j++){
                    
                    rloc_2 += pow(b_vel_2[i] - A_vel_2[i][j] * v_2_new[j], 2);
                    
                }
                
            }
            
            res = (rloc_2) / (num_body);
            res = sqrt(res);
            
            it++;
            
        }

//        /*Initialise Jacobi iteration parameters*/
//        it = 0;
//        res = 1e6; 
//
//        /*Jacobi iteration (3rd BDF)*/
//        while (it < itermax && res > eps){
//            
//            for (i = 0; i <= num_body - 1; i++){
//                
//                p = i * 2 + 1;
//                
//                /*Initialise A_vel_3 * x*/
//                A_velx_3 = 0;
//                
//                /*Pre-compute A_vel_3 * x*/
//                for (j = 0; j <= num_body - 1; j++){
//                    
//                    if (j != i){
//
//                        A_velx_3 += A[i][j] * x[p][n];
//                        
//                    }
//                    
//                }
//                
//                v_3_new[i] = (1 - omega) * x[p][n] +  omega * (b_vel_3[i] - A_velx_3) / A_vel_3[i][i];
//                
//                v_3_new[i] = 18 * v_2_new[i] / 11 - 9 * v_1_new[i] / 11 + 2 * x[p][n] / 11 + v_3_new[i];
//                
//                x[p][n] = v_3_new[i];
//                
//            }
//            
//            /*Compute residual*/
//            rloc_3 = 0;
//            
//            for (i = 0; i <= num_body - 1; i++){
//                
//                for (j = 0; j <= num_body - 1; j++){
//                    
//                    rloc_3 += pow(b_vel_3[i] - A_vel_3[i][j] * v_3_new[j], 2);
//                    
//                }
//                
//            }
//            
//            res = (rloc_3) / (num_body);
//            res = sqrt(res);
//            
//            it++;
//            
//        }
//
//        /*Initialise Jacobi iteration parameters*/
//        it = 0;
//        res = 1e6; 
//
//        /*Jacobi iteration (4th BDF)*/
//        while (it < itermax && res > eps){
//            
//            for (i = 0; i <= num_body - 1; i++){
//                
//                p = i * 2 + 1;
//                
//                /*Initialise A_vel_4 * x*/
//                A_velx_4 = 0;
//                
//                /*Pre-compute A_vel_4 * x*/
//                for (j = 0; j <= num_body - 1; j++){
//                    
//                    if (j != i){
//                        
//                        A_velx_4 += A[i][j] * x[p][n];
//                        
//                    }
//                    
//                }
//                
//                v_4_new[i] = (1 - omega) * x[p][n] +  omega * (b_vel_4[i] - A_velx_4) / A_vel_4[i][i];
//                
//                v_4_new[i] = 48 * v_3_new[i] / 25 - 36 * v_2_new[i] / 25 + 16 * v_1_new[i] / 25 - 3 * x[p][n] / 25 + v_4_new[i];
//                
//                x[p][n] = v_4_new[i];
//                
//            }
//            
//            /*Compute residual*/
//            rloc_4 = 0;
//            
//            for (i = 0; i <= num_body - 1; i++){
//                
//                for (j = 0; j <= num_body - 1; j++){
//
//                    rloc_4 += pow(b_vel_4[i] - A_vel_4[i][j] * v_4_new[j], 2);
//                    
//                }
//                
//            }
//            
//            res = (rloc_4) / (num_body);
//            res = sqrt(res);
//            
//            it++;
//            
//        }

        /*Update x vector*/
        for (i = 0; i <= num_body - 1; i++){

            p = 2 * i + 1;
            x[p][n+1] = v_2_new[i];
       
        }

        /*Compute displacement using direct implicit euler*/
        for (i = 0; i <= num_body - 1; i++){

            p = 2 * i + 1;
            x[p-1][n+1] = x[p-1][n] + dt * x[p][n+1];

        }

        //print results
        //printf("Timestep = %f,  Number of Jacobi Iteation = %d,    Residual = %f,   x[0][%d] = %f,    x[1][%d] = %f\n", t, it, res, n, x[0][n], n, x[1][n]);
        //printf("\n");

        t = t + dt_1;
        n++;
    
    }

    /*Free memory allocation*/
    free_matrix(A_vel_1, 0, num_body - 1, 0, num_body - 1);
    free_matrix(A_vel_2, 0, num_body - 1, 0, num_body - 1);
    free_matrix(A_vel_3, 0, num_body - 1, 0, num_body - 1);
    free_matrix(A_vel_4, 0, num_body - 1, 0, num_body - 1);
    free(b_vel_1);
    free(b_vel_2);
    free(b_vel_3);
    free(b_vel_4); 
    free(v_1_new);
    free(v_2_new);
    free(v_3_new);
    free(v_4_new);
    free(v_1_old);

    end_t = clock();
    total_t = (long double)(end_t - start_t) / CLOCKS_PER_SEC;
 
    printf("\n");
    printf("\n");                                                                                                                                                                
    printf("Total time taken by CPU for 4th order BDF solver: %lu\n", total_t);
    printf("\n");

}

