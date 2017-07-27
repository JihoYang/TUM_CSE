/*-----------------------------------------------------------------------------*/
/*                                                                             */
/*  Implicit Euler solver for solving State Space systems                      */
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
/*  Final update date: 04/08/2016                                              */
/*                                                                             */
/*-----------------------------------------------------------------------------*/ 

/*
  The solver generates a SLE (system of linear equations) for computing velocity (or x_2) based on the input matrices A, B and vector x
  To avoid confusion, the SLE for velocity (in form of Ax = b) is denoted as A_vel * v_1_new = b_vel
  
  Please note that relaxation methods are implemented for solving this SLE, hence the matrix A_vel must be diagonally dominant and/or positive definite
  Weighted Jacobi method is used for both simplicity and performance, but better methods can be implemented with ease (e.g. Gauss-Seidel/Successive over-relaxation)

  Over-relaxation factor (omega) of 1.3 is used
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "implicit_euler.h"
#include "helper.h"

void implicit_euler(
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
    int counter;
  
    /*Buffers*/
    double **x_buf;
    double **A_buf;
    double **B_buf;

    /*Velocity SLE variables*/
    double **A_vel;
    double *b_vel;
    double *v_1_new;
    double *v_1_old;
    int p, k;

    /*Matrix-vector multiplications*/
    double Ax;
    double Bu;
    double A_velx;

    /*Jacobi method iteration condition variables*/
    int it;
    double rloc, res, eps = 0.0000001;
    double omega = 1.3;

    /*Start CPU time measurement*/
    start_t = clock();

    /*Allocate memory for v(n+1)*/
    v_1_new = malloc((dof - 1) * sizeof(double));
    v_1_old = malloc((dof - 1) * sizeof(double));

    /*Allocate memory for buffers*/
    x_buf   = matrix(0, dof * 2 - 1, 0, t_end / dt); 
    A_buf   = matrix(0, dof * 2 - 1, 0, dof * 2 - 1);
    B_buf   = matrix(0, dof * 2 - 1, 0, dof - 1);

    /*Allocate memory for velocity SLE*/
    A_vel   = matrix(0, dof - 1, 0, dof - 1);
    b_vel   = malloc((dof - 1) * sizeof(double));

    /*Initialise time counters*/
    t = 0;
    n = 0;

    /*Create buffer space for x and u vectors*/
    for (i = 0; i <= dof * 2 - 1; i++){

        /*x vector*/
        for (j = 0; j <= t_end / dt; j++){

            x_buf[i][j] = x[i][j];
            
        }

    }

    /*Create buffer space for A and B matrices*/
    for (i = 0; i <= dof * 2 - 1; i++){

        for (j = 0; j <= dof * 2 - 1; j++){

            A_buf[i][j] = A[i][j];

        }

        for (j = 0; j <= dof - 1; j++){

            B_buf[i][j] = B[i][j];

        }

    }

    /*Adjust inputs (displacement)*/
    counter = 0;

    for (i = 0; i <= dof - 1; i++){

        p = 2 * i;

        /*x vector*/
        for (j = 0; j <= t_end/ dt; j++){

            x[p][j] = x_buf[counter][j];

        }

        /*A and B matrices*/
        for (j = 0; j <= dof * 2 - 1; j++){

            A[p][j] = A_buf[counter][j];

        }

        for (j = 0; j <= dof - 1; j++){

            B[p][j] = B_buf[counter][j];

        }

        counter ++;

    }

    counter = 0;

    for (j = 0; j <= dof - 1; j++){

        p = 2 * j;

        /*A matrix (columns)*/
        for (i = 0; i <= dof * 2 - 1; i++){

            A[i][p] = A_buf[i][counter];

        }

        counter++;

    }

    /*Adjust inputs (velocity)*/
    counter = dof;

    for (i = 0; i <= dof - 1; i++){

        p = 2 * i + 1;

        /*x vector*/
        for (j = 0; j <= t_end / dt; j++){

            x[p][j] = x_buf[counter][j];

        }

        /*A and B matrix*/
        for (j = 0; j <= dof * 2 - 1; j++){

            A[p][j] = A_buf[counter][j];

        }

        for (j = 0; j <= dof - 1; j++){

            B[p][j] = B_buf[counter][j];

        }

        counter ++;

    }

    counter = dof;

    for (j = 0; j <= dof - 1; j++){

        p = 2 * j + 1;

        /*A matrix (columns)*/
        for (i = 0; i <= dof * 2 - 1; i++){

            A[i][p] = A_buf[i][counter];

        }

        counter ++;

    }

    /*Create A_vel matrix*/
    for (i = 0; i <= dof - 1; i++){

        /*Index converter*/
        p = 2 * i + 1;

        /*Diagonal elements*/
        A_vel[i][i] = 1 - A[p][p-1] * dt * dt - A[p][p] * dt;    

        for (j = 0; j <= dof - 1; j++){
            
            if (j != i){

                /*Index converter*/
                k = (j + 1) * 2 - 1;

                /*Other elements (upper and lower triangular matrices)*/            
                A_vel[i][j] = -A[p][k-1] * dt * dt - A[p][k] * dt;

            }

        }

        /*Initialise initial guess for v(k+1)*/
        v_1_old[i] = 0;

    }

    /*Iterate over time*/
    while (t <= t_end){

        /*Create b_vel*/
        for (i = 0; i <= dof - 1; i++){

            Ax = 0;
            Bu = 0;

            p = 2 * i + 1;

            for (j = 0; j <= dof - 1; j++){

                k = (j + 1) * 2 - 1;

                /*Compute A*x(k)*/
                Ax += A[p][k-1] * x[2*j][n];

                /*Compute B*u*/
                Bu += B[p][j] * u[j];

            }

            b_vel[i] = x[p][n] + dt * (Ax + Bu);

        }

        /*Initialise Jacobi iteration parameters*/
        it = 0;
        res = 1e6;

        /*Jacobi iteration*/
        while (it < itermax && res > eps){

            for (i = 0; i <= dof - 1; i++){

                p = i * 2 + 1;

                /*Initialise A_vel * x*/
                A_velx = 0;
       
                /*Pre-compute A_vel * x*/
                for (j = 0; j <= dof - 1; j++){

                    if (j != i){

                        A_velx += A_vel[i][j] * v_1_old[j];

                    }

                }

                v_1_new[i] = (1 - omega) * v_1_old[i] +  omega * (b_vel[i] - A_velx) / A_vel[i][i];
                v_1_old[i] = v_1_new[i];

            }

            /*Compute residual*/
            rloc = 0;

            for (i = 0; i <= dof - 1; i++){

                for (j = 0; j <= dof - 1; j++){

                    rloc += pow(b_vel[i] - A_vel[i][j] * v_1_new[j], 2);
    
                }

            }

            res = rloc / dof;
            res = sqrt(res);

            it++;

        }

        /*Update x vector*/
        for (i = 0; i <= dof - 1; i++){

            p = 2 * i + 1;
            x[p][n+1] = v_1_new[i];
       
        }

        /*Compute displacement using direct implicit Euler*/
        for (i = 0; i <= dof - 1; i++){

            p = 2 * i + 1;
            x[p-1][n+1] = x[p-1][n] + dt * x[p][n+1];

        }

        /*print results*/
        /*printf("Timestep = %f,  Number of Jacobi Iteation = %d,    Residual = %f,   x[0][%d] = %f,    x[1][%d] = %f\n", t, it, res, n, x[0][n], n, x[1][n]);*/   
        /*printf("\n");*/

        t = t + dt;
        n++;
    
    }

    /*Update buffer space for inputs*/
    for (i = 0; i <= dof * 2 - 1; i++){

        /*x vector*/
        for (j = 0; j <= t_end / dt; j++){

            x_buf[i][j] = x[i][j];
            
        }

        /*A and B matrices*/
        for (j = 0; j <= dof * 2 - 1; j++){

            A_buf[i][j] = A[i][j];

        }

        for (j = 0; j <= dof - 1; j++){

            B_buf[i][j] = B[i][j];

        }

    }

    /*Adjust back the inputs (displacement)*/
    counter = 0;

    for (i = 0; i <= dof - 1; i++){

        p = 2 * i;

        /*x vector*/
        for (j = 0; j <= t_end/ dt; j++){

            x[counter][j] = x_buf[p][j];

        }

        /*A and B matrices*/
        for (j = 0; j <= dof * 2 - 1; j++){

            A[counter][j] = A_buf[p][j];

        }

        for (j = 0; j <= dof - 1; j++){

            B[counter][j] = B_buf[p][j];

        }

        counter ++;

    }
    
    counter = 0;
    
    for (j = 0; j <= dof - 1; j++){
        
        p = 2 * j;
        
        /*A matrix (columns)*/
        for (i = 0; i <= dof * 2 - 1; i++){
            
            A[i][counter] = A_buf[i][p];
            
        }
        
        counter++;
        
    }

    /*Adjust back the inputs (velocity)*/
    counter = dof;

    for (i = 0; i <= dof - 1; i++){

        p = 2 * i + 1;

        /*x vector*/
        for (j = 0; j <= t_end / dt; j++){

            x[counter][j] = x_buf[p][j];

        }

        /*A and B matrix*/
        for (j = 0; j <= dof * 2 - 1; j++){

            A[counter][j] = A_buf[p][j];

        }

        for (j = 0; j <= dof - 1; j++){

            B[counter][j] = B_buf[p][j];

        }

        counter ++;

    }
    
    counter = dof;
    
    for (j = 0; j <= dof - 1; j++){
        
        p = 2 * j + 1;
        
        /*A matrix (columns)*/
        for (i = 0; i <= dof * 2 - 1; i++){
            
            A[i][counter] = A_buf[i][p];
            
        }
        
        counter ++;
        
    }
    
    /*Free memory allocation*/
    free_matrix(A_vel, 0, dof - 1, 0, dof - 1);
    free(b_vel);
    free(v_1_new);
    free(v_1_old);

    /*Free buffer memory*/
    free_matrix(x_buf, 0, dof * 2 - 1, 0, t_end / dt);
    free_matrix(A_buf, 0, dof * 2 - 1, 0, dof * 2 - 1);
    free_matrix(B_buf, 0, dof * 2 - 1, 0, dof - 1); 

    end_t = clock();
    total_t = (long double)(end_t - start_t) / CLOCKS_PER_SEC;
 
    printf("\n");
    printf("\n");                                                                                                                                                                
    printf("Total time taken by CPU for implicit euler solver: %lu seconds\n", total_t);
    printf("\n");

}

