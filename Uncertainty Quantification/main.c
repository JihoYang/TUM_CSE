//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  2D Stochastic Galerkin Finite Difference CFD Solver for Cavity Problem  //
//                                                                          //
//  Developer: Jiho Yang (MEng)                                             //
//             M.Sc. candidate, Computational Scinece & Engineering         //
//             Technical University of Munich, Germany                      //
//                                                                          //
//  Work conducted for 2016 Seminar Uncertainty Quantification              //
//                                                                          //
//  Final update date: 23/05/2016                                           //
//                                                                          //
//  ISSUES: Diverges for Re > 1000                                          // 
//                                                                          //
////////////////////////////////////////////////////////////////////////////// 

#include <stdio.h>
#include <time.h>
#include "helper.h"
#include "init.h"
#include "boundary_val.h"
#include "uvp.h"
#include "sor.h"
#include "visual.h"


int main(int argn, char** args){
    
    clock_t start_t, end_t, total_t;

    //***U,V,P must be defined here//
    double ***U, ***V, ***P, ***F, **F_sg1, **F_sg2, ***G, **G_sg1, **G_sg2, ***RS, nu_0, nu_1, lambda;
    const char *szFileName = "cavity100.dat";
    double Re, UI, VI, PI, GX, GY, t_end, xlength, ylength, dt, dx, dy, alpha, omg, tau, eps, dt_value;
    double res = 0, t = 0, n = 0;
    int imax, jmax, itermax, it, k, mode;
   
    start_t = clock();
    
    /*Read the program configuration file using read_parameters()*/
    read_parameters(szFileName, &Re, &UI, &VI, &PI, &GX, &GY, &t_end, &xlength, &ylength, &dt, &dx, &dy, &imax, &jmax, &alpha, &omg, &tau, &itermax, &eps, &dt_value, &mode);

    printf("\n");
    printf("Simulation Start\n");
    printf("\n");

    //Define viscosity values (nu_0: average viscosity, nu_1: standard deviation of viscosity)
    nu_0 = 1/Re;
    nu_1 = nu_0 *0.3;
    lambda = nu_1 / nu_0;    

    /*Set up the matrices (arrays) needed using the matrix() command*/
    U     = matrix(0, imax  , 0, jmax+1, 0, mode);
    V     = matrix(0, imax+1, 0, jmax  , 0, mode);
    P     = matrix(0, imax+1, 0, jmax+1, 0, mode);
    F     = matrix(0, imax  , 0, jmax+1, 0, mode);
    G     = matrix(0, imax+1, 0, jmax  , 0, mode);
    RS    = matrix(0, imax+1, 0, jmax+1, 0, mode);

    F_sg1 = matrix_2d(0, imax  , 0, jmax+1);
    F_sg2 = matrix_2d(0, imax  , 0, jmax+1);
    G_sg1 = matrix_2d(0, imax+1, 0, jmax);
    G_sg2 = matrix_2d(0, imax+1, 0, jmax);

    /*Set the initial conditions by using init_uvp from init.c*/
    init_uvp(UI, VI, PI, imax, jmax, mode, U, V, P);
                
    /*Iterate over time*/ 
    while(t <= t_end){
            
        /*Compute the maximal time step size*/
        calculate_dt(k, Re, tau, &dt, dx, dy, imax, jmax, U, V);
       
        //Iterate over each polynomial 
        for(k = 0; k <= mode; k++){
    
            /*Set the boundary conditions for each time step*/
            boundaryvalues(imax, jmax, k, mode, U, V);

            /*Solve F(n) and G(n)*/
            calculate_fg(lambda, Re, GX, GY, alpha, dt, dx, dy, imax, jmax, k, mode, U, V, F, F_sg1, F_sg2, G, G_sg1, G_sg2);
        
            /*Solve the RHS of the Pressure Poisson Equation*/
            calculate_rs(dt, dx, dy, imax, jmax, k, mode, F, G, RS);
        
            /*Solve the whole Pressure Poisson Equation and compute P(n+1)*/
            it = 0;
            res = 1e6;
        
            while(it < itermax && res > eps){
                
                sor(omg, dx, dy, imax, jmax, k, mode, P, RS, &res);
                
                it++;
                
            }

            /*Update the velocities*/
            calculate_uv(dt, dx, dy, imax, jmax, k, mode, U, V, F, G, P);

            //Free the F_sg1,2 and G_sg1,2 matrices
            free_matrix_2d(F_sg1, 0, imax ,  0, jmax+1);
            free_matrix_2d(F_sg2, 0, imax ,  0, jmax+1);
            free_matrix_2d(G_sg1, 0, imax+1, 0, jmax);
            free_matrix_2d(G_sg2, 0, imax+1, 0, jmax);

            //Re-allocate memory for F_sg1,2 and G_sg1,2 matrices
            F_sg1 = matrix_2d(0, imax  , 0, jmax+1);
            F_sg2 = matrix_2d(0, imax  , 0, jmax+1);
            G_sg1 = matrix_2d(0, imax+1, 0, jmax);
            G_sg2 = matrix_2d(0, imax+1, 0, jmax);

            printf("Polynomial order k = %i finished for t = %f     |     Residual = %f\n", k, t, res);

       }

       t = t + dt;
       n++;
        
   write_vtkFile("cavity", n, mode, xlength, ylength, imax, jmax, dx, dy, U, V, P);

   }

   //Export the solutions for each k//

   //Reset the matrices// 
   free_matrix(U , 0, imax  , 0, jmax+1);
   free_matrix(V , 0, imax+1, 0, jmax  );
   free_matrix(P , 0, imax+1, 0, jmax+1);
   free_matrix(F , 0, imax  , 0, jmax+1);
   free_matrix(G , 0, imax+1, 0, jmax  );
   free_matrix(RS, 0, imax+1, 0, jmax+1);
  
   free_matrix_2d(F_sg1, 0, imax  , 0, jmax+1);
   free_matrix_2d(F_sg2, 0, imax  , 0, jmax+1);
   free_matrix_2d(G_sg1, 0, imax+1, 0, jmax  );
   free_matrix_2d(G_sg2, 0, imax+1, 0, jmax  );
        
   end_t = clock();
    
   total_t = (long double)(end_t - start_t) / CLOCKS_PER_SEC;

   printf("Total time taken by CPU: %lu\n", total_t  );
   printf("Exiting of the program...\n");
    
   return -1;
    
}

