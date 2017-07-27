/*------------------------------------------------------------------------*/
/*                                                                        */
/*  HELPER.C provides some useful auxiliary functions                     */
/*                                                                        */
/*  Developer: Jiho Yang (MEng)                                           */
/*             M.Sc. candidate, Computational Science & Engineering       */
/*             Technische Universitat Munchen                             */
/*                                                                        */
/*  Final update: 21/07/2016                                              */
/*                                                                        */
/*------------------------------------------------------------------------*/

#include <stdio.h>
#include "helper.h"

/*Generate Matrix*/
double **matrix(int nrl, int nrh, int ncl, int nch){
 
     int i;
     int nrow = nrh - nrl + 1;    /* compute number of lines */
     int ncol = nch - ncl + 1;    /* compute number of columns */
          
     double **m = (double **) malloc((size_t)( nrow * sizeof(double*)));
                  
     for (i = 0; i < nrow; i++){
          
        m[i] = (double *) malloc((size_t)( ncol * sizeof(double)));
     
     } 
     
     return m;
                                                                                                                                                                     
}

 
/*Free allocated memories from a matrix*/
void free_matrix(
                 double **m, 
                 int nrl, 
                 int nrh, 
                 int ncl, 
                 int nch
                ){  

     double **m2d = m;
     
     free( m2d );

}


/*Initialise matrix*/
void init_matrix( 
                 double **m, 
                 int nrl, 
                 int nrh, 
                 int ncl, 
                 int nch, 
                 double a
                ){

    int i,j;

    for(i = nrl; i <= nrh; i++){

        for(j = ncl; j <= nch; j++){
    
            m[i][j] = a;

        }

    }

}

