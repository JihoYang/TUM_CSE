#ifndef __HELPER_H__                                                                                                                      
#define __HELPER_H__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double **matrix( 
                int nrl, 
                int nrh, 
                int ncl, 
                int nch 
               );

void free_matrix( 
                 double **m, 
                 int nrl, 
                 int nrh, 
                 int ncl, 
                 int nch 
                );

void init_matrix( 
                 double **m, 
                 int nrl, 
                 int nrh, 
                 int ncl, 
                 int nch, 
                 double a
                );

#endif

