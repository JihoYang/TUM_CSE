#include "helper.h"
#include "visual.h"
#include <stdio.h>

/* Column convention: 1. Time   2. Displacement (or x_1)  3. Velocity (or x_2) ...*/
/* Line <23> needs to be adjusted accordingly for different number of bodies*/
void export_result(
                   const char *method,
                   double t_end,
                   double dt,
                   double **x
                  ){
     
    int  t_count;                                                                                                      
    char FileName[80];
 
    sprintf(FileName, "%s_displacement.csv", method);
         
    FILE *f = fopen(FileName, "wb");

    for (t_count = 0; t_count <= t_end / dt; t_count++){

        fprintf(f, " %f %f %f\n", dt * t_count,  x[0][t_count], x[1][t_count]);

    }
 
    fclose(f);
 
}

