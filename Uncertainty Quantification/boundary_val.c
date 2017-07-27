#include "boundary_val.h"

//Set the Boundary Conditions (velocities)//
void boundaryvalues(
                    int imax,
                    int jmax,
                    int k,
                    int mode,
                    double ***U,
                    double ***V
                    ){
    
    int i, j;
    double U_mv_wall = 1.0;
    
    //BC @ left wall//
    for(j = 0; j <= jmax; j++){
        
        U[k][0][j] = 0;
        V[k][0][j] = -V[k][1][j];
        
    }
    
    //BC @ right wall//
    for(j = 0; j <= jmax; j++) {
        
        U[k][imax][j] = 0;
        V[k][imax+1][j] = -V[k][imax][j];
        
    }
    
    //BC @ floor//
    for(i = 0; i <= imax; i++) {
        
        U[k][i][0] = -U[k][i][1];
        V[k][i][0] = 0;
        
    }
    
    //BC @ moving wall//
    for(i = 0; i <= imax; i++) {
        
        U[k][i][jmax+1] = 2 * U_mv_wall - U[k][i][jmax];
        V[k][i][jmax] = 0;

        if (k != 0){

        U[k][i][jmax+1] =  - U[k][i][jmax];
        V[k][i][jmax] = 0;
    
        }
        
    }
    
}

