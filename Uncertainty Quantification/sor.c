#include "sor.h"
#include <math.h>

void sor(
  double omg,
  double dx,
  double dy,
  int    imax,
  int    jmax,
  int    k,
  int    mode,
  double ***P,
  double ***RS,
  double *res
) {
  int i,j;
  double rloc;
  double coeff = omg/(2.0*(1.0/(dx*dx)+1.0/(dy*dy)));

  /* SOR iteration */
  for(i = 1; i <= imax; i++) {
    for(j = 1; j<=jmax; j++) {
      P[k][i][j] = (1.0-omg)*P[k][i][j]
              + coeff*(( P[k][i+1][j]+P[k][i-1][j])/(dx*dx) + ( P[k][i][j+1]+P[k][i][j-1])/(dy*dy) - RS[k][i][j]);
    }
  }

  /* compute the residual */
  rloc = 0;
  for(i = 1; i <= imax; i++) {
    for(j = 1; j <= jmax; j++) {
      rloc += ( (P[k][i+1][j]-2.0*P[k][i][j]+P[k][i-1][j])/(dx*dx) + ( P[k][i][j+1]-2.0*P[k][i][j]+P[k][i][j-1])/(dy*dy) - RS[k][i][j])*
              ( (P[k][i+1][j]-2.0*P[k][i][j]+P[k][i-1][j])/(dx*dx) + ( P[k][i][j+1]-2.0*P[k][i][j]+P[k][i][j-1])/(dy*dy) - RS[k][i][j]);
    }
  }

  rloc = rloc/(imax*jmax);
  rloc = sqrt(rloc);
    
  /* set residual */
  *res = rloc;


  /* set boundary values */
  for(i = 1; i <= imax; i++) {
    
    P[k][i][0] = P[k][i][1];
    P[k][i][jmax+1] = P[k][i][jmax];
    
  }

  for(j = 1; j <= jmax; j++) {
    
    P[k][0][j] = P[k][1][j];
    P[k][imax+1][j] = P[k][imax][j];
    
  }
}

