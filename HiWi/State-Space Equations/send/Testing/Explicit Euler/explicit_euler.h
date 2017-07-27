#ifndef __BDF_H_                                                                                                                    
#define __BDF_H_

void explicit_euler(
                    double **x,
                    double *u,
                    double **A,
                    double **B,
                    double t_end,
                    double dt,
                    int itermax,
                    int dof
                   );

#endif

