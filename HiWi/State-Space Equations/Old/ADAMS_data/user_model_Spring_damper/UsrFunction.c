#include"main.h"
#include "slv_c_utils.h"

double sign(double werte){

	double val;

	if (werte >= ZERO){
		val = ONE;
	}else{
		val = -ONE;
	}
	return val;
}

void act_line_parse(char* Value, int ValueLen){

} 

void fialafx(double* fx, double* kappa, double* mu, 
			 double* fz, double* c_slip){
/* DESCRIPTION:

   Calculates the longitudinal tire force for the Fiala tire
   model.

 
 ARGUMENT LIST:

    name    type  storage    use           description
    ======  ====  =========  ===  ====================================
    fx      D.S.     1        R   Longitudinal force

    kappa   D.S.     1        R   Longitudinal Slip

    mu      D.S.     1        R   Effective coefficient of friction

    fz      D.S.     1        R   Vertical force (+downward) in SAE
                                  contact patch axis system.

    c_slip  D.S.     1        R   Partial derivative of longitudinal
                                  force with longitudinal slip at zero
                                  longitudinal slip.
*/
	double fx1, fx2, sstar;
	double temp = 0.0e0;

	sstar = abs((*mu) * (*fz) /(2.0e0 * (*c_slip)));
	if (abs(*kappa) <= sstar){
		(*fx) = (*c_slip) * (*kappa);
	}else{
		fx1 = (*mu) * abs(*fz);
		fx2 = pow(((*mu) * (*fz)),2) / (4.0e0 * abs(*kappa) * (*c_slip));

		(*fx) = (fx1 - fx2) * sign(*kappa);
	}

}

void fialafy(double* fy, double* tz, double* alpha, 
			 double* mu, double* fz, double* tirew, 
			 double* c_alpha){
/* DESCRIPTION:

   Calculates the lateral force and aligning torque for the 
   Fiala tire model.

 ARGUMENT LIST:
    name    type  storage    use           description
    ======  ====  =========  ===  ====================================
    fy      D.S.     1        E   Lateral force in SAE contact patch
                                  axis system.

    tz      D.S.     1        E   Aligning torque in SAE contact patch
                                  axis system.

    alpha   D.S.     1        R   Slip Angle in SAE Axis System.

    mu      D.S.     1        R   Effective coefficient of friction

    fz      D.S.     1        R   Vertical force (+downward) in SAE
                                  contact patch axis system.

   tirew    D.S.     1        R   Tire width.

   c_alpha  D.S.     1        R   Partial derivative of lateral
                                  force with lateral slip at zero
                                  lateral slip.
*/
	double astar;
	double h;
	double zerlim = 1.0e-10;

	if( abs(*alpha) <= zerlim){
		*fy = ZERO;
		*tz = ZERO;
	}else{
		astar = atan(abs(3.0e0 * (*mu) * (*fz) / (*c_alpha)));
		if(abs(*alpha) <= astar){
			h = ONE - (*c_alpha) * abs(tan(*alpha)) / (3.0e0 * (*mu) * abs(*fz));
			*fy = -(*mu) * abs(*fz) * (ONE - pow(h,3)) * sign(*alpha);
			*tz = (*mu) * abs(*fz) * (*tirew) * (ONE - h) * pow(h,3) * sign(*alpha);
		}else{
			*fy = -(*mu) * abs(*fz) * sign(*alpha);
			*tz = ZERO;
		}
	}
}

void fialaty(double* ty, double* wy, double* fz, double* rolco){
/*
 DESCRIPTION:

   Calculates the rolling resistance moment (ty) for the
   Fiala tire model.

     s = step( abs(wy), 0.125, 0.0, 0.5, 1.0 )
    ty = -rolco * fz * sign( 1.0, wy) * s

 
 ARGUMENT LIST:

    name    type  storage    use           description
    ======  ====  =========  ===  ====================================
    ty      D.S.     1        E   Rolling resistance moment in SAE 
                                  contact patch axis system.

    wy      D.S.     1        R   Angular velocity about tire spin
                                  axis.

    fz      D.S.     1        R   Vertical force (+downward) in SAE
                                  contact patch axis system.

   rolco    D.S.     1        R   Rolling resistance coefficient
*/

	double wy_min = 0.125e0;
	double wy_max = 0.5e0;
	int iord =0;
	double val;
	int errflg;

	c_step(abs(*wy), wy_min, ZERO, wy_max, ONE, iord, &val, &errflg);

	*ty = -sign(*wy) * val * (*rolco) * (*fz);

}