/*******************************************************
Enumerated Types For Fiala Tire property file
*******************************************************/

#ifndef _TYR369_H
#define _TYR369_H

#define NO_FRC	448
#define	WSPNMX	5.0e-1

//int maxlen_errmsg = 256;

/**/
const int use_mode				= 1;
const int unloaded_radius		= 7;
const int width					= 8;
const int aspect_ratio			= 10;
const int vertical_stiffness	= 15;
const int vertical_damping		= 16;
const int rolling_resistance	= 17;
const int cslip					= 18;
const int calpha				= 19;
const int cgamma				= 20;
const int umin					= 21;
const int umax					= 22;
const int relaxation_length		= 23;

/*Carcuss Shape Array Parmeters*/
const int max_shape		= 10;
const int n_shape		= 30;
const int shape			= 31;

/*Scaling and Drift parameters*/
const int rrscl_ptr	= 291;
const int fyscl_ptr	= 292;
const int mxscl_ptr	= 293;
const int mzscl_ptr	= 294;
const int plyfrc_ptr	= 295;
const int confrc_ptr	= 296;
const int plytrq_ptr	= 297;
const int contrq_ptr	= 298;

#endif