#ifndef _AC_TIR_JOBFLG_H
#define	_AC_TIR_JOBFLG_H
#endif

/******************************************************************
	Holds the jobflg types for the Standard Tire Interface (STI):

     -1 = No evaluation
      0 = Normal Mode
      1 = Inquire:  Subroutine must return the actually used dimensions of
                    NTYPAR, NDEQVR, NVARS, NIWORK, NWORK
      2 = First Initialization (new model)
      3 = Re-Initialization (mid-simulation)
      4 = Successful Simulaton Step (handled from sen980)
      5 = Differencing mode
      6 = End Simulation (handled from acarSDI_fin)

 jobflg = 0   Normal Mode
 jobflg = 1   Subroutine must return the actually used dimensions of
 jobflg = 2   First initialization before starting the simulation
 jobflg = 3   Re-initialization during simulation
 jobflg = 4   Successful step (not used in A/Car)
 jobflg = 5   DFLAG is .TRUE. (A/Car deviation)
 jobflg = 99  Final call of tire model (not used in A/Car)

*****************************************************************/
	const int noeval	= -1;
	const int normal	= 0;
	const int inquire	= 1;    
	const int init		= 2;
	const int reset		= 3;
	const int sstep		= 4;
	const int diff		= 5;
	const int endsim	= 99;