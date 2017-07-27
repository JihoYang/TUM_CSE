#include"main.h"
#include"tyr369.h"
#include"slv_c_utils.h"
#include"ac_tir_jobflg.h"

extern void rpf369(int* NCHTDS, char* CHTDST, int* IDTYRE, int* NTYPAR, double* TYPARR);

void TYR369(int*	NDEV,	int*	ISWITCH,	int*	JOBFLG, 
			int*	IDTYRE,	double*	TIME,		double*	DIS, 
			double* TRAMAT,	double*	ANGTWC,		double*	VEL, 
			double*	OMEGA,	double*	OMEGAR,		int*	NDEQVR, 
			double*	DEQVAR, int*	NTYPAR,		double*	TYPARR, 
			int*	NCHTDS, char*	CHTDST,		int*	LenCHTDST, 
			void (*ROAD)(int*	JOBFLG, int*	IDTYRE,		double* TIME, 
						double* DIS,	double* TRAMAT,		int*	IDROAD, 
						int*	NROPAR, double*	ROPAR,		int*	NCHRDS, 
						char*	CHRDST, int*	LenCHRDST,	int*	N_T_SHAPE,
						double* T_SHAPE,double* UNLRAD,		double* TIREW,
						int*	NROAD,	double* EFFVOL,		double* EFFPEN, 
						double*	RCP,	double* RNORM,		double* SURFAC, 
						int*	IERR,	char*	ERRMSG,		int*	LenERRMSG), 
			int*	IDROAD, int*	NROPAR,	double*	ROPAR, 
			int*	NCHRDS, char*	CHRDST, int*	LenCHRDST, 
			double*	FORCES, double*	TORQUE, double*	DEQINI, 
			double*	DEQDER, char*	TYRMOD, int*	LenTYRMOD, 
			int*	NVARS,	double*	VARINF, int*	NWORK, 
			double*	WRKARR, int*	NIWORK, int*	IWRKAR, 
			int*	IERR){

	/*local variables */
	int i;

	double stiff_coeff;
	double damp_coeff;
	double mass;

	double x_g;
	double x_m;
	double v_g;
	double v_m;

	int umode;

	/* initial*/
	*IERR = 0;


	/*****************************************/
	/*****In case of an INITIALIZE call*******/
	if ((*JOBFLG) == init){
		c_usrmes(TRUE, "   ", 0, "INFO_NOPAD");
		c_usrmes(TRUE, "   ", 0, "INFO_NOPAD");
		c_usrmes(TRUE, "****************************************", 0, "INFO_NOPAD");
		c_usrmes(TRUE, "TYR501: Using example user tire model 369", *IDTYRE, "INFO_NOPAD");
		c_usrmes(TRUE, "****************************************",0,"INFO_NOPAD");
		c_usrmes(TRUE, "   ", 0, "INFO_NOPAD");
		c_usrmes(TRUE, "   ", 0, "INFO_NOPAD");

		printf("the logical output device number for error messages(NDEV) = %d\n", *NDEV);
		printf("the value of the USE_MODE control switch(ISWITCH) = %d\n", *ISWITCH);
		printf("value determines the action TYRSUB(JOBFLG) = %d\n", *JOBFLG);
		printf("the ID of the GFORCE statement(IDTYRE) = %d\n", *IDTYRE);
		printf("the current simulation time(TIME) = %f\n", *TIME);
		printf("wheel carrier translational displacement(DIS) \n");
		printf("size of DIS array = %d,  double = %d\n", sizeof(DIS), sizeof(DIS[0]));
		for(i=0 ; i < 3; i++)
			printf("DIS[%d] = %f\n", i, DIS[i]);
		
		printf("id tyre  = %d\n", *IDTYRE);
		printf("road file name = %s\n",CHRDST);
		printf("contact model = %d\n", *IDROAD);
		printf("dimension of the tire parameters array(TYPARR) = %d\n", *NTYPAR);
		printf("number of characters inter tire property file name(NCHTDS) = %d\n", *NCHTDS);
		printf("dimensin of ROPAR(NROPAR) = %d,  \n", *NROPAR, sizeof(*ROPAR));

		/* 
			Reading the tire property file and storing the attributes of
			tire in the arrray TYPARR. 
		*/
		rpf369(NCHTDS, CHTDST, IDTYRE, NTYPAR, TYPARR);
	}


	/**************************************/
	/*****In case of an INQUIRE call*******/
	if (*JOBFLG == inquire){
		DEQINI[0] = (double) 0.0e0;
		DEQINI[1] = (double) 0.0e0;
	}


	/*********************************************/
	/*****In case of an NORMAL or DIFF call*******/
	if( (*JOBFLG) == normal || (*JOBFLG) == diff ){
		
		// Obtaining the values from the output of tyre model
		mass = TYPARR[MASS];
		stiff_coeff = TYPARR[SPRING_CONSTANT_VERTICAL];
		damp_coeff = TYPARR[DAMPING_COEFF];

		// Calculating the displacement of ground as a step of amplitude 100 
		if(time < 0.05){
			x_g = time*200;
			v_g = 200;
		} else if (time > 0.05){
			x_g = 0.0;
			v_g = 0.0;
		}
		
		// Calculating acceleration of mass in current/next time step
		x_m = DIS[1];
		v_m = VEL[1];

		double acc_m = (damp_coeff/mass)*v_m + x_m*(damp_coeff+stiff_coeff)/mass - stiff_coeff*x_g/mass - damp_coeff*v_g/mass;

		// Calculating the forces in Y Direction.
		FORCES[0] = 0.0;
		FORCES[1] = mass*acc_m;
		FORCES[2] = 0.0;

		// Moment is zero here
		TORQUE[0] = 0.0;
		TORQUE[1] = 0.0;
		TORQUE[2] = 0.0;

	}


	/*********************************************/
	/*****In case of an END SIMULATION call*******/
	if((*JOBFLG) == endsim){
		
		printf("Mass of tyre :: %f\n", mass);
		printf("Stiffness Coefficient of spring used :: %f\n",stiff_coeff);
		printf("Damping Coefficient used :: %f\n",damp_coeff);
						
		}
	}


	/**************************************/
	/*****In case of an NO EVAL call*******/
	if((*JOBFLG) == noeval){
		return;
	}
	


	/****************************/
	/***** ERROR Handeling*******/
	if((*IERR) == 0 ){
		strcpy(TYRMOD, "TYR369 -> Example User Fiala Tire Model");
	}else{
		strcpy(TYRMOD, errmsg);
	}

}