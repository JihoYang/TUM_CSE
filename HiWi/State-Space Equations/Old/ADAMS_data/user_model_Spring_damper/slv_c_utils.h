/*

	RTO_OPEN_FILE_F2C(char* TPFName, int* FileNameLen, int* errflg)
	: Utilities for Reading TeimOrbit Format perperty file

	PFName: Tire Property File Name
*/

#ifndef SOLVER_C_UTILS_H
#define SOLVER_C_UTILS_H

/* fortran calls */
#if (!defined LINKDLL_ASUTILITY)
#  if (defined _WIN32)
#   if (defined MAKE_ASUTILITY_DLL)
#     define LINKDLL_ASUTILITY __declspec(dllexport)
#   else
#     define LINKDLL_ASUTILITY __declspec(dllimport)
#   endif
#  else
#    define LINKDLL_ASUTILITY
#  endif
#endif

#if (!defined STDCALL)
#  ifdef _WIN32
#    define STDCALL __stdcall
#  else
#    define STDCALL
#  endif
#endif

#ifdef __cplusplus
extern "C" {
#endif
/*============================================================================================
 * Structures for c user subroutines
 *===========================================================================================*/
	//struct sAdamsSensor{
	//	int ID;
	//	int NPAR;
	//	const double* PAR;
	//	double VALUE;F2C
	//	double Error;
	//	char Logic[2];
	//};
/*==========================================================================================
 *
 *   UTILITY SUBROUTINES
 *
 *==========================================================================================*/
LINKDLL_ASUTILITY void c_errmes(int errflg, const char* mesage, int id, const char* endflg);
LINKDLL_ASUTILITY void c_usrmes(int msgflg, const char *mesage, int id, const char *msgtyp);
LINKDLL_ASUTILITY void c_step(double x,			double	x0, 
							  double h0,		double	x1,
							  double h1,		int		iord,
							  double *value,	int*	errflg);
LINKDLL_ASUTILITY void RTO_OPEN_FILE_F2C(char*	TPFNAME,		int* MaxStringLen, 
										 int*	FILENAMELEN,	int* SUCCESS);
LINKDLL_ASUTILITY void RTO_CLOSE_FILE_F2C(char* TPFNAME,		int* MaxStringLen,
										  int*	FILENAMELEN,	int* SUCCESS);
LINKDLL_ASUTILITY void ATRTOU(int* ID, char* UNITS);
LINKDLL_ASUTILITY void ACUNFN(char* UNITS, double* CV2MDL, double* CV2SI);
LINKDLL_ASUTILITY void RTO_READ_REAL_F2C(char*	BLOCKNAME,			int*	MaxStringLen,
										 int*	BLOCKNAMELEN,		char*	ATTRIBUTENAME, 
										 int*	AttributeMaxLen,	int*	ATTRIBUTENAMELEN, 
										 double* RETURNVALUE,		int*	SUCCESS);
LINKDLL_ASUTILITY void RTO_START_TABLE_READ_F2C(char*	BLOCKNAME,		int*	MaxStringLen, 
												int*	BLOCKNAMELEN,	char*	FORMAT, 
												int*	FormatMaxLen,	int*	FORMATLEN, 
												int*	SUCCESS);
LINKDLL_ASUTILITY void RTO_READ_TABLE_LINE_F2C(char*	VALUE,		int* MaxValueLen, 
											   int*		VALUELEN,	int* success);
LINKDLL_ASUTILITY void ACTCLC(double* TRAMAT,	double* VEL,	double* OMEGA,
							  double* OMEGAR,	double* RADIUS, double* RNORM,
							  double* VLON,		double* VCPLON,	double* VCPLAT,
							  double* VCPRVT,	double* ALPHA,	double* GAMMA,
							  double* KAPPA,	double* URAD,	double* CPMTX);
LINKDLL_ASUTILITY void ACTFZ(double* VCPVRT,	double* RADIUS,		double* TIREK,
							 double* TIREC,		double* UNLRAD,		double* FRCRAD,
							 char*	 ERRMSG,	int*	LenERRMSG,	int* IERR);
LINKDLL_ASUTILITY void XCP2HB(double* FCP,		double* TCP,	double* RAD, 
							  double* CPMTX,	double* FORCE,	double* TORQUE); 

#ifdef __cplusplus
}
#endif
#endif