#include "main.h"
#include "slv_c_utils.h"
#include "rpf369.h"

/* subroutine RPF369()
: purpose for loading the tire property file

NCHTDS: -. FileNameLen, integer
		-. Numer of characters in the file name(not the length of the FileName array)
CHTDST: -. FileName, cahracter, 256
		-. File name with full path(for example, 'usr/peploe/smith/tire.tir)
TYPARR: -. to save the tire paramter
IDTYRE: -. An integer variable that provides the ID of the GFORCE statement 
		   that applies the tire forces and moments to the wheel
*/

void rpf369(int* NCHTDS, char* CHTDST, int* IDTYRE, int* NTYPAR, double* TYPARR){

	int		i;
	int		success;
	double	return_value;
	char	units[5][12];
	double	cv2mdl[5], cv2si[5];
	double	fcvsi, mcvsi, lcvsi, tcvsi, acvsi;

	int n_nodes, arrptr;
	/*open tire property file*/
	RTO_OPEN_FILE_F2C(CHTDST, &MaxStringLen, NCHTDS, &success);
	
	if(success == 0){
		c_errmes(TRUE, "rpf501: Error opening tire property file", *IDTYRE, "STOP");
	} else {
		printf(" successed to be opened tire property file = %d\n", success);
	}

	ATRTOU(IDTYRE, units[0]);
	for(i =0; i <5; i++){
		printf(" units block units[%d] = %s\n", i, units[i]);
	}

	ACUNFN(units[0], cv2mdl, cv2si);
	fcvsi	= cv2si[0];
	mcvsi	= cv2si[1];
	lcvsi = cv2si[2];
	tcvsi	= cv2si[3];
	acvsi	= cv2si[4];

	
	/**************************  TIRPRP POPULATION  ************************/
	/******Read [MODEL] block******/
	RTO_READ_REAL_F2C("MODEL", &MaxStringLen, &ModelBlockNameLen, 
					  "USE_MODE", &MaxAttributeNameLen, &UsemodeAttributeNameLen, 
					  &return_value, &success);
	
	if(success == 0){
		c_errmes(TRUE,"rpf369: use mode undefined\n", *IDTYRE, "STOP");
	} else {
		printf(" successed to be read the [MODEL] block =%d\n", success);
	}

	/******Read [DIMENSION]Block******/
	//unloaded radius
	RTO_READ_REAL_F2C("DIMENSION", &MaxStringLen, &DimensionBlockNameLen, 
					  "UNLOADED_RADIUS", &MaxAttributeNameLen,
					  &UnloadedradiusAttributeNameLen, &return_value, &success);
	//printf(" successed to be read the [DIMENSION] block =%d\n", success);
	if(success =0){
		c_errmes(TRUE,"rpf369: unloaded radius undefined\n", *IDTYRE, "STOP");
	}
	TYPARR[UNLOADED_RADIUS] = return_value * lcvsi;
	
	
	//width
	RTO_READ_REAL_F2C("DIMENSION", &MaxStringLen, &DimensionBlockNameLen,
					  "WIDTH", &MaxAttributeNameLen, &WidthAttributeNameLen,
					  &return_value, &success);
	if(success =0){
		c_errmes(TRUE,"rpf369: width undefined\n", *IDTYRE, "STOP");
	}
	TYPARR[WIDTH] = return_value * lcvsi;

	//aspect ratio
	RTO_READ_REAL_F2C("DIMENSION", &MaxStringLen, &DimensionBlockNameLen,
					  "ASPECT_RATIO", &MaxAttributeNameLen, &AspectratioAttributeNameLen,
					  &return_value, &success);
	if(success == 0){
		c_errmes(TRUE,"rpf369: aspect ratio undefined\n", *IDTYRE, "STOP");
	}
	TYPARR[ASPECT_RATIO] = return_value;


	/******Read [PARAMETER]Block******/
	//mass
	RTO_READ_REAL_F2C("PARAMETER", &MaxStringLen, &ParameterBlockNameLen,
					  "MASS", &MaxAttributeNameLen, &MassNameLen,
					  &return_value, &success);
	if(success == 0){
		c_errmes(TRUE,"rpf369: Mass undefined\n", *IDTYRE, "STOP");
	}
	TYPARR[MASS] = return_value;

	// Spring Constant
	RTO_READ_REAL_F2C("PARAMETER", &MaxStringLen, &ParameterBlockNameLen,
					  "SPRING_CONSTANT_VERTICAL", &MaxAttributeNameLen, &SpringConstVertNameLen,
					  &return_value, &success);
	if(success == 0){
		c_errmes(TRUE,"rpf369: Spring Constant undefined\n", *IDTYRE, "STOP");
	}
	TYPARR[SPRING_CONSTANT_VERTICAL] = return_value;

	// Damping Coefficient 
	RTO_READ_REAL_F2C("PARAMETER", &MaxStringLen, &ParameterBlockNameLen,
					  "DAMPING_COEFF", &MaxAttributeNameLen, &DampingCoefNameLen,
					  &return_value, &success);
	if(success == 0){
		c_errmes(TRUE,"rpf369: Damping Coefficient undefined\n", *IDTYRE, "STOP");
	}
	TYPARR[DAMPING_COEFF] = return_value;

	

	/******Read [SHAPE]Block if it exists******/
	n_nodes = 0;
	arrptr = SHAPE;

	RTO_START_TABLE_READ_F2C("SHAPE", &MaxStringLen, &ShapeBlockNameLen, 
							Format, &MaxAttributeNameLen, &Formatlen, &success);

	if(success == TRUE){

		int MaxValueLen = 128;
		char Value[256];
		int ValueLen ;

		RTO_READ_TABLE_LINE_F2C(Value, &MaxValueLen, &ValueLen, &success);
		/*//printf("value = %s\n", Value);
		//printf("value = %c\n", Value[0]);
		//printf("value = %c\n", Value[1]);
		//printf("value = %c\n", Value[2]);
		//printf("value = %c\n", Value[3]);
		//printf("value = %c\n", Value[4]);
		//printf("value = %c\n", Value[5]);
		//printf("value = %c\n", Value[6]);
		//printf("value = %c\n", Value[7]);
		//printf("value = %c\n", Value[8]);
		//printf("value = %c\n", Value[9]);
		//printf("value = %c\n", Value[10]);
		//printf("value = %c\n", Value[11]);
		//printf("length of value = %d\n", ValueLen);*/
		if((success == TRUE) && (ValueLen > 3)){
				act_line_parse(Value, ValueLen); // in UsrFunction.c
			}
	}else{
		c_usrmes(TRUE, 
				"rpf369: No tire carcaas shape table. Cylinder will be saved", 
				*IDTYRE, "WARN");
	}
	/*************************close tire property file*********************/
	RTO_CLOSE_FILE_F2C(CHTDST, &MaxStringLen, NCHTDS, &success);
	//printf(" successed to be closed the tire property file = %d\n", success);
	if(success == 0){
		c_errmes(TRUE, "exit final: Error closing tire property file", *IDTYRE, "STOP");
	}
	/***************printf*********************
	//printf(" UNLOADED_RADIUS = %f\n",TYPARR[UNLOADED_RADIUS]);
	//printf(" WIDTH = %f\n",TYPARR[WIDTH]); 
	//printf(" ASPECT_RATIO = %f\n\n",TYPARR[ASPECT_RATIO]);

	//printf(" VERTICAL_STIFFNESS = %f\n",TYPARR[VERTICAL_STIFFNESS]);
	//printf(" VERTICAL_DAMPING = %f\n",TYPARR[VERTICAL_DAMPING]);
	//printf(" ROLLING_RESISTANCE = %f\n",TYPARR[ROLLING_RESISTANCE]);
	//printf(" CSLIP = %f\n",TYPARR[CSLIP]);
	//printf(" CALPHA = %f\n",TYPARR[CALPHA]);
	//printf(" CGAMMA = %f\n",TYPARR[CGAMMA]);
	//printf(" UMIN = %f\n",TYPARR[UMIN]);
	//printf(" UMAX = %f\n",TYPARR[UMAX]);
	//printf(" RELAXATION_LENGTH = %f\n",TYPARR[RELAXATION_LENGTH]);*/
}