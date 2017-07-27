#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {      
	// nlhs holds the number of output arguments
	// plhs holds references to the input arguments
	// nrhs holds the number of input arguments
	// prhs holds references to the output arguments

    int i;

    /* Examine input (right-hand-side) arguments. */

    mexPrintf("\nThere are %d right-hand-side argument(s).", nrhs);
    for (i = 0; i < nrhs; i++) {
        mexPrintf("\n\tInput Arg %i is of type:\t%s ",i,mxGetClassName(prhs[i]));
    }
    
    /* Examine output (left-hand-side) arguments. */

    mexPrintf("\n\nThere are %d left-hand-side argument(s).\n", nlhs);
    if (nlhs > nrhs)
      mexErrMsgIdAndTxt( "MATLAB:mexfunction:inputOutputMismatch",
              "Cannot specify more outputs than inputs.\n");
    
    for (i = 0; i < nlhs; i++) {
        plhs[i] = mxCreateDoubleMatrix(1, 1, mxREAL);
        *mxGetPr(plhs[i]) = (double)mxGetNumberOfElements(prhs[i]);
    }
}

