#include "mex.h"
#include "stdio.h"

/*   TRIDEIGs Eigenvalues an egien vectors of symmetric tridiagonal matrix.   */

/*   V. Yordanov <v.yordanov@phys.uni-sofia.bg>   */
/*   Department of Optics and Spectroscopy, Sofia University    */
/*   May, 2009                  */
/*  [W] = trideigs(D,E) */
/*  [W,Z] = trideigs(D,E) */
/*  [W,Z] = trideigs(D,E,RANGE,VL,VU) */
/*  [W,Z] = trideigs(D,E,RANGE,IL,IU) */

/* Quick install */
/* For windows: */
/* %MATLAB\bin\mex trideigs.c %MATLAB\extern\lib\win32\lcc\libmwlapack.lib */
/* For linux: */
/* $MATLAB/bin/mex -llapack trideigs.c */
/* Add the directory that contains compiled mex file to your path. */
/* Add tridiegs.m to the path also. */

extern void dstevr_(char* JOBZ, char *RANGE, int* N, double* D, double* E, 
        double* VL, double* VU, int* IL, int* IU, double* ABSTOL, int* M, 
        double* W, double* Z, int* LDZ, int* ISUPPZ, double* WORK, 
        int* LWORK, int* IWORK, int* LIWORK, int* INFO);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
  if (nrhs>5)
    mexErrMsgTxt("Must have max 5 input arguments.");
  if (nrhs<2)
    mexErrMsgTxt("Must have min 2 input arguments.");
  if (nlhs>2)
    mexErrMsgTxt("Too many output arguments.");
  if (!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]))
    mexErrMsgTxt("Both inputs must be double arrays.");
  if (mxIsComplex(prhs[0]) || mxIsComplex(prhs[1]))
    mexErrMsgTxt("Both inputs must be real.");
  if (mxGetM(prhs[0])!=1 && mxGetN(prhs[0])!=1)
    mexErrMsgTxt("First input must be a vector.");
  if (mxGetM(prhs[1])!=1 && mxGetN(prhs[1])!=1)
    mexErrMsgTxt("Second input must be a vector.");
  if (mxGetM(prhs[1])!=1 && mxGetN(prhs[1])!=1)
    mexErrMsgTxt("Second input must be a vector.");
  
  int N=mxGetNumberOfElements(prhs[0]);
  if (mxGetNumberOfElements(prhs[1])!=N-1)
    mexErrMsgTxt("The input vectors must have length N and N-1.");
  
  /* default values */
  char *JOBZ, *RANGE="A";
  double VL=0.0, VU=1.0, ABSTOL=0.0;
  int IL=1, IU=N, M=N, LDZ=N, LWORK=20*N, LIWORK=10*N, INFO=0;
  
  if (nrhs > 3) { /* we have RANGE parameter specified */
    if (mxIsChar(prhs[2]) != 1) mexErrMsgTxt("Input must be a string.");
	if (nrhs < 5) 
		mexErrMsgTxt("If range is specified then must have 5 input arguments");
	if (!(mxGetM(prhs[3]) == 1 && mxGetN(prhs[3]) == 1)) {
		mexErrMsgTxt("Input arg 4 must be scalar");
	}
	if (!(mxGetM(prhs[4]) == 1 && mxGetN(prhs[4]) == 1)) {
		mexErrMsgTxt("Input arg 5 must be scalar");
	}
	mxGetString(prhs[2], RANGE, 2);
	
    if (RANGE[0]=='V') {
		if (!mxIsDouble(prhs[3])) {
			mexErrMsgTxt("Input arg 4 must be number of type double.");
		}
		if (!mxIsDouble(prhs[4])) {
			mexErrMsgTxt("Input arg 5 must be number of type double.");
		}
		VL=mxGetScalar(prhs[3]);   
		VU=mxGetScalar(prhs[4]);
    } else if (RANGE[0]=='I'){
		if (!mxIsDouble(prhs[3])) {
			mexErrMsgTxt("Input arg 4 must be number of type unsigned int.");
		}
		if (!mxIsDouble(prhs[4])) {
			mexErrMsgTxt("Input arg 5 must be number of type unsigned int.");
		}		
		IL=(int)mxGetScalar(prhs[3]);   
		IU=(int)mxGetScalar(prhs[4]);
		if (IL>IU) mexErrMsgTxt("Input arg 5 must be greater than arg 4.");
		M=IU-IL+1;
	} else {
		mexErrMsgTxt("Third argument mus be 'V' or 'I'");
	}
  }
  double *D=mxGetPr(prhs[0]);
  double *E=mxGetPr(prhs[1]);
  int *ISUPPZ=(int*)malloc(2*M*sizeof(int));
  
  double *W, *Z;
  mxArray *W_arr, *Z_arr;
  W_arr=mxCreateDoubleMatrix(1, M, mxREAL);
  W=mxGetPr(W_arr);
  if(nlhs==2){
	  JOBZ="V";
	  Z_arr=mxCreateDoubleMatrix(LDZ, (1<M?M:1), mxREAL);
	  Z=mxGetPr(Z_arr);
  }else{
	  JOBZ="N";
	  Z_arr=0;
	  Z=0;
  }
  
  /* Append _ for UNIX platforms */
  int *IWORK=(int*)malloc(LIWORK*sizeof(int));
  double *WORK=(double*)malloc(LWORK*sizeof(double));
  dstevr_(JOBZ,RANGE,&N,D,E,&VL,&VU,&IL,&IU,&ABSTOL,&M,W,Z,&LDZ,ISUPPZ,WORK,
          &LWORK,IWORK,&LIWORK,&INFO);
  
  free((void*)WORK);
  free((void*)IWORK);

  plhs[0]=W_arr;
  if(nlhs==2){
    plhs[1]=Z_arr;
  }
  
  if(INFO){      
      char err[64];
      sprintf(err, "dstevr INFO=%d", INFO);
      mexErrMsgTxt(err);
  }
}
