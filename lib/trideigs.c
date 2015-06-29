#include "mex.h"

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


void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
{
  double *D,*E,VL,VU,*WORK,ABSTOL,*W,*Z;
  int N,IL,IU,M,*IWORK,LDZ,*IFAIL,INFO;
  mxArray *W_arr,*Z_arr;
  char *JOBZ,RANGE[2];

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

  N=mxGetNumberOfElements(prhs[0]);

  if (mxGetNumberOfElements(prhs[1])!=N-1)
    mexErrMsgTxt("The input vectors must have length N and N-1.");

  D=mxGetPr(prhs[0]);
  E=mxGetPr(prhs[1]);
  M=N; /* default value */
  RANGE[0]='A'; /* default RANGE paremeter is "A" */
  VL=0.0; VU=1.0; /* default values */
  IL=1; IU=N; /* default values */

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
		if (IL>IU) mexErrMsgTxt("Input arg 5 must be greater then arg 4.");
		M=IU-IL+1;
	} else {
		mexErrMsgTxt("Third argument mus be 'V' or 'I'");
	}
  } 
  LDZ=N;
  W_arr = mxCreateDoubleMatrix(1, M, mxREAL);
  W=mxGetPr(W_arr);
  /* Append _ for some platforms */
  if (nlhs=2) {
	  JOBZ="V";
	  Z_arr = mxCreateDoubleMatrix(LDZ,(1<M?M:1), mxREAL);
	  Z=mxGetPr(Z_arr);
  }
  else {
	  JOBZ="N";
	  Z_arr=0;
	  Z=0;
  }

  WORK=mxCalloc(5*N,sizeof(double));
  IWORK=mxCalloc(5*N,sizeof(int));
  IFAIL=mxCalloc(N,sizeof(int));
  ABSTOL=2.0*dlamch_("S");

  dstevx_(JOBZ,"I",&N,D,E,&VL,&VU,&IL,&IU,&ABSTOL,&M,W,Z,&LDZ,WORK,IWORK,IFAIL,&INFO);

  mxFree(IFAIL);
  mxFree(IWORK);
  mxFree(WORK);
  plhs[0]=W_arr;
  if (nlhs=2){
    plhs[1]=Z_arr;
  }
  if (INFO)
    mexErrMsgTxt("No convergence in DSTEQR.");
}
