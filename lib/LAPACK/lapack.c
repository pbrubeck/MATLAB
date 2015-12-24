/**************************************************************************/
/* Copyright 2007-2009, Timothy Toolan, University of Rhode Island        */
/*                                                                        */
/* To build on Windows, run:                                              */
/*   mex lapack.c                                                         */
/*                                                                        */
/* To build on any other operating system, run:                           */
/*   mex -ldl lapack.c                                                    */
/*                                                                        */
/* For details on use and examples, in Matlab type:                       */
/*   help lapack                                                          */
/*                                                                        */
/**************************************************************************/

#include "mex.h"
#include <stdio.h>   /* for sprintf */
#include <ctype.h>   /* for tolower, toupper, isupper */
#include <string.h>  /* for strcpy, strncat, strcat, strcmp */
#include <stdlib.h>  /* for bsearch */

/* note that _WIN32 should be defined in 64 bit windows also */
#if defined(_WIN32)
#include <windows.h> /* needed for HISTANCE */
#else
#include <dlfcn.h>   /* needed for dlsym and RTLD_DEFAULT */
#endif

/* Maximum number of arguments that can be specified for a function.
   The most arguments to any function is 29 in LAPACK 3.1.1 (see
   xGGEVX), and 31 in LAPACK 3.2.1 (see xLA_GBRFSX_EXTENDED). */
#define MAX_ARGS 32

/* Maximum number of arguments of type CHARACTER for any function.
   The most is 4 in LAPACK 3.1.1 and 6 in LAPACK 3.2.1.  */
#define MAX_CHAR_ARGS 6

/* These are function pointer prototypes for functions which return
   one of seven types, and have MAX_ARGS+MAX_CHAR_ARGS void*
   arguments.  Because of the way FORTRAN handles strings, we will
   append an argument of the string length for each string argument at
   the end of the argument list.  Although this is not guaranteed to
   be compatible with all FORTRAN compilers, it seems to work with the
   ones Matlab uses. */
typedef void* v;
typedef struct { float r; float i; } cfloat;
typedef struct { double r; double i; } cdouble;

typedef void (*funcV)
  (v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v);
typedef mwSignedIndex (*funcI)
  (v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v);
typedef int (*funcL)
  (v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v);
typedef float (*funcS)
  (v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v);
typedef double (*funcD)
  (v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v);
typedef cfloat (*funcC)
  (v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v);
typedef cdouble (*funcZ)
  (v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v);
typedef char (*funcH)
  (v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v);

/****************************************************************************
  This function will look for a function whose name is the lowercase
  version of szRawFuncName in the LAPACK and BLAS libraries. 
  This is the only function that is operating system dependent.
  Arguments:
    nLen - the length of the string in szRawFuncName.
    szRawFuncName - a string containing the name of the desired
                    function.  This string is not null terminated.
  Return Value: 
    the function pointer, or 0 if no function is found.
 ****************************************************************************/
funcV GetFunctionPointer(int nLen, char *szRawFuncName)
{
  int i;
  funcV pFunc;
#if defined(_WIN32)
  HINSTANCE hinstLib; 
#endif

  /* convert function to all lowercase, then append an underscore if UNIX */
  char *szFuncName = mxCalloc(nLen+2, sizeof(char));
  for (i=0; i<nLen; i++)
    szFuncName[i] = tolower(szRawFuncName[i]);
#if !defined(_WIN32)
  szFuncName[i++] = '_';
#endif
  szFuncName[i] = '\0';

#if defined(_WIN32)
  /* on windows, explicitly look in the LAPACK and BLAS libraries */
  if ((hinstLib = LoadLibrary("libmwlapack")))
    pFunc = (funcV)GetProcAddress(hinstLib, szFuncName);
  if (!pFunc && (hinstLib = LoadLibrary("libmwblas")))
    pFunc = (funcV)GetProcAddress(hinstLib, szFuncName);
#else
  /* on UNIX, use the library search path */
  pFunc = dlsym(RTLD_DEFAULT, szFuncName);
#endif

  return pFunc;
}

/* This prototype is needed for the ParsePrototype function */
char *SearchForProto(const char *szFuncName);

/****************************************************************************
  This function will parse the LAPACK/BLAS function prototype in argument 0.
  Arguments:
    nProtoOnly - a boolean indicating that we only want to lookup the
      internal function prototype and return it.
    pszProto - a pointer to a string which contains argument 0.  If this
      is just the name of the function without argument types or return
      type, it will be updated to point to the internal prototype for
      the specified function.
  Return Arguments:
    numArgOut - a pointer to an integer which will be set to the number
      of return arguments requested by the function prototype.  This 
      will be the number of uppercase type letters.
    numArgIn - a pointer to an integer which will be set to the number
      of arguments that are required according to the function prototype.
    retType - a pointer to a char which will be set to the type of the 
      return value, or 'V' if there is no return value.
    argType - an array of numArgOut chars which will be set to the type 
      of the corresponding argument.
    indArgOut - an array of numArgOut integers which will be set to 
      either the index within the output arguments that a given input 
      argument corresponds to, or -1 if that argument should not be
      an output argument.
  Return Value:
    a pointer to the requested BLAS/LAPACK function, or NULL if
    only the internal prototype is requested by nProtoOnly = 1.
 ****************************************************************************/
funcV ParsePrototype(int nProtoOnly, char **pszProto, 
                     int *numArgOut, int *numArgIn, 
                     char *retType, char *argType, int indArgOut[])
{
  funcV pFunc;
  char c, *szTmpProto, *szRawFuncName, *szFoundProto, *szProto;
  int nLen = 0, nRawFuncLen = 0, argOutCt, i;
  char msg[128];
  int nOpenParen = 0, nEqualSign = 0;

  /* squeeze out all space */
  szProto = *pszProto;
  szTmpProto = szProto;
  while (*szTmpProto){
    if (!isspace(*szTmpProto)){
      szProto[nLen++] = *szTmpProto;
      if (*szTmpProto == '(')
        nOpenParen = 1;
      if (*szTmpProto == '=')
        nEqualSign = 1;
    }
    szTmpProto++;
  }
  szProto[nLen] = '\0';
  if (nLen < 3)
    mexErrMsgTxt("Invalid function prototype.");
  
  /* if no parentheses, we should search for the built in prototype string */
  if (!nEqualSign && !nOpenParen && (szProto[nLen-1] != ')')){
    for (i=0; i<nLen; i++)
      szProto[i] = toupper(szProto[i]);
    szFoundProto = SearchForProto(szProto);
    if (!szFoundProto)
      mexErrMsgTxt("Prototype not found.");
    nLen = strlen(szFoundProto);
    *pszProto = szProto = mxCalloc(nLen+1, sizeof(char));
    strcpy(szProto, szFoundProto);

    /* check if all we want is the prototype string */
    if (nProtoOnly)
      return 0;
  }

  if ((nLen < 3) || (szProto[nLen-1] != ')'))
    mexErrMsgTxt("Invalid function prototype.");
  
  /* check if there is a return value */
  *retType = 'V';
  if (szProto[1] == '='){
    c = *szProto;
    if ((c == 'I') || (c == 'S') || (c == 'D') || 
        (c == 'C') || (c == 'Z') || (c == 'L') || (c == 'H')){
      *retType = *szProto;
      szProto += 2;
    }
    else
      mexErrMsgTxt("Invalid return type in function prototype, "
                   "valid return types are the letters: SDCZILH.");
  }
  
  /* determine function name string location, then get pointer */
  szRawFuncName = szProto;
  while (*szProto && (*szProto != '(')){
    nRawFuncLen++;
    szProto++;
  }
  if (*szProto++ != '(')
    mexErrMsgTxt("No open parenthesis found in function prototype.");
  szRawFuncName[nRawFuncLen] = '\0';
  if (!(pFunc = GetFunctionPointer(nRawFuncLen, szRawFuncName))){
    strcpy(msg, "Cannot find BLAS/LAPACK library function \"");
    strncat(msg, szRawFuncName, sizeof(msg) - 50);
    strcat(msg, "\".");
    mexErrMsgTxt(msg);
  }

  /* check the argument types */
  *numArgOut = *numArgIn = 0;
  argOutCt = (*retType != 'V') ? 1 : 0;
  while (*szProto){
    c = toupper(*szProto);
    if ((c == 'I') || (c == 'S') || (c == 'D') || 
        (c == 'C') || (c == 'Z') || (c == 'L') || (c == 'H')){
      indArgOut[*numArgIn] = -1;
      if (isupper(*szProto))
        indArgOut[*numArgIn] = argOutCt++;
      argType[(*numArgIn)++] = *szProto;
      if (isupper(*szProto))
        (*numArgOut)++;
      szProto++;
      if (*szProto == ')')
        break;
      else if (*szProto == ',')
        szProto++;
      else{
        sprintf(msg, "Argument %d is invalid in function prototype, "
                "valid types are the letters: sSdDcCzZiIlLhH.", *numArgIn);
        mexErrMsgTxt(msg);
      }
    }
    else{
      sprintf(msg, "Argument %d is invalid in function prototype, "
              "valid types are the letters: sSdDcCzZiIlLhH.", *numArgIn+1);
      mexErrMsgTxt(msg);
    }
  }
  
  return pFunc;
}

/****************************************************************************
  This function makes sure that the arguments are of the correct
  type for the BLAS/LAPACK routine, and converts them if necessary.
  All arguments will be copied to new memory.  Complex arguments that
  are returned will require double the space allocated due to the
  difference in storage between FORTRAN and Matlab.
  Arguments:
    nNumArgs - the number of input arguments that need to be copied.
    argType - a length nNumArgs array of chars, which are the type
      of the corresponding argument in "args".
    args - a length nNumArgs array of pointers to mxArray objects, 
      which are the input arguments.
    indArgOut - a length nNumArgs array of integers which are
      either the index within the output arguments that a given input 
      argument corresponds to, or -1 if that input argument should not 
      be an output argument.
  Return Arguments:
    a - a length nNumArgs+(the number of type CHARACTER elements) array 
      of pointers, which are the actual arguments that get sent to the 
      BLAS/LAPACK function.  The first nNumArgs elements are the
      converted elements from "args".   The length of each string 
      argument is appended to the end.
    plhs - this is the array of return values for this mex function.
      These will be mxArray objects that share the memory with
      the corresponding element in "a" when possible.  All return 
      memory is allocated at this time.  If plhs[0] is a cell
      object, then all the mxArray objects will be placed in this cell,
      otherwise they will be placed in the elements of plhs.
 ****************************************************************************/
void CopyArgs(int nNumArgs, char *argType, const mxArray *args[],
              void *a[], mxArray *plhs[], int indArgOut[])
{
  char msg[128];
  int i, nNumChars = 0;
  mwSignedIndex j, m, n, mn;
  void *srcr, *srci;
  char c;
  mxArray *value;

  /* copy all of the argument pointers to an array and convert if necessary */
  for (i=0; i<nNumArgs; i++){
    m = mxGetM(args[i]);
    n = mxGetN(args[i]);
    mn = m*n;
    c = toupper(argType[i]);
    if (c == 'I'){
      /* INTEGER *******************************************************/
      /* can be: int32, int64, double, logical */
      if (isupper(argType[i])){
        value = mxCreateNumericMatrix(m, n, 
             ((sizeof(mwSignedIndex) == 8) ? mxINT64_CLASS : mxINT32_CLASS), 
                 mxREAL);
        a[i] = mxGetPr(value);
        if (plhs[0] && mxIsCell(plhs[0]))
          mxSetCell(plhs[0], indArgOut[i], value);
        else
          plhs[indArgOut[i]] = value;
      }
      else
        a[i] = mxCalloc(mn, sizeof(mwSignedIndex));
      srcr = mxGetPr(args[i]);
      if (mxIsInt32(args[i]))
        for (j=0; j<mn; j++)
          ((mwSignedIndex*)(a[i]))[j] = ((int*)srcr)[j];
      else if (mxIsInt64(args[i]))
        for (j=0; j<mn; j++)
          ((mwSignedIndex*)(a[i]))[j] = ((long long*)srcr)[j];
      else if (mxIsDouble(args[i]) && !mxIsComplex(args[i]))
        for (j=0; j<mn; j++)
          ((mwSignedIndex*)(a[i]))[j] = ((double*)srcr)[j];
      else if (mxIsLogical(args[i]))
        for (j=0; j<mn; j++)
          ((mwSignedIndex*)(a[i]))[j] = ((char*)srcr)[j];
      else{
        sprintf(msg, "Can't convert argument %d to FORTRAN type "
                "INTEGER, use double, int32, int64 or logical.", i+1);
        mexErrMsgTxt(msg);
      }
    }
    else if (c == 'L'){
      /* LOGICAL *******************************************************/
      /* can be: int32, int64, double, logical */
      if (isupper(argType[i])){
        value = mxCreateNumericMatrix(m, n, mxINT32_CLASS, mxREAL);
        a[i] = mxGetPr(value);
        if (plhs[0] && mxIsCell(plhs[0]))
          mxSetCell(plhs[0], indArgOut[i], value);
        else
          plhs[indArgOut[i]] = value;
      }
      else
        a[i] = mxCalloc(mn, sizeof(int));
      srcr = mxGetPr(args[i]);
      if (mxIsInt32(args[i]))
        for (j=0; j<mn; j++)
          ((int*)(a[i]))[j] = ((int*)srcr)[j];
      else if (mxIsInt64(args[i]))
        for (j=0; j<mn; j++)
          ((int*)(a[i]))[j] = ((long long*)srcr)[j];
      else if (mxIsDouble(args[i]) && !mxIsComplex(args[i]))
        for (j=0; j<mn; j++)
          ((int*)(a[i]))[j] = ((double*)srcr)[j];
      else if (mxIsLogical(args[i]))
        for (j=0; j<mn; j++)
          ((int*)(a[i]))[j] = ((char*)srcr)[j];
      else{
        sprintf(msg, "Can't convert argument %d to FORTRAN type "
                "LOGICAL, use double, int32, int64 or logical.", i+1);
        mexErrMsgTxt(msg);
      }
    }
    else if (c == 'D'){
      /* DOUBLE PRECISION **********************************************/
      /* can be: double */
      if (isupper(argType[i])){
        value = mxCreateNumericMatrix(m, n, mxDOUBLE_CLASS, mxREAL);
        a[i] = mxGetPr(value);
        if (plhs[0] && mxIsCell(plhs[0]))
          mxSetCell(plhs[0], indArgOut[i], value);
        else
          plhs[indArgOut[i]] = value;
      }
      else
        a[i] = mxCalloc(mn, sizeof(double));
      srcr = mxGetPr(args[i]);
      if (mxIsDouble(args[i]) && !mxIsComplex(args[i]))
        for (j=0; j<mn; j++)
          ((double*)a[i])[j] = ((double*)srcr)[j];
      else{
        sprintf(msg, "Can't convert argument %d to FORTRAN type "
                "DOUBLE PRECISION, use double.", i+1);
        mexErrMsgTxt(msg);
      }
    }
    else if (c == 'Z'){
      /* COMPLEX*16 ****************************************************/
      /* can be: double, complex double */
      if (isupper(argType[i])){
        value = mxCreateNumericMatrix(m, n, mxDOUBLE_CLASS, mxCOMPLEX);
        if (plhs[0] && mxIsCell(plhs[0]))
          mxSetCell(plhs[0], indArgOut[i], value);
        else
          plhs[indArgOut[i]] = value;
      }
      a[i] = mxCalloc(2*mn, sizeof(double));
      srcr = mxGetPr(args[i]);
      srci = mxGetPi(args[i]);
      if (mxIsDouble(args[i]) && !mxIsComplex(args[i]))
        for (j=0; j<mn; j++){
          ((double*)a[i])[2*j] = ((double*)srcr)[j];
          ((double*)a[i])[2*j+1] = 0;
        }
      else if (mxIsDouble(args[i]) && mxIsComplex(args[i]))
        for (j=0; j<mn; j++){
          ((double*)a[i])[2*j] = ((double*)srcr)[j];
          ((double*)a[i])[2*j+1] = ((double*)srci)[j];
        }
      else{
        sprintf(msg, "Can't convert argument %d to FORTRAN type "
                "COMPLEX*16, use double or complex double", i+1);
        mexErrMsgTxt(msg);
      }
    }
    else if (c == 'S'){
      /* REAL **********************************************************/
      /* can be: single, double */
      if (isupper(argType[i])){
        value = mxCreateNumericMatrix(m, n, mxSINGLE_CLASS, mxREAL);
        a[i] = mxGetPr(value);
        if (plhs[0] && mxIsCell(plhs[0]))
          mxSetCell(plhs[0], indArgOut[i], value);
        else
          plhs[indArgOut[i]] = value;
      }
      else
        a[i] = mxCalloc(mn, sizeof(float));
      srcr = mxGetPr(args[i]);
      if (mxIsSingle(args[i]) && !mxIsComplex(args[i]))
        for (j=0; j<mn; j++)
          ((float*)a[i])[j] = ((float*)srcr)[j];
      else if (mxIsDouble(args[i]) && !mxIsComplex(args[i]))
        for (j=0; j<mn; j++)
          ((float*)a[i])[j] = ((double*)srcr)[j];
      else{
        sprintf(msg, "Can't convert argument %d to FORTRAN type "
                "REAL, use double or single.", i+1);
        mexErrMsgTxt(msg);
      }
    }
    else if (c == 'C'){
      /* COMPLEX *******************************************************/
      /* can be: single, complex single, double, complex double */
      if (isupper(argType[i])){
        value = mxCreateNumericMatrix(m, n, mxSINGLE_CLASS, mxCOMPLEX);
        if (plhs[0] && mxIsCell(plhs[0]))
          mxSetCell(plhs[0], indArgOut[i], value);
        else
          plhs[indArgOut[i]] = value;
      }
      a[i] = mxCalloc(2*mn, sizeof(float));
      srcr = mxGetPr(args[i]);
      srci = mxGetPi(args[i]);
      if (mxIsSingle(args[i]) && !mxIsComplex(args[i]))
        for (j=0; j<mn; j++){
          ((float*)a[i])[2*j] = ((float*)srcr)[j];
          ((float*)a[i])[2*j+1] = 0;
        }
      else if (mxIsSingle(args[i]) && mxIsComplex(args[i]))
        for (j=0; j<mn; j++){
          ((float*)a[i])[2*j] = ((float*)srcr)[j];
          ((float*)a[i])[2*j+1] = ((float*)srci)[j];
        }
      else if (mxIsDouble(args[i]) && !mxIsComplex(args[i]))
        for (j=0; j<mn; j++){
          ((float*)a[i])[2*j] = ((double*)srcr)[j];
          ((float*)a[i])[2*j+1] = 0;
        }
      else if (mxIsDouble(args[i]) && mxIsComplex(args[i]))
        for (j=0; j<mn; j++){
          ((float*)a[i])[2*j] = ((double*)srcr)[j];
          ((float*)a[i])[2*j+1] = ((double*)srci)[j];
        }
      else{
        sprintf(msg, "Can't convert argument %d to FORTRAN type "
                "COMPLEX, use double or single", i+1);
        mexErrMsgTxt(msg);
      }
    }
    else if (c == 'H'){
      /* CHARACTER *****************************************************/
      /* can be: char */
      if (mxIsChar(args[i])){
        a[i] = mxArrayToString(args[i]);
        if (isupper(argType[i])){
          value = mxCreateString(a[i]);
          if (plhs[0] && mxIsCell(plhs[0]))
            mxSetCell(plhs[0], indArgOut[i], value);
          else
            plhs[indArgOut[i]] = value;
        }
        /* Append the lengths of string to the end of the argument list */
        if (nNumChars >= MAX_CHAR_ARGS){
          sprintf(msg, "There cannot be more than %d arguments of type "
                  "CHARACTER (type 'h' or 'H') in the function prototype.",
                  MAX_CHAR_ARGS);
          mexErrMsgTxt(msg);
         }
        *(mwSignedIndex*)&(a[nNumArgs+nNumChars]) = strlen(a[i]);
        nNumChars++;
      }
      else{
        sprintf(msg, "Can't convert argument %d to FORTRAN type "
                "CHARACTER, use a string.", i+1);
        mexErrMsgTxt(msg);
      }
    }
  }
}

/****************************************************************************
   This function will call the BLAS/LAPACK function that is desired.
   Arguments:
     retType - the type that is returned by the BLAS/LAPACK function.
     pFunc - a pointer to the BLAS/LAPACK function.
     a - an array of arguments to be passed to the BLAS/LAPACK function.
   Return Value:
     a mxArray that contains the return value.
 ****************************************************************************/
mxArray *EvaluateFunction(char retType, funcV pFunc, void *a[])
{
  mxArray *retVal = 0;

  if (retType == 'V'){
    pFunc
      (a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9],
       a[10],a[11],a[12],a[13],a[14],a[15],a[16],a[17],a[18],a[19],
       a[20],a[21],a[22],a[23],a[24],a[25],a[26],a[27],a[28],a[29],
       a[30],a[31],a[32],a[33],a[34],a[35],a[36],a[37]);
  }
  else if (retType == 'I'){
    mwSignedIndex result = ((funcI)pFunc)
      (a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9],
       a[10],a[11],a[12],a[13],a[14],a[15],a[16],a[17],a[18],a[19],
       a[20],a[21],a[22],a[23],a[24],a[25],a[26],a[27],a[28],a[29],
       a[30],a[31],a[32],a[33],a[34],a[35],a[36],a[37]);

    retVal = mxCreateNumericMatrix(1, 1, 
       ((sizeof(mwSignedIndex) == 8) ? mxINT64_CLASS : mxINT32_CLASS), mxREAL);
    *(mwSignedIndex*)(mxGetPr(retVal)) = result;
  }
  else if (retType == 'L'){
    int result = ((funcL)pFunc)
      (a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9],
       a[10],a[11],a[12],a[13],a[14],a[15],a[16],a[17],a[18],a[19],
       a[20],a[21],a[22],a[23],a[24],a[25],a[26],a[27],a[28],a[29],
       a[30],a[31],a[32],a[33],a[34],a[35],a[36],a[37]);

    retVal = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
    *(int*)(mxGetPr(retVal)) = result;
  }
  else if (retType == 'D'){
    double result = ((funcD)pFunc)
      (a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9],
       a[10],a[11],a[12],a[13],a[14],a[15],a[16],a[17],a[18],a[19],
       a[20],a[21],a[22],a[23],a[24],a[25],a[26],a[27],a[28],a[29],
       a[30],a[31],a[32],a[33],a[34],a[35],a[36],a[37]);

    retVal = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
    *(mxGetPr(retVal)) = result;
  }
  else if (retType == 'Z'){
    cdouble result = ((funcZ)pFunc)
      (a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9],
       a[10],a[11],a[12],a[13],a[14],a[15],a[16],a[17],a[18],a[19],
       a[20],a[21],a[22],a[23],a[24],a[25],a[26],a[27],a[28],a[29],
       a[30],a[31],a[32],a[33],a[34],a[35],a[36],a[37]);

    retVal = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxCOMPLEX);
    *(mxGetPr(retVal)) = result.r;
    *(mxGetPi(retVal)) = result.i;
  }
  else if (retType == 'S'){
    float result = ((funcS)pFunc)
      (a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9],
       a[10],a[11],a[12],a[13],a[14],a[15],a[16],a[17],a[18],a[19],
       a[20],a[21],a[22],a[23],a[24],a[25],a[26],a[27],a[28],a[29],
       a[30],a[31],a[32],a[33],a[34],a[35],a[36],a[37]);

    retVal = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
    *(float*)(mxGetPr(retVal)) = result;
  }
  else if (retType == 'C'){
#ifndef __LCC__
    cfloat result = ((funcC)pFunc) 
      (a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9],
       a[10],a[11],a[12],a[13],a[14],a[15],a[16],a[17],a[18],a[19],
       a[20],a[21],a[22],a[23],a[24],a[25],a[26],a[27],a[28],a[29],
       a[30],a[31],a[32],a[33],a[34],a[35],a[36],a[37]);
#else
    /* Because calling a routine that returns a fortran COMPLEX  */
    /* crashes the lcc compiler on windows, this is an incredible */
    /* hack job which seems to get around this problem. */
    cdouble resultx = ((funcZ)pFunc) 
      (a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9],
       a[10],a[11],a[12],a[13],a[14],a[15],a[16],a[17],a[18],a[19],
       a[20],a[21],a[22],a[23],a[24],a[25],a[26],a[27],a[28],a[29],
       a[30],a[31],a[32],a[33],a[34],a[35],a[36],a[37]);
    cfloat result = *(cfloat*)&(resultx.r);
#endif

    retVal = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxCOMPLEX);
    *(float*)(mxGetPr(retVal)) = result.r;
    *(float*)(mxGetPi(retVal)) = result.i;
  } 
  else if (retType == 'H'){
    char result[2];
    result[0] = ((funcH)pFunc)
      (a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9],
       a[10],a[11],a[12],a[13],a[14],a[15],a[16],a[17],a[18],a[19],
       a[20],a[21],a[22],a[23],a[24],a[25],a[26],a[27],a[28],a[29],
       a[30],a[31],a[32],a[33],a[34],a[35],a[36],a[37]);

    result[1] = '\0';
    retVal = mxCreateString(result);
  }

  return retVal;
}

/****************************************************************************
   This will convert complex arguments back to Matlab format and insert
   the first character of any strings into the return string.
   Arguments: 
     nNumArgs - the number of input arguments in 'a'.
     argType - a length nNumArgs array of chars, which are the type
       of the corresponding argument in "a".
     a - a length nNumArgs array of pointers that are the possibly
       modified results from the BLAS/LAPACK function.
     plhs - this is the array of return values for this mex function.
     indArgOut - a length nNumArgs array of integers which are
       either the index within the output arguments that a given input 
       argument corresponds to, or -1 if that input argument should not 
       be an output argument.
 ****************************************************************************/
void RestoreArgs(int nNumArgs, char *argType,
                 void *a[], mxArray *plhs[], int indArgOut[])
{
  int i;
  mwSignedIndex j, mn;
  void *destr, *desti;
  mxArray *value;
  mxChar *pChar;

  for (i=0; i<nNumArgs; i++){
    if (argType[i] == 'Z'){
      value = (plhs[0] && mxIsCell(plhs[0])) ? 
        mxGetCell(plhs[0], indArgOut[i]) : plhs[indArgOut[i]];
      destr = mxGetPr(value);
      desti = mxGetPi(value);
      mn = mxGetM(value) * mxGetN(value);
      for (j=0; j<mn; j++){
        ((double*)destr)[j] = ((double*)a[i])[2*j];
        ((double*)desti)[j] = ((double*)a[i])[2*j+1];
      }
    }
    else if (argType[i] == 'C'){
      value = (plhs[0] && mxIsCell(plhs[0])) ? 
        mxGetCell(plhs[0], indArgOut[i]) : plhs[indArgOut[i]];
      destr = mxGetPr(value);
      desti = mxGetPi(value);
      mn = mxGetM(value) * mxGetN(value);
      for (j=0; j<mn; j++){
        ((float*)destr)[j] = ((float*)a[i])[2*j];
        ((float*)desti)[j] = ((float*)a[i])[2*j+1];
      }
    }
    else if (argType[i] == 'H'){
      value = (plhs[0] && mxIsCell(plhs[0])) ? 
        mxGetCell(plhs[0], indArgOut[i]) : plhs[indArgOut[i]];
      if ((pChar = mxGetChars(value)))
        *pChar = ((char*)a[i])[0];
    }
  }
}

/****************************************************************************
  This is the main function for this mex program.
 ****************************************************************************/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  funcV pFunc;
  int numArgOut, numArgIn, indArgOut[MAX_ARGS], numReturn;
  char retType, argType[MAX_ARGS], *szProto;
  void *a[MAX_ARGS+MAX_CHAR_ARGS];
  mxArray *retVal;
  char msg[128];

  if ((nrhs < 1) || (mxIsChar(prhs[0]) != 1) || (mxGetM(prhs[0]) != 1))
    mexErrMsgTxt("First argument must be the BLAS/LAPACK function prototype.");

  /* parse the function prototype and get function pointer */
  szProto = mxArrayToString(prhs[0]);
  pFunc = ParsePrototype((nrhs == 1), &szProto, &numArgOut, &numArgIn, 
                         &retType, argType, indArgOut);

  /* if pFunc is NULL, all we wanted is the prototype string, so return it */
  if (!pFunc){
    if (nlhs > 1){
      sprintf(msg, "Only the function prototype was requested, "
              "but %d output arguments specified.", nlhs);
      mexErrMsgTxt(msg);
    }
    plhs[0] = mxCreateString(szProto);
    return;
  }

  /* check if number of input arguments is okay */
  if (numArgIn != nrhs-1){
    sprintf(msg, "Function prototype specifies %d arguments, "
            "but %d arguments given.", numArgIn, nrhs-1);
    mexErrMsgTxt(msg);
  }
  if (numArgIn > MAX_ARGS){
    sprintf(msg, "There cannot be more than %d BLAS/LAPACK arguments.", 
            MAX_ARGS);
    mexErrMsgTxt(msg);
  }

  /* check if number of output arguments is okay */
  numReturn = numArgOut + ((retType == 'V') ? 0 : 1);
  if ((nlhs > 1) && (numReturn != nlhs)){
    if (retType == 'V')
      sprintf(msg, "Function prototype specifies %d output arguments, "
              "but %d output arguments specified.", numArgOut, nlhs);
    else
      sprintf(msg, "Function prototype specifies a return value plus "
              "%d output arguments, but %d output arguments specified.", 
              numArgOut, nlhs);
    mexErrMsgTxt(msg);
  }
  /* there is more than one return value requested but only one
     return value so use a cell array */
  if ((numReturn > 1) && (nlhs <= 1))
    plhs[0] = mxCreateCellMatrix(numReturn, 1);

  /* copy all of the argument pointers to an array and convert if necessary */
  CopyArgs(numArgIn, argType, &prhs[1], a, plhs, indArgOut);

  /* call the BLAS/LAPACK function */
  retVal = EvaluateFunction(retType, pFunc, a);
  if (retType != 'V'){
    if (plhs[0] && mxIsCell(plhs[0]))
      mxSetCell(plhs[0], 0, retVal);
    else
      plhs[0] = retVal;
  }

  /* copy complex output values back to plhs */
  RestoreArgs(numArgIn, argType, a, plhs, indArgOut);

  return;
}

/****************************************************************************
  This is the comparison function for the bsearch call in SearchForProto.
  Arguments:
    pa - pointer to a null terminated string.
    pb - pointer to a null or open parenthesis terminated string,
         which starts at *pb[0] if *pb[1] is not an equal sign and
         at *pb[2] if *pb[1] is an equal sign.
 ****************************************************************************/
static int Compare(const void *pa, const void *pb)
{
  int result = 0;
  const unsigned char *a = *(unsigned char**)pa, *b = *(unsigned char**)pb;
  if (b[1] == '=')
    b += 2;
  
  while (!(result = *a - *b) && *a){
    a += 1;
    b += 1;
  }

  return (*a == '\0' && *b == '(') ? 0 : result;
}

/****************************************************************************
  This function will search for a prototype of the given BLAS/LAPACK
  function.
  Arguments:
    szFuncName - an uppercase string containing the name of the
      BLAS/LAPACK function to search for.
  Returns:
    If the function is found, a string containing the prototype, 
    otherwise NULL.
 ****************************************************************************/
char *SearchForProto(const char *szFuncName)
{
  /* BLAS and LAPACK prototypes generated from LAPACK version 3.1.1.
     The variable pszProtos is an array of strings sorted
     alphabetically by function name, which are prototypes for all of
     the BLAS and LAPACK functions from a given LAPACK source
     distribution. */
  static char *pszProtos[] = {
    "CAXPY(I,C,C,I,C,I)",
    "CBDSQR(H,I,I,I,I,S,S,C,I,C,I,C,I,S,I)",
    "CCOPY(I,C,I,C,I)",
    "C=CDOTC(I,C,I,C,I)",
    "C=CDOTU(I,C,I,C,I)",
    "CGBBRD(H,I,I,I,I,I,C,I,S,S,C,I,C,I,C,I,C,S,I)",
    "CGBCON(H,I,I,I,C,I,I,S,S,C,S,I)",
    "CGBEQU(I,I,I,I,C,I,S,S,S,S,S,I)",
    "CGBMV(H,I,I,I,I,C,C,I,C,I,C,C,I)",
    "CGBRFS(H,I,I,I,I,C,I,C,I,I,C,I,C,I,S,S,C,S,I)",
    "CGBSV(I,I,I,I,C,I,I,C,I,I)",
    "CGBSVX(H,H,I,I,I,I,C,I,C,I,I,H,S,S,C,I,C,I,S,S,S,C,S,I)",
    "CGBTF2(I,I,I,I,C,I,I,I)",
    "CGBTRF(I,I,I,I,C,I,I,I)",
    "CGBTRS(H,I,I,I,I,C,I,I,C,I,I)",
    "CGEBAK(H,H,I,I,I,S,I,C,I,I)",
    "CGEBAL(H,I,C,I,I,I,S,I)",
    "CGEBD2(I,I,C,I,S,S,C,C,C,I)",
    "CGEBRD(I,I,C,I,S,S,C,C,C,I,I)",
    "CGECON(H,I,C,I,S,S,C,S,I)",
    "CGEEQU(I,I,C,I,S,S,S,S,S,I)",
    "CGEES(H,H,X,I,C,I,I,C,C,I,C,I,S,L,I)",
    "CGEESX(H,H,X,H,I,C,I,I,C,C,I,S,S,C,I,S,L,I)",
    "CGEEV(H,H,I,C,I,C,C,I,C,I,C,I,S,I)",
    "CGEEVX(H,H,H,H,I,C,I,C,C,I,C,I,I,I,S,S,S,S,C,I,S,I)",
    "CGEGS(H,H,I,C,I,C,I,C,C,C,I,C,I,C,I,S,I)",
    "CGEGV(H,H,I,C,I,C,I,C,C,C,I,C,I,C,I,S,I)",
    "CGEHD2(I,I,I,C,I,C,C,I)",
    "CGEHRD(I,I,I,C,I,C,C,I,I)",
    "CGELQ2(I,I,C,I,C,C,I)",
    "CGELQF(I,I,C,I,C,C,I,I)",
    "CGELS(H,I,I,I,C,I,C,I,C,I,I)",
    "CGELSD(I,I,I,C,I,C,I,S,S,I,C,I,S,I,I)",
    "CGELSS(I,I,I,C,I,C,I,S,S,I,C,I,S,I)",
    "CGELSX(I,I,I,C,I,C,I,I,S,I,C,S,I)",
    "CGELSY(I,I,I,C,I,C,I,I,S,I,C,I,S,I)",
    "CGEMM(H,H,I,I,I,C,C,I,C,I,C,C,I)",
    "CGEMV(H,I,I,C,C,I,C,I,C,C,I)",
    "CGEQL2(I,I,C,I,C,C,I)",
    "CGEQLF(I,I,C,I,C,C,I,I)",
    "CGEQP3(I,I,C,I,I,C,C,I,S,I)",
    "CGEQPF(I,I,C,I,I,C,C,S,I)",
    "CGEQR2(I,I,C,I,C,C,I)",
    "CGEQRF(I,I,C,I,C,C,I,I)",
    "CGERC(I,I,C,C,I,C,I,C,I)",
    "CGERFS(H,I,I,C,I,C,I,I,C,I,C,I,S,S,C,S,I)",
    "CGERQ2(I,I,C,I,C,C,I)",
    "CGERQF(I,I,C,I,C,C,I,I)",
    "CGERU(I,I,C,C,I,C,I,C,I)",
    "CGESC2(I,C,I,C,I,I,S)",
    "CGESDD(H,I,I,C,I,S,C,I,C,I,C,I,S,I,I)",
    "CGESV(I,I,C,I,I,C,I,I)",
    "CGESVD(H,H,I,I,C,I,S,C,I,C,I,C,I,S,I)",
    "CGESVX(H,H,I,I,C,I,C,I,I,H,S,S,C,I,C,I,S,S,S,C,S,I)",
    "CGETC2(I,C,I,I,I,I)",
    "CGETF2(I,I,C,I,I,I)",
    "CGETRF(I,I,C,I,I,I)",
    "CGETRI(I,C,I,I,C,I,I)",
    "CGETRS(H,I,I,C,I,I,C,I,I)",
    "CGGBAK(H,H,I,I,I,S,S,I,C,I,I)",
    "CGGBAL(H,I,C,I,C,I,I,I,S,S,S,I)",
    "CGGES(H,H,H,X,I,C,I,C,I,I,C,C,C,I,C,I,C,I,S,L,I)",
    "CGGESX(H,H,H,X,H,I,C,I,C,I,I,C,C,C,I,C,I,S,S,C,I,S,I,I,L,I)",
    "CGGEV(H,H,I,C,I,C,I,C,C,C,I,C,I,C,I,S,I)",
    "CGGEVX(H,H,H,H,I,C,I,C,I,C,C,C,I,C,I,I,I,S,S,S,S,S,S,C,I,S,I,L,I)",
    "CGGGLM(I,I,I,C,I,C,I,C,C,C,C,I,I)",
    "CGGHRD(H,H,I,I,I,C,I,C,I,C,I,C,I,I)",
    "CGGLSE(I,I,I,C,I,C,I,C,C,C,C,I,I)",
    "CGGQRF(I,I,I,C,I,C,C,I,C,C,I,I)",
    "CGGRQF(I,I,I,C,I,C,C,I,C,C,I,I)",
    "CGGSVD(H,H,H,I,I,I,I,I,C,I,C,I,S,S,C,I,C,I,C,I,C,S,I,I)",
    "CGGSVP(H,H,H,I,I,I,C,I,C,I,S,S,I,I,C,I,C,I,C,I,I,S,C,C,I)",
    "CGTCON(H,I,C,C,C,C,I,S,S,C,I)",
    "CGTRFS(H,I,I,C,C,C,C,C,C,C,I,C,I,C,I,S,S,C,S,I)",
    "CGTSV(I,I,C,C,C,C,I,I)",
    "CGTSVX(H,H,I,I,C,C,C,C,C,C,C,I,C,I,C,I,S,S,S,C,S,I)",
    "CGTTRF(I,C,C,C,C,I,I)",
    "CGTTRS(H,I,I,C,C,C,C,I,C,I,I)",
    "CGTTS2(I,I,I,C,C,C,C,I,C,I)",
    "CHBEV(H,H,I,I,C,I,S,C,I,C,S,I)",
    "CHBEVD(H,H,I,I,C,I,S,C,I,C,I,S,I,I,I,I)",
    "CHBEVX(H,H,H,I,I,C,I,C,I,S,S,I,I,S,I,S,C,I,C,S,I,I,I)",
    "CHBGST(H,H,I,I,I,C,I,C,I,C,I,C,S,I)",
    "CHBGV(H,H,I,I,I,C,I,C,I,S,C,I,C,S,I)",
    "CHBGVD(H,H,I,I,I,C,I,C,I,S,C,I,C,I,S,I,I,I,I)",
    "CHBGVX(H,H,H,I,I,I,C,I,C,I,C,I,S,S,I,I,S,I,S,C,I,C,S,I,I,I)",
    "CHBMV(H,I,I,C,C,I,C,I,C,C,I)",
    "CHBTRD(H,H,I,I,C,I,S,S,C,I,C,I)",
    "CHECON(H,I,C,I,I,S,S,C,I)",
    "CHEEV(H,H,I,C,I,S,C,I,S,I)",
    "CHEEVD(H,H,I,C,I,S,C,I,S,I,I,I,I)",
    "CHEEVR(H,H,H,I,C,I,S,S,I,I,S,I,S,C,I,I,C,I,S,I,I,I,I)",
    "CHEEVX(H,H,H,I,C,I,S,S,I,I,S,I,S,C,I,C,I,S,I,I,I)",
    "CHEGS2(I,H,I,C,I,C,I,I)",
    "CHEGST(I,H,I,C,I,C,I,I)",
    "CHEGV(I,H,H,I,C,I,C,I,S,C,I,S,I)",
    "CHEGVD(I,H,H,I,C,I,C,I,S,C,I,S,I,I,I,I)",
    "CHEGVX(I,H,H,H,I,C,I,C,I,S,S,I,I,S,I,S,C,I,C,I,S,I,I,I)",
    "CHEMM(H,H,I,I,C,C,I,C,I,C,C,I)",
    "CHEMV(H,I,C,C,I,C,I,C,C,I)",
    "CHER(H,I,S,C,I,C,I)",
    "CHER2(H,I,C,C,I,C,I,C,I)",
    "CHER2K(H,H,I,I,C,C,I,C,I,S,C,I)",
    "CHERFS(H,I,I,C,I,C,I,I,C,I,C,I,S,S,C,S,I)",
    "CHERK(H,H,I,I,S,C,I,S,C,I)",
    "CHESV(H,I,I,C,I,I,C,I,C,I,I)",
    "CHESVX(H,H,I,I,C,I,C,I,I,C,I,C,I,S,S,S,C,I,S,I)",
    "CHETD2(H,I,C,I,S,S,C,I)",
    "CHETF2(H,I,C,I,I,I)",
    "CHETRD(H,I,C,I,S,S,C,C,I,I)",
    "CHETRF(H,I,C,I,I,C,I,I)",
    "CHETRI(H,I,C,I,I,C,I)",
    "CHETRS(H,I,I,C,I,I,C,I,I)",
    "CHGEQZ(H,H,H,I,I,I,C,I,C,I,C,C,C,I,C,I,C,I,S,I)",
    "CHPCON(H,I,C,I,S,S,C,I)",
    "CHPEV(H,H,I,C,S,C,I,C,S,I)",
    "CHPEVD(H,H,I,C,S,C,I,C,I,S,I,I,I,I)",
    "CHPEVX(H,H,H,I,C,S,S,I,I,S,I,S,C,I,C,S,I,I,I)",
    "CHPGST(I,H,I,C,C,I)",
    "CHPGV(I,H,H,I,C,C,S,C,I,C,S,I)",
    "CHPGVD(I,H,H,I,C,C,S,C,I,C,I,S,I,I,I,I)",
    "CHPGVX(I,H,H,H,I,C,C,S,S,I,I,S,I,S,C,I,C,S,I,I,I)",
    "CHPMV(H,I,C,C,C,I,C,C,I)",
    "CHPR(H,I,S,C,I,C)",
    "CHPR2(H,I,C,C,I,C,I,C)",
    "CHPRFS(H,I,I,C,C,I,C,I,C,I,S,S,C,S,I)",
    "CHPSV(H,I,I,C,I,C,I,I)",
    "CHPSVX(H,H,I,I,C,C,I,C,I,C,I,S,S,S,C,S,I)",
    "CHPTRD(H,I,C,S,S,C,I)",
    "CHPTRF(H,I,C,I,I)",
    "CHPTRI(H,I,C,I,C,I)",
    "CHPTRS(H,I,I,C,I,C,I,I)",
    "CHSEIN(H,H,H,L,I,C,I,C,C,I,C,I,I,I,C,S,I,I,I)",
    "CHSEQR(H,H,I,I,I,C,I,C,C,I,C,I,I)",
    "CLABRD(I,I,I,C,I,S,S,C,C,C,I,C,I)",
    "CLACGV(I,C,I)",
    "CLACN2(I,C,C,S,I,I)",
    "CLACON(I,C,C,S,I)",
    "CLACP2(H,I,I,S,I,C,I)",
    "CLACPY(H,I,I,C,I,C,I)",
    "CLACRM(I,I,C,I,S,I,C,I,S)",
    "CLACRT(I,C,I,C,I,C,C)",
    "C=CLADIV(C,C)",
    "CLAED0(I,I,S,S,C,I,C,I,S,I,I)",
    "CLAED7(I,I,I,I,I,I,S,C,I,S,I,S,I,I,I,I,I,S,C,S,I,I)",
    "CLAED8(I,I,I,C,I,S,S,I,S,S,C,I,S,I,I,I,I,I,I,S,I)",
    "CLAEIN(L,L,I,C,I,C,C,C,I,S,S,S,I)",
    "CLAESY(C,C,C,C,C,C,C,C)",
    "CLAEV2(C,C,C,S,S,S,C)",
    "CLAG2Z(I,I,C,I,Z,I,I)",
    "CLAGS2(L,S,C,S,S,C,S,S,C,S,C,S,C)",
    "CLAGTM(H,I,I,S,C,C,C,C,I,S,C,I)",
    "CLAHEF(H,I,I,I,C,I,I,C,I,I)",
    "CLAHQR(L,L,I,I,I,C,I,C,I,I,C,I,I)",
    "CLAHR2(I,I,I,C,I,C,C,I,C,I)",
    "CLAHRD(I,I,I,C,I,C,C,I,C,I)",
    "CLAIC1(I,I,C,S,C,C,S,C,C)",
    "CLALS0(I,I,I,I,I,C,I,C,I,I,I,I,I,S,I,S,S,S,S,I,S,S,S,I)",
    "CLALSA(I,I,I,I,C,I,C,I,S,I,S,I,S,S,S,S,I,I,I,I,S,S,S,S,I,I)",
    "CLALSD(H,I,I,I,S,S,C,I,S,I,C,S,I,I)",
    "S=CLANGB(H,I,I,I,C,I,S)",
    "S=CLANGE(H,I,I,C,I,S)",
    "S=CLANGT(H,I,C,C,C)",
    "S=CLANHB(H,H,I,I,C,I,S)",
    "S=CLANHE(H,H,I,C,I,S)",
    "S=CLANHP(H,H,I,C,S)",
    "S=CLANHS(H,I,C,I,S)",
    "S=CLANHT(H,I,S,C)",
    "S=CLANSB(H,H,I,I,C,I,S)",
    "S=CLANSP(H,H,I,C,S)",
    "S=CLANSY(H,H,I,C,I,S)",
    "S=CLANTB(H,H,H,I,I,C,I,S)",
    "S=CLANTP(H,H,H,I,C,S)",
    "S=CLANTR(H,H,H,I,I,C,I,S)",
    "CLAPLL(I,C,I,C,I,S)",
    "CLAPMT(L,I,I,C,I,I)",
    "CLAQGB(I,I,I,I,C,I,S,S,S,S,S,H)",
    "CLAQGE(I,I,C,I,S,S,S,S,S,H)",
    "CLAQHB(H,I,I,C,I,S,S,S,H)",
    "CLAQHE(H,I,C,I,S,S,S,H)",
    "CLAQHP(H,I,C,S,S,S,H)",
    "CLAQP2(I,I,I,C,I,I,C,S,S,C)",
    "CLAQPS(I,I,I,I,I,C,I,I,C,S,S,C,C,I)",
    "CLAQR0(L,L,I,I,I,C,I,C,I,I,C,I,C,I,I)",
    "CLAQR1(I,C,I,C,C,C)",
    "CLAQR2(L,L,I,I,I,I,C,I,I,I,C,I,I,I,C,C,I,I,C,I,I,C,I,C,I)",
    "CLAQR3(L,L,I,I,I,I,C,I,I,I,C,I,I,I,C,C,I,I,C,I,I,C,I,C,I)",
    "CLAQR4(L,L,I,I,I,C,I,C,I,I,C,I,C,I,I)",
    "CLAQR5(L,L,I,I,I,I,I,C,C,I,I,I,C,I,C,I,C,I,I,C,I,I,C,I)",
    "CLAQSB(H,I,I,C,I,S,S,S,H)",
    "CLAQSP(H,I,C,S,S,S,H)",
    "CLAQSY(H,I,C,I,S,S,S,H)",
    "CLAR1V(I,I,I,S,S,S,S,S,S,S,C,L,I,S,S,I,I,S,S,S,S)",
    "CLAR2V(I,C,C,C,I,S,C,I)",
    "CLARCM(I,I,S,I,C,I,C,I,S)",
    "CLARF(H,I,I,C,I,C,C,I,C)",
    "CLARFB(H,H,H,H,I,I,I,C,I,C,I,C,I,C,I)",
    "CLARFG(I,C,C,I,C)",
    "CLARFT(H,H,I,I,C,I,C,C,I)",
    "CLARFX(H,I,I,C,C,C,I,C)",
    "CLARGV(I,C,I,C,I,S,I)",
    "CLARNV(I,I,I,C)",
    "CLARRV(I,S,S,S,S,S,I,I,I,I,S,S,S,S,S,S,I,I,S,C,I,I,S,I,I)",
    "CLARTG(C,C,S,C,C)",
    "CLARTV(I,C,I,C,I,S,C,I)",
    "CLARZ(H,I,I,I,C,I,C,C,I,C)",
    "CLARZB(H,H,H,H,I,I,I,I,C,I,C,I,C,I,C,I)",
    "CLARZT(H,H,I,I,C,I,C,C,I)",
    "CLASCL(H,I,I,S,S,I,I,C,I,I)",
    "CLASET(H,I,I,C,C,C,I)",
    "CLASR(H,H,H,I,I,S,S,C,I)",
    "CLASSQ(I,C,I,S,S)",
    "CLASWP(I,C,I,I,I,I,I)",
    "CLASYF(H,I,I,I,C,I,I,C,I,I)",
    "CLATBS(H,H,H,H,I,I,C,I,C,S,S,I)",
    "CLATDF(I,I,C,I,C,S,S,I,I)",
    "CLATPS(H,H,H,H,I,C,C,S,S,I)",
    "CLATRD(H,I,I,C,I,S,C,C,I)",
    "CLATRS(H,H,H,H,I,C,I,C,S,S,I)",
    "CLATRZ(I,I,I,C,I,C,C)",
    "CLATZM(H,I,I,C,I,C,C,C,I,C)",
    "CLAUU2(H,I,C,I,I)",
    "CLAUUM(H,I,C,I,I)",
    "CPBCON(H,I,I,C,I,S,S,C,S,I)",
    "CPBEQU(H,I,I,C,I,S,S,S,I)",
    "CPBRFS(H,I,I,I,C,I,C,I,C,I,C,I,S,S,C,S,I)",
    "CPBSTF(H,I,I,C,I,I)",
    "CPBSV(H,I,I,I,C,I,C,I,I)",
    "CPBSVX(H,H,I,I,I,C,I,C,I,H,S,C,I,C,I,S,S,S,C,S,I)",
    "CPBTF2(H,I,I,C,I,I)",
    "CPBTRF(H,I,I,C,I,I)",
    "CPBTRS(H,I,I,I,C,I,C,I,I)",
    "CPOCON(H,I,C,I,S,S,C,S,I)",
    "CPOEQU(I,C,I,S,S,S,I)",
    "CPORFS(H,I,I,C,I,C,I,C,I,C,I,S,S,C,S,I)",
    "CPOSV(H,I,I,C,I,C,I,I)",
    "CPOSVX(H,H,I,I,C,I,C,I,H,S,C,I,C,I,S,S,S,C,S,I)",
    "CPOTF2(H,I,C,I,I)",
    "CPOTRF(H,I,C,I,I)",
    "CPOTRI(H,I,C,I,I)",
    "CPOTRS(H,I,I,C,I,C,I,I)",
    "CPPCON(H,I,C,S,S,C,S,I)",
    "CPPEQU(H,I,C,S,S,S,I)",
    "CPPRFS(H,I,I,C,C,C,I,C,I,S,S,C,S,I)",
    "CPPSV(H,I,I,C,C,I,I)",
    "CPPSVX(H,H,I,I,C,C,H,S,C,I,C,I,S,S,S,C,S,I)",
    "CPPTRF(H,I,C,I)",
    "CPPTRI(H,I,C,I)",
    "CPPTRS(H,I,I,C,C,I,I)",
    "CPTCON(I,S,C,S,S,S,I)",
    "CPTEQR(H,I,S,S,C,I,S,I)",
    "CPTRFS(H,I,I,S,C,S,C,C,I,C,I,S,S,C,S,I)",
    "CPTSV(I,I,S,C,C,I,I)",
    "CPTSVX(H,I,I,S,C,S,C,C,I,C,I,S,S,S,C,S,I)",
    "CPTTRF(I,S,C,I)",
    "CPTTRS(H,I,I,S,C,C,I,I)",
    "CPTTS2(I,I,I,S,C,C,I)",
    "CROT(I,C,I,C,I,S,C)",
    "CROTG(C,C,S,C)",
    "CSCAL(I,C,C,I)",
    "CSPCON(H,I,C,I,S,S,C,I)",
    "CSPMV(H,I,C,C,C,I,C,C,I)",
    "CSPR(H,I,C,C,I,C)",
    "CSPRFS(H,I,I,C,C,I,C,I,C,I,S,S,C,S,I)",
    "CSPSV(H,I,I,C,I,C,I,I)",
    "CSPSVX(H,H,I,I,C,C,I,C,I,C,I,S,S,S,C,S,I)",
    "CSPTRF(H,I,C,I,I)",
    "CSPTRI(H,I,C,I,C,I)",
    "CSPTRS(H,I,I,C,I,C,I,I)",
    "CSROT(I,C,I,C,I,S,S)",
    "CSRSCL(I,S,C,I)",
    "CSSCAL(I,S,C,I)",
    "CSTEDC(H,I,S,S,C,I,C,I,S,I,I,I,I)",
    "CSTEGR(H,H,I,S,S,S,S,I,I,S,I,S,C,I,I,S,I,I,I,I)",
    "CSTEIN(I,S,S,I,S,I,I,C,I,S,I,I,I)",
    "CSTEMR(H,H,I,S,S,S,S,I,I,I,S,C,I,I,I,L,S,I,I,I,I)",
    "CSTEQR(H,I,S,S,C,I,S,I)",
    "CSWAP(I,C,I,C,I)",
    "CSYCON(H,I,C,I,I,S,S,C,I)",
    "CSYMM(H,H,I,I,C,C,I,C,I,C,C,I)",
    "CSYMV(H,I,C,C,I,C,I,C,C,I)",
    "CSYR(H,I,C,C,I,C,I)",
    "CSYR2K(H,H,I,I,C,C,I,C,I,C,C,I)",
    "CSYRFS(H,I,I,C,I,C,I,I,C,I,C,I,S,S,C,S,I)",
    "CSYRK(H,H,I,I,C,C,I,C,C,I)",
    "CSYSV(H,I,I,C,I,I,C,I,C,I,I)",
    "CSYSVX(H,H,I,I,C,I,C,I,I,C,I,C,I,S,S,S,C,I,S,I)",
    "CSYTF2(H,I,C,I,I,I)",
    "CSYTRF(H,I,C,I,I,C,I,I)",
    "CSYTRI(H,I,C,I,I,C,I)",
    "CSYTRS(H,I,I,C,I,I,C,I,I)",
    "CTBCON(H,H,H,I,I,C,I,S,C,S,I)",
    "CTBMV(H,H,H,I,I,C,I,C,I)",
    "CTBRFS(H,H,H,I,I,I,C,I,C,I,C,I,S,S,C,S,I)",
    "CTBSV(H,H,H,I,I,C,I,C,I)",
    "CTBTRS(H,H,H,I,I,I,C,I,C,I,I)",
    "CTGEVC(H,H,L,I,C,I,C,I,C,I,C,I,I,I,C,S,I)",
    "CTGEX2(L,L,I,C,I,C,I,C,I,C,I,I,I)",
    "CTGEXC(L,L,I,C,I,C,I,C,I,C,I,I,I,I)",
    "CTGSEN(I,L,L,L,I,C,I,C,I,C,C,C,I,C,I,I,S,S,S,C,I,I,I,I)",
    "CTGSJA(H,H,H,I,I,I,I,I,C,I,C,I,S,S,S,S,C,I,C,I,C,I,C,I,I)",
    "CTGSNA(H,H,L,I,C,I,C,I,C,I,C,I,S,S,I,I,C,I,I,I)",
    "CTGSY2(H,I,I,I,C,I,C,I,C,I,C,I,C,I,C,I,S,S,S,I)",
    "CTGSYL(H,I,I,I,C,I,C,I,C,I,C,I,C,I,C,I,S,S,C,I,I,I)",
    "CTPCON(H,H,H,I,C,S,C,S,I)",
    "CTPMV(H,H,H,I,C,C,I)",
    "CTPRFS(H,H,H,I,I,C,C,I,C,I,S,S,C,S,I)",
    "CTPSV(H,H,H,I,C,C,I)",
    "CTPTRI(H,H,I,C,I)",
    "CTPTRS(H,H,H,I,I,C,C,I,I)",
    "CTRCON(H,H,H,I,C,I,S,C,S,I)",
    "CTREVC(H,H,L,I,C,I,C,I,C,I,I,I,C,S,I)",
    "CTREXC(H,I,C,I,C,I,I,I,I)",
    "CTRMM(H,H,H,H,I,I,C,C,I,C,I)",
    "CTRMV(H,H,H,I,C,I,C,I)",
    "CTRRFS(H,H,H,I,I,C,I,C,I,C,I,S,S,C,S,I)",
    "CTRSEN(H,H,L,I,C,I,C,I,C,I,S,S,C,I,I)",
    "CTRSM(H,H,H,H,I,I,C,C,I,C,I)",
    "CTRSNA(H,H,L,I,C,I,C,I,C,I,S,S,I,I,C,I,S,I)",
    "CTRSV(H,H,H,I,C,I,C,I)",
    "CTRSYL(H,H,I,I,I,C,I,C,I,C,I,S,I)",
    "CTRTI2(H,H,I,C,I,I)",
    "CTRTRI(H,H,I,C,I,I)",
    "CTRTRS(H,H,H,I,I,C,I,C,I,I)",
    "CTZRQF(I,I,C,I,C,I)",
    "CTZRZF(I,I,C,I,C,C,I,I)",
    "CUNG2L(I,I,I,C,I,C,C,I)",
    "CUNG2R(I,I,I,C,I,C,C,I)",
    "CUNGBR(H,I,I,I,C,I,C,C,I,I)",
    "CUNGHR(I,I,I,C,I,C,C,I,I)",
    "CUNGL2(I,I,I,C,I,C,C,I)",
    "CUNGLQ(I,I,I,C,I,C,C,I,I)",
    "CUNGQL(I,I,I,C,I,C,C,I,I)",
    "CUNGQR(I,I,I,C,I,C,C,I,I)",
    "CUNGR2(I,I,I,C,I,C,C,I)",
    "CUNGRQ(I,I,I,C,I,C,C,I,I)",
    "CUNGTR(H,I,C,I,C,C,I,I)",
    "CUNM2L(H,H,I,I,I,C,I,C,C,I,C,I)",
    "CUNM2R(H,H,I,I,I,C,I,C,C,I,C,I)",
    "CUNMBR(H,H,H,I,I,I,C,I,C,C,I,C,I,I)",
    "CUNMHR(H,H,I,I,I,I,C,I,C,C,I,C,I,I)",
    "CUNML2(H,H,I,I,I,C,I,C,C,I,C,I)",
    "CUNMLQ(H,H,I,I,I,C,I,C,C,I,C,I,I)",
    "CUNMQL(H,H,I,I,I,C,I,C,C,I,C,I,I)",
    "CUNMQR(H,H,I,I,I,C,I,C,C,I,C,I,I)",
    "CUNMR2(H,H,I,I,I,C,I,C,C,I,C,I)",
    "CUNMR3(H,H,I,I,I,I,C,I,C,C,I,C,I)",
    "CUNMRQ(H,H,I,I,I,C,I,C,C,I,C,I,I)",
    "CUNMRZ(H,H,I,I,I,I,C,I,C,C,I,C,I,I)",
    "CUNMTR(H,H,H,I,I,C,I,C,C,I,C,I,I)",
    "CUPGTR(H,I,C,C,C,I,C,I)",
    "CUPMTR(H,H,H,I,I,C,C,C,I,C,I)",
    "D=DASUM(I,D,I)",
    "DAXPY(I,D,D,I,D,I)",
    "DBDSDC(H,H,I,D,D,D,I,D,I,D,I,D,I,I)",
    "DBDSQR(H,I,I,I,I,D,D,D,I,D,I,D,I,D,I)",
    "D=DCABS1(Z)",
    "DCOPY(I,D,I,D,I)",
    "DDISNA(H,I,I,D,D,I)",
    "D=DDOT(I,D,I,D,I)",
    "DGBBRD(H,I,I,I,I,I,D,I,D,D,D,I,D,I,D,I,D,I)",
    "DGBCON(H,I,I,I,D,I,I,D,D,D,I,I)",
    "DGBEQU(I,I,I,I,D,I,D,D,D,D,D,I)",
    "DGBMV(H,I,I,I,I,D,D,I,D,I,D,D,I)",
    "DGBRFS(H,I,I,I,I,D,I,D,I,I,D,I,D,I,D,D,D,I,I)",
    "DGBSV(I,I,I,I,D,I,I,D,I,I)",
    "DGBSVX(H,H,I,I,I,I,D,I,D,I,I,H,D,D,D,I,D,I,D,D,D,D,I,I)",
    "DGBTF2(I,I,I,I,D,I,I,I)",
    "DGBTRF(I,I,I,I,D,I,I,I)",
    "DGBTRS(H,I,I,I,I,D,I,I,D,I,I)",
    "DGEBAK(H,H,I,I,I,D,I,D,I,I)",
    "DGEBAL(H,I,D,I,I,I,D,I)",
    "DGEBD2(I,I,D,I,D,D,D,D,D,I)",
    "DGEBRD(I,I,D,I,D,D,D,D,D,I,I)",
    "DGECON(H,I,D,I,D,D,D,I,I)",
    "DGEEQU(I,I,D,I,D,D,D,D,D,I)",
    "DGEES(H,H,X,I,D,I,I,D,D,D,I,D,I,L,I)",
    "DGEESX(H,H,X,H,I,D,I,I,D,D,D,I,D,D,D,I,I,I,L,I)",
    "DGEEV(H,H,I,D,I,D,D,D,I,D,I,D,I,I)",
    "DGEEVX(H,H,H,H,I,D,I,D,D,D,I,D,I,I,I,D,D,D,D,D,I,I,I)",
    "DGEGS(H,H,I,D,I,D,I,D,D,D,D,I,D,I,D,I,I)",
    "DGEGV(H,H,I,D,I,D,I,D,D,D,D,I,D,I,D,I,I)",
    "DGEHD2(I,I,I,D,I,D,D,I)",
    "DGEHRD(I,I,I,D,I,D,D,I,I)",
    "DGELQ2(I,I,D,I,D,D,I)",
    "DGELQF(I,I,D,I,D,D,I,I)",
    "DGELS(H,I,I,I,D,I,D,I,D,I,I)",
    "DGELSD(I,I,I,D,I,D,I,D,D,I,D,I,I,I)",
    "DGELSS(I,I,I,D,I,D,I,D,D,I,D,I,I)",
    "DGELSX(I,I,I,D,I,D,I,I,D,I,D,I)",
    "DGELSY(I,I,I,D,I,D,I,I,D,I,D,I,I)",
    "DGEMM(H,H,I,I,I,D,D,I,D,I,D,D,I)",
    "DGEMV(H,I,I,D,D,I,D,I,D,D,I)",
    "DGEQL2(I,I,D,I,D,D,I)",
    "DGEQLF(I,I,D,I,D,D,I,I)",
    "DGEQP3(I,I,D,I,I,D,D,I,I)",
    "DGEQPF(I,I,D,I,I,D,D,I)",
    "DGEQR2(I,I,D,I,D,D,I)",
    "DGEQRF(I,I,D,I,D,D,I,I)",
    "DGER(I,I,D,D,I,D,I,D,I)",
    "DGERFS(H,I,I,D,I,D,I,I,D,I,D,I,D,D,D,I,I)",
    "DGERQ2(I,I,D,I,D,D,I)",
    "DGERQF(I,I,D,I,D,D,I,I)",
    "DGESC2(I,D,I,D,I,I,D)",
    "DGESDD(H,I,I,D,I,D,D,I,D,I,D,I,I,I)",
    "DGESV(I,I,D,I,I,D,I,I)",
    "DGESVD(H,H,I,I,D,I,D,D,I,D,I,D,I,I)",
    "DGESVX(H,H,I,I,D,I,D,I,I,H,D,D,D,I,D,I,D,D,D,D,I,I)",
    "DGETC2(I,D,I,I,I,I)",
    "DGETF2(I,I,D,I,I,I)",
    "DGETRF(I,I,D,I,I,I)",
    "DGETRI(I,D,I,I,D,I,I)",
    "DGETRS(H,I,I,D,I,I,D,I,I)",
    "DGGBAK(H,H,I,I,I,D,D,I,D,I,I)",
    "DGGBAL(H,I,D,I,D,I,I,I,D,D,D,I)",
    "DGGES(H,H,H,X,I,D,I,D,I,I,D,D,D,D,I,D,I,D,I,L,I)",
    "DGGESX(H,H,H,X,H,I,D,I,D,I,I,D,D,D,D,I,D,I,D,D,D,I,I,I,L,I)",
    "DGGEV(H,H,I,D,I,D,I,D,D,D,D,I,D,I,D,I,I)",
    "DGGEVX(H,H,H,H,I,D,I,D,I,D,D,D,D,I,D,I,I,I,D,D,D,D,D,D,D,I,I,L,I)",
    "DGGGLM(I,I,I,D,I,D,I,D,D,D,D,I,I)",
    "DGGHRD(H,H,I,I,I,D,I,D,I,D,I,D,I,I)",
    "DGGLSE(I,I,I,D,I,D,I,D,D,D,D,I,I)",
    "DGGQRF(I,I,I,D,I,D,D,I,D,D,I,I)",
    "DGGRQF(I,I,I,D,I,D,D,I,D,D,I,I)",
    "DGGSVD(H,H,H,I,I,I,I,I,D,I,D,I,D,D,D,I,D,I,D,I,D,I,I)",
    "DGGSVP(H,H,H,I,I,I,D,I,D,I,D,D,I,I,D,I,D,I,D,I,I,D,D,I)",
    "DGTCON(H,I,D,D,D,D,I,D,D,D,I,I)",
    "DGTRFS(H,I,I,D,D,D,D,D,D,D,I,D,I,D,I,D,D,D,I,I)",
    "DGTSV(I,I,D,D,D,D,I,I)",
    "DGTSVX(H,H,I,I,D,D,D,D,D,D,D,I,D,I,D,I,D,D,D,D,I,I)",
    "DGTTRF(I,D,D,D,D,I,I)",
    "DGTTRS(H,I,I,D,D,D,D,I,D,I,I)",
    "DGTTS2(I,I,I,D,D,D,D,I,D,I)",
    "DHGEQZ(H,H,H,I,I,I,D,I,D,I,D,D,D,D,I,D,I,D,I,I)",
    "DHSEIN(H,H,H,L,I,D,I,D,D,D,I,D,I,I,I,D,I,I,I)",
    "DHSEQR(H,H,I,I,I,D,I,D,D,D,I,D,I,I)",
    "L=DISNAN(D)",
    "DLABAD(D,D)",
    "DLABRD(I,I,I,D,I,D,D,D,D,D,I,D,I)",
    "DLACN2(I,D,D,I,D,I,I)",
    "DLACON(I,D,D,I,D,I)",
    "DLACPY(H,I,I,D,I,D,I)",
    "DLADIV(D,D,D,D,D,D)",
    "DLAE2(D,D,D,D,D)",
    "DLAEBZ(I,I,I,I,I,I,D,D,D,D,D,D,I,D,D,I,I,D,I,I)",
    "DLAED0(I,I,I,D,D,D,I,D,I,D,I,I)",
    "DLAED1(I,D,D,I,I,D,I,D,I,I)",
    "DLAED2(I,I,I,D,D,I,I,D,D,D,D,D,I,I,I,I,I)",
    "DLAED3(I,I,I,D,D,I,D,D,D,I,I,D,D,I)",
    "DLAED4(I,I,D,D,D,D,D,I)",
    "DLAED5(I,D,D,D,D,D)",
    "DLAED6(I,L,D,D,D,D,D,I)",
    "DLAED7(I,I,I,I,I,I,D,D,I,I,D,I,D,I,I,I,I,I,D,D,I,I)",
    "DLAED8(I,I,I,I,D,D,I,I,D,I,D,D,D,I,D,I,I,I,D,I,I,I)",
    "DLAED9(I,I,I,I,D,D,I,D,D,D,D,I,I)",
    "DLAEDA(I,I,I,I,I,I,I,I,D,D,I,D,D,I)",
    "DLAEIN(L,L,I,D,I,D,D,D,D,D,I,D,D,D,D,I)",
    "DLAEV2(D,D,D,D,D,D,D)",
    "DLAEXC(L,I,D,I,D,I,I,I,I,D,I)",
    "DLAG2(D,I,D,I,D,D,D,D,D,D)",
    "DLAG2S(I,I,D,I,S,I,I)",
    "DLAGS2(L,D,D,D,D,D,D,D,D,D,D,D,D)",
    "DLAGTF(I,D,D,D,D,D,D,I,I)",
    "DLAGTM(H,I,I,D,D,D,D,D,I,D,D,I)",
    "DLAGTS(I,I,D,D,D,D,I,D,D,I)",
    "DLAGV2(D,I,D,I,D,D,D,D,D,D,D)",
    "DLAHQR(L,L,I,I,I,D,I,D,D,I,I,D,I,I)",
    "DLAHR2(I,I,I,D,I,D,D,I,D,I)",
    "DLAHRD(I,I,I,D,I,D,D,I,D,I)",
    "DLAIC1(I,I,D,D,D,D,D,D,D)",
    "L=DLAISNAN(D,D)",
    "DLALN2(L,I,I,D,D,D,I,D,D,D,I,D,D,D,I,D,D,I)",
    "DLALS0(I,I,I,I,I,D,I,D,I,I,I,I,I,D,I,D,D,D,D,I,D,D,D,I)",
    "DLALSA(I,I,I,I,D,I,D,I,D,I,D,I,D,D,D,D,I,I,I,I,D,D,D,D,I,I)",
    "DLALSD(H,I,I,I,D,D,D,I,D,I,D,I,I)",
    "DLAMRG(I,I,D,I,I,I)",
    "I=DLANEG(I,D,D,D,D,I)",
    "D=DLANGB(H,I,I,I,D,I,D)",
    "D=DLANGE(H,I,I,D,I,D)",
    "D=DLANGT(H,I,D,D,D)",
    "D=DLANHS(H,I,D,I,D)",
    "D=DLANSB(H,H,I,I,D,I,D)",
    "D=DLANSP(H,H,I,D,D)",
    "D=DLANST(H,I,D,D)",
    "D=DLANSY(H,H,I,D,I,D)",
    "D=DLANTB(H,H,H,I,I,D,I,D)",
    "D=DLANTP(H,H,H,I,D,D)",
    "D=DLANTR(H,H,H,I,I,D,I,D)",
    "DLANV2(D,D,D,D,D,D,D,D,D,D)",
    "DLAPLL(I,D,I,D,I,D)",
    "DLAPMT(L,I,I,D,I,I)",
    "D=DLAPY2(D,D)",
    "D=DLAPY3(D,D,D)",
    "DLAQGB(I,I,I,I,D,I,D,D,D,D,D,H)",
    "DLAQGE(I,I,D,I,D,D,D,D,D,H)",
    "DLAQP2(I,I,I,D,I,I,D,D,D,D)",
    "DLAQPS(I,I,I,I,I,D,I,I,D,D,D,D,D,I)",
    "DLAQR0(L,L,I,I,I,D,I,D,D,I,I,D,I,D,I,I)",
    "DLAQR1(I,D,I,D,D,D,D,D)",
    "DLAQR2(L,L,I,I,I,I,D,I,I,I,D,I,I,I,D,D,D,I,I,D,I,I,D,I,D,I)",
    "DLAQR3(L,L,I,I,I,I,D,I,I,I,D,I,I,I,D,D,D,I,I,D,I,I,D,I,D,I)",
    "DLAQR4(L,L,I,I,I,D,I,D,D,I,I,D,I,D,I,I)",
    "DLAQR5(L,L,I,I,I,I,I,D,D,D,I,I,I,D,I,D,I,D,I,I,D,I,I,D,I)",
    "DLAQSB(H,I,I,D,I,D,D,D,H)",
    "DLAQSP(H,I,D,D,D,D,H)",
    "DLAQSY(H,I,D,I,D,D,D,H)",
    "DLAQTR(L,L,I,D,I,D,D,D,D,D,I)",
    "DLAR1V(I,I,I,D,D,D,D,D,D,D,D,L,I,D,D,I,I,D,D,D,D)",
    "DLAR2V(I,D,D,D,I,D,D,I)",
    "DLARF(H,I,I,D,I,D,D,I,D)",
    "DLARFB(H,H,H,H,I,I,I,D,I,D,I,D,I,D,I)",
    "DLARFG(I,D,D,I,D)",
    "DLARFT(H,H,I,I,D,I,D,D,I)",
    "DLARFX(H,I,I,D,D,D,I,D)",
    "DLARGV(I,D,I,D,I,D,I)",
    "DLARNV(I,I,I,D)",
    "DLARRA(I,D,D,D,D,D,I,I,I)",
    "DLARRB(I,D,D,I,I,D,D,I,D,D,D,D,I,D,D,I,I)",
    "DLARRC(H,I,D,D,D,D,D,I,I,I,I)",
    "DLARRD(H,H,I,D,D,I,I,D,D,D,D,D,D,I,I,I,D,D,D,D,I,I,D,I,I)",
    "DLARRE(H,I,D,D,I,I,D,D,D,D,D,D,I,I,I,D,D,D,I,I,D,D,D,I,I)",
    "DLARRF(I,D,D,D,I,I,D,D,D,D,D,D,D,D,D,D,D,I)",
    "DLARRJ(I,D,D,I,I,D,I,D,D,D,I,D,D,I)",
    "DLARRK(I,I,D,D,D,D,D,D,D,D,I)",
    "DLARRR(I,D,D,I)",
    "DLARRV(I,D,D,D,D,D,I,I,I,I,D,D,D,D,D,D,I,I,D,D,I,I,D,I,I)",
    "DLARTG(D,D,D,D,D)",
    "DLARTV(I,D,I,D,I,D,D,I)",
    "DLARUV(I,I,D)",
    "DLARZ(H,I,I,I,D,I,D,D,I,D)",
    "DLARZB(H,H,H,H,I,I,I,I,D,I,D,I,D,I,D,I)",
    "DLARZT(H,H,I,I,D,I,D,D,I)",
    "DLAS2(D,D,D,D,D)",
    "DLASCL(H,I,I,D,D,I,I,D,I,I)",
    "DLASD0(I,I,D,D,D,I,D,I,I,I,D,I)",
    "DLASD1(I,I,I,D,D,D,D,I,D,I,I,I,D,I)",
    "DLASD2(I,I,I,I,D,D,D,D,D,I,D,I,D,D,I,D,I,I,I,I,I,I,I)",
    "DLASD3(I,I,I,I,D,D,I,D,D,I,D,I,D,I,D,I,I,I,D,I)",
    "DLASD4(I,I,D,D,D,D,D,D,I)",
    "DLASD5(I,D,D,D,D,D,D)",
    "DLASD6(I,I,I,I,D,D,D,D,D,I,I,I,I,I,D,I,D,D,D,D,I,D,D,D,I,I)",
    "DLASD7(I,I,I,I,I,D,D,D,D,D,D,D,D,D,D,I,I,I,I,I,I,I,D,I,D,D,I)",
    "DLASD8(I,I,D,D,D,D,D,D,I,D,D,I)",
    "DLASDA(I,I,I,I,D,D,D,I,D,I,D,D,D,D,I,I,I,I,D,D,D,D,I,I)",
    "DLASDQ(H,I,I,I,I,I,D,D,D,I,D,I,D,I,D,I)",
    "DLASDT(I,I,I,I,I,I,I)",
    "DLASET(H,I,I,D,D,D,I)",
    "DLASQ1(I,D,D,D,I)",
    "DLASQ2(I,D,I)",
    "DLASQ3(I,I,D,I,D,D,D,D,I,I,I,L)",
    "DLASQ4(I,I,D,I,I,D,D,D,D,D,D,D,I)",
    "DLASQ5(I,I,D,I,D,D,D,D,D,D,D,L)",
    "DLASQ6(I,I,D,I,D,D,D,D,D,D)",
    "DLASR(H,H,H,I,I,D,D,D,I)",
    "DLASRT(H,I,D,I)",
    "DLASSQ(I,D,I,D,D)",
    "DLASV2(D,D,D,D,D,D,D,D,D)",
    "DLASWP(I,D,I,I,I,I,I)",
    "DLASY2(L,L,I,I,I,D,I,D,I,D,I,D,D,I,D,I)",
    "DLASYF(H,I,I,I,D,I,I,D,I,I)",
    "DLATBS(H,H,H,H,I,I,D,I,D,D,D,I)",
    "DLATDF(I,I,D,I,D,D,D,I,I)",
    "DLATPS(H,H,H,H,I,D,D,D,D,I)",
    "DLATRD(H,I,I,D,I,D,D,D,I)",
    "DLATRS(H,H,H,H,I,D,I,D,D,D,I)",
    "DLATRZ(I,I,I,D,I,D,D)",
    "DLATZM(H,I,I,D,I,D,D,D,I,D)",
    "DLAUU2(H,I,D,I,I)",
    "DLAUUM(H,I,D,I,I)",
    "DLAZQ3(I,I,D,I,D,D,D,D,I,I,I,L,I,D,D,D,D,D,D)",
    "DLAZQ4(I,I,D,I,I,D,D,D,D,D,D,D,I,D)",
    "D=DNRM2(I,D,I)",
    "DOPGTR(H,I,D,D,D,I,D,I)",
    "DOPMTR(H,H,H,I,I,D,D,D,I,D,I)",
    "DORG2L(I,I,I,D,I,D,D,I)",
    "DORG2R(I,I,I,D,I,D,D,I)",
    "DORGBR(H,I,I,I,D,I,D,D,I,I)",
    "DORGHR(I,I,I,D,I,D,D,I,I)",
    "DORGL2(I,I,I,D,I,D,D,I)",
    "DORGLQ(I,I,I,D,I,D,D,I,I)",
    "DORGQL(I,I,I,D,I,D,D,I,I)",
    "DORGQR(I,I,I,D,I,D,D,I,I)",
    "DORGR2(I,I,I,D,I,D,D,I)",
    "DORGRQ(I,I,I,D,I,D,D,I,I)",
    "DORGTR(H,I,D,I,D,D,I,I)",
    "DORM2L(H,H,I,I,I,D,I,D,D,I,D,I)",
    "DORM2R(H,H,I,I,I,D,I,D,D,I,D,I)",
    "DORMBR(H,H,H,I,I,I,D,I,D,D,I,D,I,I)",
    "DORMHR(H,H,I,I,I,I,D,I,D,D,I,D,I,I)",
    "DORML2(H,H,I,I,I,D,I,D,D,I,D,I)",
    "DORMLQ(H,H,I,I,I,D,I,D,D,I,D,I,I)",
    "DORMQL(H,H,I,I,I,D,I,D,D,I,D,I,I)",
    "DORMQR(H,H,I,I,I,D,I,D,D,I,D,I,I)",
    "DORMR2(H,H,I,I,I,D,I,D,D,I,D,I)",
    "DORMR3(H,H,I,I,I,I,D,I,D,D,I,D,I)",
    "DORMRQ(H,H,I,I,I,D,I,D,D,I,D,I,I)",
    "DORMRZ(H,H,I,I,I,I,D,I,D,D,I,D,I,I)",
    "DORMTR(H,H,H,I,I,D,I,D,D,I,D,I,I)",
    "DPBCON(H,I,I,D,I,D,D,D,I,I)",
    "DPBEQU(H,I,I,D,I,D,D,D,I)",
    "DPBRFS(H,I,I,I,D,I,D,I,D,I,D,I,D,D,D,I,I)",
    "DPBSTF(H,I,I,D,I,I)",
    "DPBSV(H,I,I,I,D,I,D,I,I)",
    "DPBSVX(H,H,I,I,I,D,I,D,I,H,D,D,I,D,I,D,D,D,D,I,I)",
    "DPBTF2(H,I,I,D,I,I)",
    "DPBTRF(H,I,I,D,I,I)",
    "DPBTRS(H,I,I,I,D,I,D,I,I)",
    "DPOCON(H,I,D,I,D,D,D,I,I)",
    "DPOEQU(I,D,I,D,D,D,I)",
    "DPORFS(H,I,I,D,I,D,I,D,I,D,I,D,D,D,I,I)",
    "DPOSV(H,I,I,D,I,D,I,I)",
    "DPOSVX(H,H,I,I,D,I,D,I,H,D,D,I,D,I,D,D,D,D,I,I)",
    "DPOTF2(H,I,D,I,I)",
    "DPOTRF(H,I,D,I,I)",
    "DPOTRI(H,I,D,I,I)",
    "DPOTRS(H,I,I,D,I,D,I,I)",
    "DPPCON(H,I,D,D,D,D,I,I)",
    "DPPEQU(H,I,D,D,D,D,I)",
    "DPPRFS(H,I,I,D,D,D,I,D,I,D,D,D,I,I)",
    "DPPSV(H,I,I,D,D,I,I)",
    "DPPSVX(H,H,I,I,D,D,H,D,D,I,D,I,D,D,D,D,I,I)",
    "DPPTRF(H,I,D,I)",
    "DPPTRI(H,I,D,I)",
    "DPPTRS(H,I,I,D,D,I,I)",
    "DPTCON(I,D,D,D,D,D,I)",
    "DPTEQR(H,I,D,D,D,I,D,I)",
    "DPTRFS(I,I,D,D,D,D,D,I,D,I,D,D,D,I)",
    "DPTSV(I,I,D,D,D,I,I)",
    "DPTSVX(H,I,I,D,D,D,D,D,I,D,I,D,D,D,D,I)",
    "DPTTRF(I,D,D,I)",
    "DPTTRS(I,I,D,D,D,I,I)",
    "DPTTS2(I,I,D,D,D,I)",
    "DROT(I,D,I,D,I,D,D)",
    "DROTG(D,D,D,D)",
    "DROTM(I,D,I,D,I,D)",
    "DROTMG(D,D,D,D,D)",
    "DRSCL(I,D,D,I)",
    "DSBEV(H,H,I,I,D,I,D,D,I,D,I)",
    "DSBEVD(H,H,I,I,D,I,D,D,I,D,I,I,I,I)",
    "DSBEVX(H,H,H,I,I,D,I,D,I,D,D,I,I,D,I,D,D,I,D,I,I,I)",
    "DSBGST(H,H,I,I,I,D,I,D,I,D,I,D,I)",
    "DSBGV(H,H,I,I,I,D,I,D,I,D,D,I,D,I)",
    "DSBGVD(H,H,I,I,I,D,I,D,I,D,D,I,D,I,I,I,I)",
    "DSBGVX(H,H,H,I,I,I,D,I,D,I,D,I,D,D,I,I,D,I,D,D,I,D,I,I,I)",
    "DSBMV(H,I,I,D,D,I,D,I,D,D,I)",
    "DSBTRD(H,H,I,I,D,I,D,D,D,I,D,I)",
    "DSCAL(I,D,D,I)",
    "D=DSDOT(I,S,I,S,I)",
    "DSGESV(I,I,D,I,I,D,I,D,I,D,S,I,I)",
    "DSPCON(H,I,D,I,D,D,D,I,I)",
    "DSPEV(H,H,I,D,D,D,I,D,I)",
    "DSPEVD(H,H,I,D,D,D,I,D,I,I,I,I)",
    "DSPEVX(H,H,H,I,D,D,D,I,I,D,I,D,D,I,D,I,I,I)",
    "DSPGST(I,H,I,D,D,I)",
    "DSPGV(I,H,H,I,D,D,D,D,I,D,I)",
    "DSPGVD(I,H,H,I,D,D,D,D,I,D,I,I,I,I)",
    "DSPGVX(I,H,H,H,I,D,D,D,D,I,I,D,I,D,D,I,D,I,I,I)",
    "DSPMV(H,I,D,D,D,I,D,D,I)",
    "DSPR(H,I,D,D,I,D)",
    "DSPR2(H,I,D,D,I,D,I,D)",
    "DSPRFS(H,I,I,D,D,I,D,I,D,I,D,D,D,I,I)",
    "DSPSV(H,I,I,D,I,D,I,I)",
    "DSPSVX(H,H,I,I,D,D,I,D,I,D,I,D,D,D,D,I,I)",
    "DSPTRD(H,I,D,D,D,D,I)",
    "DSPTRF(H,I,D,I,I)",
    "DSPTRI(H,I,D,I,D,I)",
    "DSPTRS(H,I,I,D,I,D,I,I)",
    "DSTEBZ(H,H,I,D,D,I,I,D,D,D,I,I,D,I,I,D,I,I)",
    "DSTEDC(H,I,D,D,D,I,D,I,I,I,I)",
    "DSTEGR(H,H,I,D,D,D,D,I,I,D,I,D,D,I,I,D,I,I,I,I)",
    "DSTEIN(I,D,D,I,D,I,I,D,I,D,I,I,I)",
    "DSTEMR(H,H,I,D,D,D,D,I,I,I,D,D,I,I,I,L,D,I,I,I,I)",
    "DSTEQR(H,I,D,D,D,I,D,I)",
    "DSTERF(I,D,D,I)",
    "DSTEV(H,I,D,D,D,I,D,I)",
    "DSTEVD(H,I,D,D,D,I,D,I,I,I,I)",
    "DSTEVR(H,H,I,D,D,D,D,I,I,D,I,D,D,I,I,D,I,I,I,I)",
    "DSTEVX(H,H,I,D,D,D,D,I,I,D,I,D,D,I,D,I,I,I)",
    "DSWAP(I,D,I,D,I)",
    "DSYCON(H,I,D,I,I,D,D,D,I,I)",
    "DSYEV(H,H,I,D,I,D,D,I,I)",
    "DSYEVD(H,H,I,D,I,D,D,I,I,I,I)",
    "DSYEVR(H,H,H,I,D,I,D,D,I,I,D,I,D,D,I,I,D,I,I,I,I)",
    "DSYEVX(H,H,H,I,D,I,D,D,I,I,D,I,D,D,I,D,I,I,I,I)",
    "DSYGS2(I,H,I,D,I,D,I,I)",
    "DSYGST(I,H,I,D,I,D,I,I)",
    "DSYGV(I,H,H,I,D,I,D,I,D,D,I,I)",
    "DSYGVD(I,H,H,I,D,I,D,I,D,D,I,I,I,I)",
    "DSYGVX(I,H,H,H,I,D,I,D,I,D,D,I,I,D,I,D,D,I,D,I,I,I,I)",
    "DSYMM(H,H,I,I,D,D,I,D,I,D,D,I)",
    "DSYMV(H,I,D,D,I,D,I,D,D,I)",
    "DSYR(H,I,D,D,I,D,I)",
    "DSYR2(H,I,D,D,I,D,I,D,I)",
    "DSYR2K(H,H,I,I,D,D,I,D,I,D,D,I)",
    "DSYRFS(H,I,I,D,I,D,I,I,D,I,D,I,D,D,D,I,I)",
    "DSYRK(H,H,I,I,D,D,I,D,D,I)",
    "DSYSV(H,I,I,D,I,I,D,I,D,I,I)",
    "DSYSVX(H,H,I,I,D,I,D,I,I,D,I,D,I,D,D,D,D,I,I,I)",
    "DSYTD2(H,I,D,I,D,D,D,I)",
    "DSYTF2(H,I,D,I,I,I)",
    "DSYTRD(H,I,D,I,D,D,D,D,I,I)",
    "DSYTRF(H,I,D,I,I,D,I,I)",
    "DSYTRI(H,I,D,I,I,D,I)",
    "DSYTRS(H,I,I,D,I,I,D,I,I)",
    "DTBCON(H,H,H,I,I,D,I,D,D,I,I)",
    "DTBMV(H,H,H,I,I,D,I,D,I)",
    "DTBRFS(H,H,H,I,I,I,D,I,D,I,D,I,D,D,D,I,I)",
    "DTBSV(H,H,H,I,I,D,I,D,I)",
    "DTBTRS(H,H,H,I,I,I,D,I,D,I,I)",
    "DTGEVC(H,H,L,I,D,I,D,I,D,I,D,I,I,I,D,I)",
    "DTGEX2(L,L,I,D,I,D,I,D,I,D,I,I,I,I,D,I,I)",
    "DTGEXC(L,L,I,D,I,D,I,D,I,D,I,I,I,D,I,I)",
    "DTGSEN(I,L,L,L,I,D,I,D,I,D,D,D,D,I,D,I,I,D,D,D,D,I,I,I,I)",
    "DTGSJA(H,H,H,I,I,I,I,I,D,I,D,I,D,D,D,D,D,I,D,I,D,I,D,I,I)",
    "DTGSNA(H,H,L,I,D,I,D,I,D,I,D,I,D,D,I,I,D,I,I,I)",
    "DTGSY2(H,I,I,I,D,I,D,I,D,I,D,I,D,I,D,I,D,D,D,I,I,I)",
    "DTGSYL(H,I,I,I,D,I,D,I,D,I,D,I,D,I,D,I,D,D,D,I,I,I)",
    "DTPCON(H,H,H,I,D,D,D,I,I)",
    "DTPMV(H,H,H,I,D,D,I)",
    "DTPRFS(H,H,H,I,I,D,D,I,D,I,D,D,D,I,I)",
    "DTPSV(H,H,H,I,D,D,I)",
    "DTPTRI(H,H,I,D,I)",
    "DTPTRS(H,H,H,I,I,D,D,I,I)",
    "DTRCON(H,H,H,I,D,I,D,D,I,I)",
    "DTREVC(H,H,L,I,D,I,D,I,D,I,I,I,D,I)",
    "DTREXC(H,I,D,I,D,I,I,I,D,I)",
    "DTRMM(H,H,H,H,I,I,D,D,I,D,I)",
    "DTRMV(H,H,H,I,D,I,D,I)",
    "DTRRFS(H,H,H,I,I,D,I,D,I,D,I,D,D,D,I,I)",
    "DTRSEN(H,H,L,I,D,I,D,I,D,D,I,D,D,D,I,I,I,I)",
    "DTRSM(H,H,H,H,I,I,D,D,I,D,I)",
    "DTRSNA(H,H,L,I,D,I,D,I,D,I,D,D,I,I,D,I,I,I)",
    "DTRSV(H,H,H,I,D,I,D,I)",
    "DTRSYL(H,H,I,I,I,D,I,D,I,D,I,D,I)",
    "DTRTI2(H,H,I,D,I,I)",
    "DTRTRI(H,H,I,D,I,I)",
    "DTRTRS(H,H,H,I,I,D,I,D,I,I)",
    "DTZRQF(I,I,D,I,D,I)",
    "DTZRZF(I,I,D,I,D,D,I,I)",
    "D=DZASUM(I,Z,I)",
    "D=DZNRM2(I,Z,I)",
    "D=DZSUM1(I,Z,I)",
    "I=ICAMAX(I,C,I)",
    "I=ICMAX1(I,C,I)",
    "I=IDAMAX(I,D,I)",
    "I=IEEECK(I,S,S)",
    "I=ILAENV(I,H,H,I,I,I,I)",
    "ILAVER(I,I,I)",
    "I=IPARMQ(I,H,H,I,I,I,I)",
    "I=ISAMAX(I,S,I)",
    "I=IZAMAX(I,Z,I)",
    "I=IZMAX1(I,Z,I)",
    "L=LSAME(H,H)",
    "L=LSAMEN(I,H,H)",
    "S=SASUM(I,S,I)",
    "SAXPY(I,S,S,I,S,I)",
    "SBDSDC(H,H,I,S,S,S,I,S,I,S,I,S,I,I)",
    "SBDSQR(H,I,I,I,I,S,S,S,I,S,I,S,I,S,I)",
    "S=SCABS1(C)",
    "S=SCASUM(I,C,I)",
    "S=SCNRM2(I,C,I)",
    "SCOPY(I,S,I,S,I)",
    "S=SCSUM1(I,C,I)",
    "SDISNA(H,I,I,S,S,I)",
    "S=SDOT(I,S,I,S,I)",
    "S=SDSDOT(I,S,S,I,S,I)",
    "SGBBRD(H,I,I,I,I,I,S,I,S,S,S,I,S,I,S,I,S,I)",
    "SGBCON(H,I,I,I,S,I,I,S,S,S,I,I)",
    "SGBEQU(I,I,I,I,S,I,S,S,S,S,S,I)",
    "SGBMV(H,I,I,I,I,S,S,I,S,I,S,S,I)",
    "SGBRFS(H,I,I,I,I,S,I,S,I,I,S,I,S,I,S,S,S,I,I)",
    "SGBSV(I,I,I,I,S,I,I,S,I,I)",
    "SGBSVX(H,H,I,I,I,I,S,I,S,I,I,H,S,S,S,I,S,I,S,S,S,S,I,I)",
    "SGBTF2(I,I,I,I,S,I,I,I)",
    "SGBTRF(I,I,I,I,S,I,I,I)",
    "SGBTRS(H,I,I,I,I,S,I,I,S,I,I)",
    "SGEBAK(H,H,I,I,I,S,I,S,I,I)",
    "SGEBAL(H,I,S,I,I,I,S,I)",
    "SGEBD2(I,I,S,I,S,S,S,S,S,I)",
    "SGEBRD(I,I,S,I,S,S,S,S,S,I,I)",
    "SGECON(H,I,S,I,S,S,S,I,I)",
    "SGEEQU(I,I,S,I,S,S,S,S,S,I)",
    "SGEES(H,H,X,I,S,I,I,S,S,S,I,S,I,L,I)",
    "SGEESX(H,H,X,H,I,S,I,I,S,S,S,I,S,S,S,I,I,I,L,I)",
    "SGEEV(H,H,I,S,I,S,S,S,I,S,I,S,I,I)",
    "SGEEVX(H,H,H,H,I,S,I,S,S,S,I,S,I,I,I,S,S,S,S,S,I,I,I)",
    "SGEGS(H,H,I,S,I,S,I,S,S,S,S,I,S,I,S,I,I)",
    "SGEGV(H,H,I,S,I,S,I,S,S,S,S,I,S,I,S,I,I)",
    "SGEHD2(I,I,I,S,I,S,S,I)",
    "SGEHRD(I,I,I,S,I,S,S,I,I)",
    "SGELQ2(I,I,S,I,S,S,I)",
    "SGELQF(I,I,S,I,S,S,I,I)",
    "SGELS(H,I,I,I,S,I,S,I,S,I,I)",
    "SGELSD(I,I,I,S,I,S,I,S,S,I,S,I,I,I)",
    "SGELSS(I,I,I,S,I,S,I,S,S,I,S,I,I)",
    "SGELSX(I,I,I,S,I,S,I,I,S,I,S,I)",
    "SGELSY(I,I,I,S,I,S,I,I,S,I,S,I,I)",
    "SGEMM(H,H,I,I,I,S,S,I,S,I,S,S,I)",
    "SGEMV(H,I,I,S,S,I,S,I,S,S,I)",
    "SGEQL2(I,I,S,I,S,S,I)",
    "SGEQLF(I,I,S,I,S,S,I,I)",
    "SGEQP3(I,I,S,I,I,S,S,I,I)",
    "SGEQPF(I,I,S,I,I,S,S,I)",
    "SGEQR2(I,I,S,I,S,S,I)",
    "SGEQRF(I,I,S,I,S,S,I,I)",
    "SGER(I,I,S,S,I,S,I,S,I)",
    "SGERFS(H,I,I,S,I,S,I,I,S,I,S,I,S,S,S,I,I)",
    "SGERQ2(I,I,S,I,S,S,I)",
    "SGERQF(I,I,S,I,S,S,I,I)",
    "SGESC2(I,S,I,S,I,I,S)",
    "SGESDD(H,I,I,S,I,S,S,I,S,I,S,I,I,I)",
    "SGESV(I,I,S,I,I,S,I,I)",
    "SGESVD(H,H,I,I,S,I,S,S,I,S,I,S,I,I)",
    "SGESVX(H,H,I,I,S,I,S,I,I,H,S,S,S,I,S,I,S,S,S,S,I,I)",
    "SGETC2(I,S,I,I,I,I)",
    "SGETF2(I,I,S,I,I,I)",
    "SGETRF(I,I,S,I,I,I)",
    "SGETRI(I,S,I,I,S,I,I)",
    "SGETRS(H,I,I,S,I,I,S,I,I)",
    "SGGBAK(H,H,I,I,I,S,S,I,S,I,I)",
    "SGGBAL(H,I,S,I,S,I,I,I,S,S,S,I)",
    "SGGES(H,H,H,X,I,S,I,S,I,I,S,S,S,S,I,S,I,S,I,L,I)",
    "SGGESX(H,H,H,X,H,I,S,I,S,I,I,S,S,S,S,I,S,I,S,S,S,I,I,I,L,I)",
    "SGGEV(H,H,I,S,I,S,I,S,S,S,S,I,S,I,S,I,I)",
    "SGGEVX(H,H,H,H,I,S,I,S,I,S,S,S,S,I,S,I,I,I,S,S,S,S,S,S,S,I,I,L,I)",
    "SGGGLM(I,I,I,S,I,S,I,S,S,S,S,I,I)",
    "SGGHRD(H,H,I,I,I,S,I,S,I,S,I,S,I,I)",
    "SGGLSE(I,I,I,S,I,S,I,S,S,S,S,I,I)",
    "SGGQRF(I,I,I,S,I,S,S,I,S,S,I,I)",
    "SGGRQF(I,I,I,S,I,S,S,I,S,S,I,I)",
    "SGGSVD(H,H,H,I,I,I,I,I,S,I,S,I,S,S,S,I,S,I,S,I,S,I,I)",
    "SGGSVP(H,H,H,I,I,I,S,I,S,I,S,S,I,I,S,I,S,I,S,I,I,S,S,I)",
    "SGTCON(H,I,S,S,S,S,I,S,S,S,I,I)",
    "SGTRFS(H,I,I,S,S,S,S,S,S,S,I,S,I,S,I,S,S,S,I,I)",
    "SGTSV(I,I,S,S,S,S,I,I)",
    "SGTSVX(H,H,I,I,S,S,S,S,S,S,S,I,S,I,S,I,S,S,S,S,I,I)",
    "SGTTRF(I,S,S,S,S,I,I)",
    "SGTTRS(H,I,I,S,S,S,S,I,S,I,I)",
    "SGTTS2(I,I,I,S,S,S,S,I,S,I)",
    "SHGEQZ(H,H,H,I,I,I,S,I,S,I,S,S,S,S,I,S,I,S,I,I)",
    "SHSEIN(H,H,H,L,I,S,I,S,S,S,I,S,I,I,I,S,I,I,I)",
    "SHSEQR(H,H,I,I,I,S,I,S,S,S,I,S,I,I)",
    "L=SISNAN(S)",
    "SLABAD(S,S)",
    "SLABRD(I,I,I,S,I,S,S,S,S,S,I,S,I)",
    "SLACN2(I,S,S,I,S,I,I)",
    "SLACON(I,S,S,I,S,I)",
    "SLACPY(H,I,I,S,I,S,I)",
    "SLADIV(S,S,S,S,S,S)",
    "SLAE2(S,S,S,S,S)",
    "SLAEBZ(I,I,I,I,I,I,S,S,S,S,S,S,I,S,S,I,I,S,I,I)",
    "SLAED0(I,I,I,S,S,S,I,S,I,S,I,I)",
    "SLAED1(I,S,S,I,I,S,I,S,I,I)",
    "SLAED2(I,I,I,S,S,I,I,S,S,S,S,S,I,I,I,I,I)",
    "SLAED3(I,I,I,S,S,I,S,S,S,I,I,S,S,I)",
    "SLAED4(I,I,S,S,S,S,S,I)",
    "SLAED5(I,S,S,S,S,S)",
    "SLAED6(I,L,S,S,S,S,S,I)",
    "SLAED7(I,I,I,I,I,I,S,S,I,I,S,I,S,I,I,I,I,I,S,S,I,I)",
    "SLAED8(I,I,I,I,S,S,I,I,S,I,S,S,S,I,S,I,I,I,S,I,I,I)",
    "SLAED9(I,I,I,I,S,S,I,S,S,S,S,I,I)",
    "SLAEDA(I,I,I,I,I,I,I,I,S,S,I,S,S,I)",
    "SLAEIN(L,L,I,S,I,S,S,S,S,S,I,S,S,S,S,I)",
    "SLAEV2(S,S,S,S,S,S,S)",
    "SLAEXC(L,I,S,I,S,I,I,I,I,S,I)",
    "SLAG2(S,I,S,I,S,S,S,S,S,S)",
    "SLAG2D(I,I,S,I,D,I,I)",
    "SLAGS2(L,S,S,S,S,S,S,S,S,S,S,S,S)",
    "SLAGTF(I,S,S,S,S,S,S,I,I)",
    "SLAGTM(H,I,I,S,S,S,S,S,I,S,S,I)",
    "SLAGTS(I,I,S,S,S,S,I,S,S,I)",
    "SLAGV2(S,I,S,I,S,S,S,S,S,S,S)",
    "SLAHQR(L,L,I,I,I,S,I,S,S,I,I,S,I,I)",
    "SLAHR2(I,I,I,S,I,S,S,I,S,I)",
    "SLAHRD(I,I,I,S,I,S,S,I,S,I)",
    "SLAIC1(I,I,S,S,S,S,S,S,S)",
    "L=SLAISNAN(S,S)",
    "SLALN2(L,I,I,S,S,S,I,S,S,S,I,S,S,S,I,S,S,I)",
    "SLALS0(I,I,I,I,I,S,I,S,I,I,I,I,I,S,I,S,S,S,S,I,S,S,S,I)",
    "SLALSA(I,I,I,I,S,I,S,I,S,I,S,I,S,S,S,S,I,I,I,I,S,S,S,S,I,I)",
    "SLALSD(H,I,I,I,S,S,S,I,S,I,S,I,I)",
    "SLAMRG(I,I,S,I,I,I)",
    "I=SLANEG(I,S,S,S,S,I)",
    "S=SLANGB(H,I,I,I,S,I,S)",
    "S=SLANGE(H,I,I,S,I,S)",
    "S=SLANGT(H,I,S,S,S)",
    "S=SLANHS(H,I,S,I,S)",
    "S=SLANSB(H,H,I,I,S,I,S)",
    "S=SLANSP(H,H,I,S,S)",
    "S=SLANST(H,I,S,S)",
    "S=SLANSY(H,H,I,S,I,S)",
    "S=SLANTB(H,H,H,I,I,S,I,S)",
    "S=SLANTP(H,H,H,I,S,S)",
    "S=SLANTR(H,H,H,I,I,S,I,S)",
    "SLANV2(S,S,S,S,S,S,S,S,S,S)",
    "SLAPLL(I,S,I,S,I,S)",
    "SLAPMT(L,I,I,S,I,I)",
    "S=SLAPY2(S,S)",
    "S=SLAPY3(S,S,S)",
    "SLAQGB(I,I,I,I,S,I,S,S,S,S,S,H)",
    "SLAQGE(I,I,S,I,S,S,S,S,S,H)",
    "SLAQP2(I,I,I,S,I,I,S,S,S,S)",
    "SLAQPS(I,I,I,I,I,S,I,I,S,S,S,S,S,I)",
    "SLAQR0(L,L,I,I,I,S,I,S,S,I,I,S,I,S,I,I)",
    "SLAQR1(I,S,I,S,S,S,S,S)",
    "SLAQR2(L,L,I,I,I,I,S,I,I,I,S,I,I,I,S,S,S,I,I,S,I,I,S,I,S,I)",
    "SLAQR3(L,L,I,I,I,I,S,I,I,I,S,I,I,I,S,S,S,I,I,S,I,I,S,I,S,I)",
    "SLAQR4(L,L,I,I,I,S,I,S,S,I,I,S,I,S,I,I)",
    "SLAQR5(L,L,I,I,I,I,I,S,S,S,I,I,I,S,I,S,I,S,I,I,S,I,I,S,I)",
    "SLAQSB(H,I,I,S,I,S,S,S,H)",
    "SLAQSP(H,I,S,S,S,S,H)",
    "SLAQSY(H,I,S,I,S,S,S,H)",
    "SLAQTR(L,L,I,S,I,S,S,S,S,S,I)",
    "SLAR1V(I,I,I,S,S,S,S,S,S,S,S,L,I,S,S,I,I,S,S,S,S)",
    "SLAR2V(I,S,S,S,I,S,S,I)",
    "SLARF(H,I,I,S,I,S,S,I,S)",
    "SLARFB(H,H,H,H,I,I,I,S,I,S,I,S,I,S,I)",
    "SLARFG(I,S,S,I,S)",
    "SLARFT(H,H,I,I,S,I,S,S,I)",
    "SLARFX(H,I,I,S,S,S,I,S)",
    "SLARGV(I,S,I,S,I,S,I)",
    "SLARNV(I,I,I,S)",
    "SLARRA(I,S,S,S,S,S,I,I,I)",
    "SLARRB(I,S,S,I,I,S,S,I,S,S,S,S,I,S,S,I,I)",
    "SLARRC(H,I,S,S,S,S,S,I,I,I,I)",
    "SLARRD(H,H,I,S,S,I,I,S,S,S,S,S,S,I,I,I,S,S,S,S,I,I,S,I,I)",
    "SLARRE(H,I,S,S,I,I,S,S,S,S,S,S,I,I,I,S,S,S,I,I,S,S,S,I,I)",
    "SLARRF(I,S,S,S,I,I,S,S,S,S,S,S,S,S,S,S,S,I)",
    "SLARRJ(I,S,S,I,I,S,I,S,S,S,I,S,S,I)",
    "SLARRK(I,I,S,S,S,S,S,S,S,S,I)",
    "SLARRR(I,S,S,I)",
    "SLARRV(I,S,S,S,S,S,I,I,I,I,S,S,S,S,S,S,I,I,S,S,I,I,S,I,I)",
    "SLARTG(S,S,S,S,S)",
    "SLARTV(I,S,I,S,I,S,S,I)",
    "SLARUV(I,I,S)",
    "SLARZ(H,I,I,I,S,I,S,S,I,S)",
    "SLARZB(H,H,H,H,I,I,I,I,S,I,S,I,S,I,S,I)",
    "SLARZT(H,H,I,I,S,I,S,S,I)",
    "SLAS2(S,S,S,S,S)",
    "SLASCL(H,I,I,S,S,I,I,S,I,I)",
    "SLASD0(I,I,S,S,S,I,S,I,I,I,S,I)",
    "SLASD1(I,I,I,S,S,S,S,I,S,I,I,I,S,I)",
    "SLASD2(I,I,I,I,S,S,S,S,S,I,S,I,S,S,I,S,I,I,I,I,I,I,I)",
    "SLASD3(I,I,I,I,S,S,I,S,S,I,S,I,S,I,S,I,I,I,S,I)",
    "SLASD4(I,I,S,S,S,S,S,S,I)",
    "SLASD5(I,S,S,S,S,S,S)",
    "SLASD6(I,I,I,I,S,S,S,S,S,I,I,I,I,I,S,I,S,S,S,S,I,S,S,S,I,I)",
    "SLASD7(I,I,I,I,I,S,S,S,S,S,S,S,S,S,S,I,I,I,I,I,I,I,S,I,S,S,I)",
    "SLASD8(I,I,S,S,S,S,S,S,I,S,S,I)",
    "SLASDA(I,I,I,I,S,S,S,I,S,I,S,S,S,S,I,I,I,I,S,S,S,S,I,I)",
    "SLASDQ(H,I,I,I,I,I,S,S,S,I,S,I,S,I,S,I)",
    "SLASDT(I,I,I,I,I,I,I)",
    "SLASET(H,I,I,S,S,S,I)",
    "SLASQ1(I,S,S,S,I)",
    "SLASQ2(I,S,I)",
    "SLASQ3(I,I,S,I,S,S,S,S,I,I,I,L)",
    "SLASQ4(I,I,S,I,I,S,S,S,S,S,S,S,I)",
    "SLASQ5(I,I,S,I,S,S,S,S,S,S,S,L)",
    "SLASQ6(I,I,S,I,S,S,S,S,S,S)",
    "SLASR(H,H,H,I,I,S,S,S,I)",
    "SLASRT(H,I,S,I)",
    "SLASSQ(I,S,I,S,S)",
    "SLASV2(S,S,S,S,S,S,S,S,S)",
    "SLASWP(I,S,I,I,I,I,I)",
    "SLASY2(L,L,I,I,I,S,I,S,I,S,I,S,S,I,S,I)",
    "SLASYF(H,I,I,I,S,I,I,S,I,I)",
    "SLATBS(H,H,H,H,I,I,S,I,S,S,S,I)",
    "SLATDF(I,I,S,I,S,S,S,I,I)",
    "SLATPS(H,H,H,H,I,S,S,S,S,I)",
    "SLATRD(H,I,I,S,I,S,S,S,I)",
    "SLATRS(H,H,H,H,I,S,I,S,S,S,I)",
    "SLATRZ(I,I,I,S,I,S,S)",
    "SLATZM(H,I,I,S,I,S,S,S,I,S)",
    "SLAUU2(H,I,S,I,I)",
    "SLAUUM(H,I,S,I,I)",
    "SLAZQ3(I,I,S,I,S,S,S,S,I,I,I,L,I,S,S,S,S,S,S)",
    "SLAZQ4(I,I,S,I,I,S,S,S,S,S,S,S,I,S)",
    "S=SNRM2(I,S,I)",
    "SOPGTR(H,I,S,S,S,I,S,I)",
    "SOPMTR(H,H,H,I,I,S,S,S,I,S,I)",
    "SORG2L(I,I,I,S,I,S,S,I)",
    "SORG2R(I,I,I,S,I,S,S,I)",
    "SORGBR(H,I,I,I,S,I,S,S,I,I)",
    "SORGHR(I,I,I,S,I,S,S,I,I)",
    "SORGL2(I,I,I,S,I,S,S,I)",
    "SORGLQ(I,I,I,S,I,S,S,I,I)",
    "SORGQL(I,I,I,S,I,S,S,I,I)",
    "SORGQR(I,I,I,S,I,S,S,I,I)",
    "SORGR2(I,I,I,S,I,S,S,I)",
    "SORGRQ(I,I,I,S,I,S,S,I,I)",
    "SORGTR(H,I,S,I,S,S,I,I)",
    "SORM2L(H,H,I,I,I,S,I,S,S,I,S,I)",
    "SORM2R(H,H,I,I,I,S,I,S,S,I,S,I)",
    "SORMBR(H,H,H,I,I,I,S,I,S,S,I,S,I,I)",
    "SORMHR(H,H,I,I,I,I,S,I,S,S,I,S,I,I)",
    "SORML2(H,H,I,I,I,S,I,S,S,I,S,I)",
    "SORMLQ(H,H,I,I,I,S,I,S,S,I,S,I,I)",
    "SORMQL(H,H,I,I,I,S,I,S,S,I,S,I,I)",
    "SORMQR(H,H,I,I,I,S,I,S,S,I,S,I,I)",
    "SORMR2(H,H,I,I,I,S,I,S,S,I,S,I)",
    "SORMR3(H,H,I,I,I,I,S,I,S,S,I,S,I)",
    "SORMRQ(H,H,I,I,I,S,I,S,S,I,S,I,I)",
    "SORMRZ(H,H,I,I,I,I,S,I,S,S,I,S,I,I)",
    "SORMTR(H,H,H,I,I,S,I,S,S,I,S,I,I)",
    "SPBCON(H,I,I,S,I,S,S,S,I,I)",
    "SPBEQU(H,I,I,S,I,S,S,S,I)",
    "SPBRFS(H,I,I,I,S,I,S,I,S,I,S,I,S,S,S,I,I)",
    "SPBSTF(H,I,I,S,I,I)",
    "SPBSV(H,I,I,I,S,I,S,I,I)",
    "SPBSVX(H,H,I,I,I,S,I,S,I,H,S,S,I,S,I,S,S,S,S,I,I)",
    "SPBTF2(H,I,I,S,I,I)",
    "SPBTRF(H,I,I,S,I,I)",
    "SPBTRS(H,I,I,I,S,I,S,I,I)",
    "SPOCON(H,I,S,I,S,S,S,I,I)",
    "SPOEQU(I,S,I,S,S,S,I)",
    "SPORFS(H,I,I,S,I,S,I,S,I,S,I,S,S,S,I,I)",
    "SPOSV(H,I,I,S,I,S,I,I)",
    "SPOSVX(H,H,I,I,S,I,S,I,H,S,S,I,S,I,S,S,S,S,I,I)",
    "SPOTF2(H,I,S,I,I)",
    "SPOTRF(H,I,S,I,I)",
    "SPOTRI(H,I,S,I,I)",
    "SPOTRS(H,I,I,S,I,S,I,I)",
    "SPPCON(H,I,S,S,S,S,I,I)",
    "SPPEQU(H,I,S,S,S,S,I)",
    "SPPRFS(H,I,I,S,S,S,I,S,I,S,S,S,I,I)",
    "SPPSV(H,I,I,S,S,I,I)",
    "SPPSVX(H,H,I,I,S,S,H,S,S,I,S,I,S,S,S,S,I,I)",
    "SPPTRF(H,I,S,I)",
    "SPPTRI(H,I,S,I)",
    "SPPTRS(H,I,I,S,S,I,I)",
    "SPTCON(I,S,S,S,S,S,I)",
    "SPTEQR(H,I,S,S,S,I,S,I)",
    "SPTRFS(I,I,S,S,S,S,S,I,S,I,S,S,S,I)",
    "SPTSV(I,I,S,S,S,I,I)",
    "SPTSVX(H,I,I,S,S,S,S,S,I,S,I,S,S,S,S,I)",
    "SPTTRF(I,S,S,I)",
    "SPTTRS(I,I,S,S,S,I,I)",
    "SPTTS2(I,I,S,S,S,I)",
    "SROT(I,S,I,S,I,S,S)",
    "SROTG(S,S,S,S)",
    "SROTM(I,S,I,S,I,S)",
    "SROTMG(S,S,S,S,S)",
    "SRSCL(I,S,S,I)",
    "SSBEV(H,H,I,I,S,I,S,S,I,S,I)",
    "SSBEVD(H,H,I,I,S,I,S,S,I,S,I,I,I,I)",
    "SSBEVX(H,H,H,I,I,S,I,S,I,S,S,I,I,S,I,S,S,I,S,I,I,I)",
    "SSBGST(H,H,I,I,I,S,I,S,I,S,I,S,I)",
    "SSBGV(H,H,I,I,I,S,I,S,I,S,S,I,S,I)",
    "SSBGVD(H,H,I,I,I,S,I,S,I,S,S,I,S,I,I,I,I)",
    "SSBGVX(H,H,H,I,I,I,S,I,S,I,S,I,S,S,I,I,S,I,S,S,I,S,I,I,I)",
    "SSBMV(H,I,I,S,S,I,S,I,S,S,I)",
    "SSBTRD(H,H,I,I,S,I,S,S,S,I,S,I)",
    "SSCAL(I,S,S,I)",
    "SSPCON(H,I,S,I,S,S,S,I,I)",
    "SSPEV(H,H,I,S,S,S,I,S,I)",
    "SSPEVD(H,H,I,S,S,S,I,S,I,I,I,I)",
    "SSPEVX(H,H,H,I,S,S,S,I,I,S,I,S,S,I,S,I,I,I)",
    "SSPGST(I,H,I,S,S,I)",
    "SSPGV(I,H,H,I,S,S,S,S,I,S,I)",
    "SSPGVD(I,H,H,I,S,S,S,S,I,S,I,I,I,I)",
    "SSPGVX(I,H,H,H,I,S,S,S,S,I,I,S,I,S,S,I,S,I,I,I)",
    "SSPMV(H,I,S,S,S,I,S,S,I)",
    "SSPR(H,I,S,S,I,S)",
    "SSPR2(H,I,S,S,I,S,I,S)",
    "SSPRFS(H,I,I,S,S,I,S,I,S,I,S,S,S,I,I)",
    "SSPSV(H,I,I,S,I,S,I,I)",
    "SSPSVX(H,H,I,I,S,S,I,S,I,S,I,S,S,S,S,I,I)",
    "SSPTRD(H,I,S,S,S,S,I)",
    "SSPTRF(H,I,S,I,I)",
    "SSPTRI(H,I,S,I,S,I)",
    "SSPTRS(H,I,I,S,I,S,I,I)",
    "SSTEBZ(H,H,I,S,S,I,I,S,S,S,I,I,S,I,I,S,I,I)",
    "SSTEDC(H,I,S,S,S,I,S,I,I,I,I)",
    "SSTEGR(H,H,I,S,S,S,S,I,I,S,I,S,S,I,I,S,I,I,I,I)",
    "SSTEIN(I,S,S,I,S,I,I,S,I,S,I,I,I)",
    "SSTEMR(H,H,I,S,S,S,S,I,I,I,S,S,I,I,I,L,S,I,I,I,I)",
    "SSTEQR(H,I,S,S,S,I,S,I)",
    "SSTERF(I,S,S,I)",
    "SSTEV(H,I,S,S,S,I,S,I)",
    "SSTEVD(H,I,S,S,S,I,S,I,I,I,I)",
    "SSTEVR(H,H,I,S,S,S,S,I,I,S,I,S,S,I,I,S,I,I,I,I)",
    "SSTEVX(H,H,I,S,S,S,S,I,I,S,I,S,S,I,S,I,I,I)",
    "SSWAP(I,S,I,S,I)",
    "SSYCON(H,I,S,I,I,S,S,S,I,I)",
    "SSYEV(H,H,I,S,I,S,S,I,I)",
    "SSYEVD(H,H,I,S,I,S,S,I,I,I,I)",
    "SSYEVR(H,H,H,I,S,I,S,S,I,I,S,I,S,S,I,I,S,I,I,I,I)",
    "SSYEVX(H,H,H,I,S,I,S,S,I,I,S,I,S,S,I,S,I,I,I,I)",
    "SSYGS2(I,H,I,S,I,S,I,I)",
    "SSYGST(I,H,I,S,I,S,I,I)",
    "SSYGV(I,H,H,I,S,I,S,I,S,S,I,I)",
    "SSYGVD(I,H,H,I,S,I,S,I,S,S,I,I,I,I)",
    "SSYGVX(I,H,H,H,I,S,I,S,I,S,S,I,I,S,I,S,S,I,S,I,I,I,I)",
    "SSYMM(H,H,I,I,S,S,I,S,I,S,S,I)",
    "SSYMV(H,I,S,S,I,S,I,S,S,I)",
    "SSYR(H,I,S,S,I,S,I)",
    "SSYR2(H,I,S,S,I,S,I,S,I)",
    "SSYR2K(H,H,I,I,S,S,I,S,I,S,S,I)",
    "SSYRFS(H,I,I,S,I,S,I,I,S,I,S,I,S,S,S,I,I)",
    "SSYRK(H,H,I,I,S,S,I,S,S,I)",
    "SSYSV(H,I,I,S,I,I,S,I,S,I,I)",
    "SSYSVX(H,H,I,I,S,I,S,I,I,S,I,S,I,S,S,S,S,I,I,I)",
    "SSYTD2(H,I,S,I,S,S,S,I)",
    "SSYTF2(H,I,S,I,I,I)",
    "SSYTRD(H,I,S,I,S,S,S,S,I,I)",
    "SSYTRF(H,I,S,I,I,S,I,I)",
    "SSYTRI(H,I,S,I,I,S,I)",
    "SSYTRS(H,I,I,S,I,I,S,I,I)",
    "STBCON(H,H,H,I,I,S,I,S,S,I,I)",
    "STBMV(H,H,H,I,I,S,I,S,I)",
    "STBRFS(H,H,H,I,I,I,S,I,S,I,S,I,S,S,S,I,I)",
    "STBSV(H,H,H,I,I,S,I,S,I)",
    "STBTRS(H,H,H,I,I,I,S,I,S,I,I)",
    "STGEVC(H,H,L,I,S,I,S,I,S,I,S,I,I,I,S,I)",
    "STGEX2(L,L,I,S,I,S,I,S,I,S,I,I,I,I,S,I,I)",
    "STGEXC(L,L,I,S,I,S,I,S,I,S,I,I,I,S,I,I)",
    "STGSEN(I,L,L,L,I,S,I,S,I,S,S,S,S,I,S,I,I,S,S,S,S,I,I,I,I)",
    "STGSJA(H,H,H,I,I,I,I,I,S,I,S,I,S,S,S,S,S,I,S,I,S,I,S,I,I)",
    "STGSNA(H,H,L,I,S,I,S,I,S,I,S,I,S,S,I,I,S,I,I,I)",
    "STGSY2(H,I,I,I,S,I,S,I,S,I,S,I,S,I,S,I,S,S,S,I,I,I)",
    "STGSYL(H,I,I,I,S,I,S,I,S,I,S,I,S,I,S,I,S,S,S,I,I,I)",
    "STPCON(H,H,H,I,S,S,S,I,I)",
    "STPMV(H,H,H,I,S,S,I)",
    "STPRFS(H,H,H,I,I,S,S,I,S,I,S,S,S,I,I)",
    "STPSV(H,H,H,I,S,S,I)",
    "STPTRI(H,H,I,S,I)",
    "STPTRS(H,H,H,I,I,S,S,I,I)",
    "STRCON(H,H,H,I,S,I,S,S,I,I)",
    "STREVC(H,H,L,I,S,I,S,I,S,I,I,I,S,I)",
    "STREXC(H,I,S,I,S,I,I,I,S,I)",
    "STRMM(H,H,H,H,I,I,S,S,I,S,I)",
    "STRMV(H,H,H,I,S,I,S,I)",
    "STRRFS(H,H,H,I,I,S,I,S,I,S,I,S,S,S,I,I)",
    "STRSEN(H,H,L,I,S,I,S,I,S,S,I,S,S,S,I,I,I,I)",
    "STRSM(H,H,H,H,I,I,S,S,I,S,I)",
    "STRSNA(H,H,L,I,S,I,S,I,S,I,S,S,I,I,S,I,I,I)",
    "STRSV(H,H,H,I,S,I,S,I)",
    "STRSYL(H,H,I,I,I,S,I,S,I,S,I,S,I)",
    "STRTI2(H,H,I,S,I,I)",
    "STRTRI(H,H,I,S,I,I)",
    "STRTRS(H,H,H,I,I,S,I,S,I,I)",
    "STZRQF(I,I,S,I,S,I)",
    "STZRZF(I,I,S,I,S,S,I,I)",
    "XERBLA(H,I)",
    "ZAXPY(I,Z,Z,I,Z,I)",
    "ZBDSQR(H,I,I,I,I,D,D,Z,I,Z,I,Z,I,D,I)",
    "ZCGESV(I,I,Z,I,I,Z,I,Z,I,Z,C,I,I)",
    "ZCOPY(I,Z,I,Z,I)",
    "Z=ZDOTC(I,Z,I,Z,I)",
    "Z=ZDOTU(I,Z,I,Z,I)",
    "ZDROT(I,Z,I,Z,I,D,D)",
    "ZDRSCL(I,D,Z,I)",
    "ZDSCAL(I,D,Z,I)",
    "ZGBBRD(H,I,I,I,I,I,Z,I,D,D,Z,I,Z,I,Z,I,Z,D,I)",
    "ZGBCON(H,I,I,I,Z,I,I,D,D,Z,D,I)",
    "ZGBEQU(I,I,I,I,Z,I,D,D,D,D,D,I)",
    "ZGBMV(H,I,I,I,I,Z,Z,I,Z,I,Z,Z,I)",
    "ZGBRFS(H,I,I,I,I,Z,I,Z,I,I,Z,I,Z,I,D,D,Z,D,I)",
    "ZGBSV(I,I,I,I,Z,I,I,Z,I,I)",
    "ZGBSVX(H,H,I,I,I,I,Z,I,Z,I,I,H,D,D,Z,I,Z,I,D,D,D,Z,D,I)",
    "ZGBTF2(I,I,I,I,Z,I,I,I)",
    "ZGBTRF(I,I,I,I,Z,I,I,I)",
    "ZGBTRS(H,I,I,I,I,Z,I,I,Z,I,I)",
    "ZGEBAK(H,H,I,I,I,D,I,Z,I,I)",
    "ZGEBAL(H,I,Z,I,I,I,D,I)",
    "ZGEBD2(I,I,Z,I,D,D,Z,Z,Z,I)",
    "ZGEBRD(I,I,Z,I,D,D,Z,Z,Z,I,I)",
    "ZGECON(H,I,Z,I,D,D,Z,D,I)",
    "ZGEEQU(I,I,Z,I,D,D,D,D,D,I)",
    "ZGEES(H,H,X,I,Z,I,I,Z,Z,I,Z,I,D,L,I)",
    "ZGEESX(H,H,X,H,I,Z,I,I,Z,Z,I,D,D,Z,I,D,L,I)",
    "ZGEEV(H,H,I,Z,I,Z,Z,I,Z,I,Z,I,D,I)",
    "ZGEEVX(H,H,H,H,I,Z,I,Z,Z,I,Z,I,I,I,D,D,D,D,Z,I,D,I)",
    "ZGEGS(H,H,I,Z,I,Z,I,Z,Z,Z,I,Z,I,Z,I,D,I)",
    "ZGEGV(H,H,I,Z,I,Z,I,Z,Z,Z,I,Z,I,Z,I,D,I)",
    "ZGEHD2(I,I,I,Z,I,Z,Z,I)",
    "ZGEHRD(I,I,I,Z,I,Z,Z,I,I)",
    "ZGELQ2(I,I,Z,I,Z,Z,I)",
    "ZGELQF(I,I,Z,I,Z,Z,I,I)",
    "ZGELS(H,I,I,I,Z,I,Z,I,Z,I,I)",
    "ZGELSD(I,I,I,Z,I,Z,I,D,D,I,Z,I,D,I,I)",
    "ZGELSS(I,I,I,Z,I,Z,I,D,D,I,Z,I,D,I)",
    "ZGELSX(I,I,I,Z,I,Z,I,I,D,I,Z,D,I)",
    "ZGELSY(I,I,I,Z,I,Z,I,I,D,I,Z,I,D,I)",
    "ZGEMM(H,H,I,I,I,Z,Z,I,Z,I,Z,Z,I)",
    "ZGEMV(H,I,I,Z,Z,I,Z,I,Z,Z,I)",
    "ZGEQL2(I,I,Z,I,Z,Z,I)",
    "ZGEQLF(I,I,Z,I,Z,Z,I,I)",
    "ZGEQP3(I,I,Z,I,I,Z,Z,I,D,I)",
    "ZGEQPF(I,I,Z,I,I,Z,Z,D,I)",
    "ZGEQR2(I,I,Z,I,Z,Z,I)",
    "ZGEQRF(I,I,Z,I,Z,Z,I,I)",
    "ZGERC(I,I,Z,Z,I,Z,I,Z,I)",
    "ZGERFS(H,I,I,Z,I,Z,I,I,Z,I,Z,I,D,D,Z,D,I)",
    "ZGERQ2(I,I,Z,I,Z,Z,I)",
    "ZGERQF(I,I,Z,I,Z,Z,I,I)",
    "ZGERU(I,I,Z,Z,I,Z,I,Z,I)",
    "ZGESC2(I,Z,I,Z,I,I,D)",
    "ZGESDD(H,I,I,Z,I,D,Z,I,Z,I,Z,I,D,I,I)",
    "ZGESV(I,I,Z,I,I,Z,I,I)",
    "ZGESVD(H,H,I,I,Z,I,D,Z,I,Z,I,Z,I,D,I)",
    "ZGESVX(H,H,I,I,Z,I,Z,I,I,H,D,D,Z,I,Z,I,D,D,D,Z,D,I)",
    "ZGETC2(I,Z,I,I,I,I)",
    "ZGETF2(I,I,Z,I,I,I)",
    "ZGETRF(I,I,Z,I,I,I)",
    "ZGETRI(I,Z,I,I,Z,I,I)",
    "ZGETRS(H,I,I,Z,I,I,Z,I,I)",
    "ZGGBAK(H,H,I,I,I,D,D,I,Z,I,I)",
    "ZGGBAL(H,I,Z,I,Z,I,I,I,D,D,D,I)",
    "ZGGES(H,H,H,X,I,Z,I,Z,I,I,Z,Z,Z,I,Z,I,Z,I,D,L,I)",
    "ZGGESX(H,H,H,X,H,I,Z,I,Z,I,I,Z,Z,Z,I,Z,I,D,D,Z,I,D,I,I,L,I)",
    "ZGGEV(H,H,I,Z,I,Z,I,Z,Z,Z,I,Z,I,Z,I,D,I)",
    "ZGGEVX(H,H,H,H,I,Z,I,Z,I,Z,Z,Z,I,Z,I,I,I,D,D,D,D,D,D,Z,I,D,I,L,I)",
    "ZGGGLM(I,I,I,Z,I,Z,I,Z,Z,Z,Z,I,I)",
    "ZGGHRD(H,H,I,I,I,Z,I,Z,I,Z,I,Z,I,I)",
    "ZGGLSE(I,I,I,Z,I,Z,I,Z,Z,Z,Z,I,I)",
    "ZGGQRF(I,I,I,Z,I,Z,Z,I,Z,Z,I,I)",
    "ZGGRQF(I,I,I,Z,I,Z,Z,I,Z,Z,I,I)",
    "ZGGSVD(H,H,H,I,I,I,I,I,Z,I,Z,I,D,D,Z,I,Z,I,Z,I,Z,D,I,I)",
    "ZGGSVP(H,H,H,I,I,I,Z,I,Z,I,D,D,I,I,Z,I,Z,I,Z,I,I,D,Z,Z,I)",
    "ZGTCON(H,I,Z,Z,Z,Z,I,D,D,Z,I)",
    "ZGTRFS(H,I,I,Z,Z,Z,Z,Z,Z,Z,I,Z,I,Z,I,D,D,Z,D,I)",
    "ZGTSV(I,I,Z,Z,Z,Z,I,I)",
    "ZGTSVX(H,H,I,I,Z,Z,Z,Z,Z,Z,Z,I,Z,I,Z,I,D,D,D,Z,D,I)",
    "ZGTTRF(I,Z,Z,Z,Z,I,I)",
    "ZGTTRS(H,I,I,Z,Z,Z,Z,I,Z,I,I)",
    "ZGTTS2(I,I,I,Z,Z,Z,Z,I,Z,I)",
    "ZHBEV(H,H,I,I,Z,I,D,Z,I,Z,D,I)",
    "ZHBEVD(H,H,I,I,Z,I,D,Z,I,Z,I,D,I,I,I,I)",
    "ZHBEVX(H,H,H,I,I,Z,I,Z,I,D,D,I,I,D,I,D,Z,I,Z,D,I,I,I)",
    "ZHBGST(H,H,I,I,I,Z,I,Z,I,Z,I,Z,D,I)",
    "ZHBGV(H,H,I,I,I,Z,I,Z,I,D,Z,I,Z,D,I)",
    "ZHBGVD(H,H,I,I,I,Z,I,Z,I,D,Z,I,Z,I,D,I,I,I,I)",
    "ZHBGVX(H,H,H,I,I,I,Z,I,Z,I,Z,I,D,D,I,I,D,I,D,Z,I,Z,D,I,I,I)",
    "ZHBMV(H,I,I,Z,Z,I,Z,I,Z,Z,I)",
    "ZHBTRD(H,H,I,I,Z,I,D,D,Z,I,Z,I)",
    "ZHECON(H,I,Z,I,I,D,D,Z,I)",
    "ZHEEV(H,H,I,Z,I,D,Z,I,D,I)",
    "ZHEEVD(H,H,I,Z,I,D,Z,I,D,I,I,I,I)",
    "ZHEEVR(H,H,H,I,Z,I,D,D,I,I,D,I,D,Z,I,I,Z,I,D,I,I,I,I)",
    "ZHEEVX(H,H,H,I,Z,I,D,D,I,I,D,I,D,Z,I,Z,I,D,I,I,I)",
    "ZHEGS2(I,H,I,Z,I,Z,I,I)",
    "ZHEGST(I,H,I,Z,I,Z,I,I)",
    "ZHEGV(I,H,H,I,Z,I,Z,I,D,Z,I,D,I)",
    "ZHEGVD(I,H,H,I,Z,I,Z,I,D,Z,I,D,I,I,I,I)",
    "ZHEGVX(I,H,H,H,I,Z,I,Z,I,D,D,I,I,D,I,D,Z,I,Z,I,D,I,I,I)",
    "ZHEMM(H,H,I,I,Z,Z,I,Z,I,Z,Z,I)",
    "ZHEMV(H,I,Z,Z,I,Z,I,Z,Z,I)",
    "ZHER(H,I,D,Z,I,Z,I)",
    "ZHER2(H,I,Z,Z,I,Z,I,Z,I)",
    "ZHER2K(H,H,I,I,Z,Z,I,Z,I,D,Z,I)",
    "ZHERFS(H,I,I,Z,I,Z,I,I,Z,I,Z,I,D,D,Z,D,I)",
    "ZHERK(H,H,I,I,D,Z,I,D,Z,I)",
    "ZHESV(H,I,I,Z,I,I,Z,I,Z,I,I)",
    "ZHESVX(H,H,I,I,Z,I,Z,I,I,Z,I,Z,I,D,D,D,Z,I,D,I)",
    "ZHETD2(H,I,Z,I,D,D,Z,I)",
    "ZHETF2(H,I,Z,I,I,I)",
    "ZHETRD(H,I,Z,I,D,D,Z,Z,I,I)",
    "ZHETRF(H,I,Z,I,I,Z,I,I)",
    "ZHETRI(H,I,Z,I,I,Z,I)",
    "ZHETRS(H,I,I,Z,I,I,Z,I,I)",
    "ZHGEQZ(H,H,H,I,I,I,Z,I,Z,I,Z,Z,Z,I,Z,I,Z,I,D,I)",
    "ZHPCON(H,I,Z,I,D,D,Z,I)",
    "ZHPEV(H,H,I,Z,D,Z,I,Z,D,I)",
    "ZHPEVD(H,H,I,Z,D,Z,I,Z,I,D,I,I,I,I)",
    "ZHPEVX(H,H,H,I,Z,D,D,I,I,D,I,D,Z,I,Z,D,I,I,I)",
    "ZHPGST(I,H,I,Z,Z,I)",
    "ZHPGV(I,H,H,I,Z,Z,D,Z,I,Z,D,I)",
    "ZHPGVD(I,H,H,I,Z,Z,D,Z,I,Z,I,D,I,I,I,I)",
    "ZHPGVX(I,H,H,H,I,Z,Z,D,D,I,I,D,I,D,Z,I,Z,D,I,I,I)",
    "ZHPMV(H,I,Z,Z,Z,I,Z,Z,I)",
    "ZHPR(H,I,D,Z,I,Z)",
    "ZHPR2(H,I,Z,Z,I,Z,I,Z)",
    "ZHPRFS(H,I,I,Z,Z,I,Z,I,Z,I,D,D,Z,D,I)",
    "ZHPSV(H,I,I,Z,I,Z,I,I)",
    "ZHPSVX(H,H,I,I,Z,Z,I,Z,I,Z,I,D,D,D,Z,D,I)",
    "ZHPTRD(H,I,Z,D,D,Z,I)",
    "ZHPTRF(H,I,Z,I,I)",
    "ZHPTRI(H,I,Z,I,Z,I)",
    "ZHPTRS(H,I,I,Z,I,Z,I,I)",
    "ZHSEIN(H,H,H,L,I,Z,I,Z,Z,I,Z,I,I,I,Z,D,I,I,I)",
    "ZHSEQR(H,H,I,I,I,Z,I,Z,Z,I,Z,I,I)",
    "ZLABRD(I,I,I,Z,I,D,D,Z,Z,Z,I,Z,I)",
    "ZLACGV(I,Z,I)",
    "ZLACN2(I,Z,Z,D,I,I)",
    "ZLACON(I,Z,Z,D,I)",
    "ZLACP2(H,I,I,D,I,Z,I)",
    "ZLACPY(H,I,I,Z,I,Z,I)",
    "ZLACRM(I,I,Z,I,D,I,Z,I,D)",
    "ZLACRT(I,Z,I,Z,I,Z,Z)",
    "Z=ZLADIV(Z,Z)",
    "ZLAED0(I,I,D,D,Z,I,Z,I,D,I,I)",
    "ZLAED7(I,I,I,I,I,I,D,Z,I,D,I,D,I,I,I,I,I,D,Z,D,I,I)",
    "ZLAED8(I,I,I,Z,I,D,D,I,D,D,Z,I,D,I,I,I,I,I,I,D,I)",
    "ZLAEIN(L,L,I,Z,I,Z,Z,Z,I,D,D,D,I)",
    "ZLAESY(Z,Z,Z,Z,Z,Z,Z,Z)",
    "ZLAEV2(Z,Z,Z,D,D,D,Z)",
    "ZLAG2C(I,I,Z,I,C,I,I)",
    "ZLAGS2(L,D,Z,D,D,Z,D,D,Z,D,Z,D,Z)",
    "ZLAGTM(H,I,I,D,Z,Z,Z,Z,I,D,Z,I)",
    "ZLAHEF(H,I,I,I,Z,I,I,Z,I,I)",
    "ZLAHQR(L,L,I,I,I,Z,I,Z,I,I,Z,I,I)",
    "ZLAHR2(I,I,I,Z,I,Z,Z,I,Z,I)",
    "ZLAHRD(I,I,I,Z,I,Z,Z,I,Z,I)",
    "ZLAIC1(I,I,Z,D,Z,Z,D,Z,Z)",
    "ZLALS0(I,I,I,I,I,Z,I,Z,I,I,I,I,I,D,I,D,D,D,D,I,D,D,D,I)",
    "ZLALSA(I,I,I,I,Z,I,Z,I,D,I,D,I,D,D,D,D,I,I,I,I,D,D,D,D,I,I)",
    "ZLALSD(H,I,I,I,D,D,Z,I,D,I,Z,D,I,I)",
    "D=ZLANGB(H,I,I,I,Z,I,D)",
    "D=ZLANGE(H,I,I,Z,I,D)",
    "D=ZLANGT(H,I,Z,Z,Z)",
    "D=ZLANHB(H,H,I,I,Z,I,D)",
    "D=ZLANHE(H,H,I,Z,I,D)",
    "D=ZLANHP(H,H,I,Z,D)",
    "D=ZLANHS(H,I,Z,I,D)",
    "D=ZLANHT(H,I,D,Z)",
    "D=ZLANSB(H,H,I,I,Z,I,D)",
    "D=ZLANSP(H,H,I,Z,D)",
    "D=ZLANSY(H,H,I,Z,I,D)",
    "D=ZLANTB(H,H,H,I,I,Z,I,D)",
    "D=ZLANTP(H,H,H,I,Z,D)",
    "D=ZLANTR(H,H,H,I,I,Z,I,D)",
    "ZLAPLL(I,Z,I,Z,I,D)",
    "ZLAPMT(L,I,I,Z,I,I)",
    "ZLAQGB(I,I,I,I,Z,I,D,D,D,D,D,H)",
    "ZLAQGE(I,I,Z,I,D,D,D,D,D,H)",
    "ZLAQHB(H,I,I,Z,I,D,D,D,H)",
    "ZLAQHE(H,I,Z,I,D,D,D,H)",
    "ZLAQHP(H,I,Z,D,D,D,H)",
    "ZLAQP2(I,I,I,Z,I,I,Z,D,D,Z)",
    "ZLAQPS(I,I,I,I,I,Z,I,I,Z,D,D,Z,Z,I)",
    "ZLAQR0(L,L,I,I,I,Z,I,Z,I,I,Z,I,Z,I,I)",
    "ZLAQR1(I,Z,I,Z,Z,Z)",
    "ZLAQR2(L,L,I,I,I,I,Z,I,I,I,Z,I,I,I,Z,Z,I,I,Z,I,I,Z,I,Z,I)",
    "ZLAQR3(L,L,I,I,I,I,Z,I,I,I,Z,I,I,I,Z,Z,I,I,Z,I,I,Z,I,Z,I)",
    "ZLAQR4(L,L,I,I,I,Z,I,Z,I,I,Z,I,Z,I,I)",
    "ZLAQR5(L,L,I,I,I,I,I,Z,Z,I,I,I,Z,I,Z,I,Z,I,I,Z,I,I,Z,I)",
    "ZLAQSB(H,I,I,Z,I,D,D,D,H)",
    "ZLAQSP(H,I,Z,D,D,D,H)",
    "ZLAQSY(H,I,Z,I,D,D,D,H)",
    "ZLAR1V(I,I,I,D,D,D,D,D,D,D,Z,L,I,D,D,I,I,D,D,D,D)",
    "ZLAR2V(I,Z,Z,Z,I,D,Z,I)",
    "ZLARCM(I,I,D,I,Z,I,Z,I,D)",
    "ZLARF(H,I,I,Z,I,Z,Z,I,Z)",
    "ZLARFB(H,H,H,H,I,I,I,Z,I,Z,I,Z,I,Z,I)",
    "ZLARFG(I,Z,Z,I,Z)",
    "ZLARFT(H,H,I,I,Z,I,Z,Z,I)",
    "ZLARFX(H,I,I,Z,Z,Z,I,Z)",
    "ZLARGV(I,Z,I,Z,I,D,I)",
    "ZLARNV(I,I,I,Z)",
    "ZLARRV(I,D,D,D,D,D,I,I,I,I,D,D,D,D,D,D,I,I,D,Z,I,I,D,I,I)",
    "ZLARTG(Z,Z,D,Z,Z)",
    "ZLARTV(I,Z,I,Z,I,D,Z,I)",
    "ZLARZ(H,I,I,I,Z,I,Z,Z,I,Z)",
    "ZLARZB(H,H,H,H,I,I,I,I,Z,I,Z,I,Z,I,Z,I)",
    "ZLARZT(H,H,I,I,Z,I,Z,Z,I)",
    "ZLASCL(H,I,I,D,D,I,I,Z,I,I)",
    "ZLASET(H,I,I,Z,Z,Z,I)",
    "ZLASR(H,H,H,I,I,D,D,Z,I)",
    "ZLASSQ(I,Z,I,D,D)",
    "ZLASWP(I,Z,I,I,I,I,I)",
    "ZLASYF(H,I,I,I,Z,I,I,Z,I,I)",
    "ZLATBS(H,H,H,H,I,I,Z,I,Z,D,D,I)",
    "ZLATDF(I,I,Z,I,Z,D,D,I,I)",
    "ZLATPS(H,H,H,H,I,Z,Z,D,D,I)",
    "ZLATRD(H,I,I,Z,I,D,Z,Z,I)",
    "ZLATRS(H,H,H,H,I,Z,I,Z,D,D,I)",
    "ZLATRZ(I,I,I,Z,I,Z,Z)",
    "ZLATZM(H,I,I,Z,I,Z,Z,Z,I,Z)",
    "ZLAUU2(H,I,Z,I,I)",
    "ZLAUUM(H,I,Z,I,I)",
    "ZPBCON(H,I,I,Z,I,D,D,Z,D,I)",
    "ZPBEQU(H,I,I,Z,I,D,D,D,I)",
    "ZPBRFS(H,I,I,I,Z,I,Z,I,Z,I,Z,I,D,D,Z,D,I)",
    "ZPBSTF(H,I,I,Z,I,I)",
    "ZPBSV(H,I,I,I,Z,I,Z,I,I)",
    "ZPBSVX(H,H,I,I,I,Z,I,Z,I,H,D,Z,I,Z,I,D,D,D,Z,D,I)",
    "ZPBTF2(H,I,I,Z,I,I)",
    "ZPBTRF(H,I,I,Z,I,I)",
    "ZPBTRS(H,I,I,I,Z,I,Z,I,I)",
    "ZPOCON(H,I,Z,I,D,D,Z,D,I)",
    "ZPOEQU(I,Z,I,D,D,D,I)",
    "ZPORFS(H,I,I,Z,I,Z,I,Z,I,Z,I,D,D,Z,D,I)",
    "ZPOSV(H,I,I,Z,I,Z,I,I)",
    "ZPOSVX(H,H,I,I,Z,I,Z,I,H,D,Z,I,Z,I,D,D,D,Z,D,I)",
    "ZPOTF2(H,I,Z,I,I)",
    "ZPOTRF(H,I,Z,I,I)",
    "ZPOTRI(H,I,Z,I,I)",
    "ZPOTRS(H,I,I,Z,I,Z,I,I)",
    "ZPPCON(H,I,Z,D,D,Z,D,I)",
    "ZPPEQU(H,I,Z,D,D,D,I)",
    "ZPPRFS(H,I,I,Z,Z,Z,I,Z,I,D,D,Z,D,I)",
    "ZPPSV(H,I,I,Z,Z,I,I)",
    "ZPPSVX(H,H,I,I,Z,Z,H,D,Z,I,Z,I,D,D,D,Z,D,I)",
    "ZPPTRF(H,I,Z,I)",
    "ZPPTRI(H,I,Z,I)",
    "ZPPTRS(H,I,I,Z,Z,I,I)",
    "ZPTCON(I,D,Z,D,D,D,I)",
    "ZPTEQR(H,I,D,D,Z,I,D,I)",
    "ZPTRFS(H,I,I,D,Z,D,Z,Z,I,Z,I,D,D,Z,D,I)",
    "ZPTSV(I,I,D,Z,Z,I,I)",
    "ZPTSVX(H,I,I,D,Z,D,Z,Z,I,Z,I,D,D,D,Z,D,I)",
    "ZPTTRF(I,D,Z,I)",
    "ZPTTRS(H,I,I,D,Z,Z,I,I)",
    "ZPTTS2(I,I,I,D,Z,Z,I)",
    "ZROT(I,Z,I,Z,I,D,Z)",
    "ZROTG(Z,Z,D,Z)",
    "ZSCAL(I,Z,Z,I)",
    "ZSPCON(H,I,Z,I,D,D,Z,I)",
    "ZSPMV(H,I,Z,Z,Z,I,Z,Z,I)",
    "ZSPR(H,I,Z,Z,I,Z)",
    "ZSPRFS(H,I,I,Z,Z,I,Z,I,Z,I,D,D,Z,D,I)",
    "ZSPSV(H,I,I,Z,I,Z,I,I)",
    "ZSPSVX(H,H,I,I,Z,Z,I,Z,I,Z,I,D,D,D,Z,D,I)",
    "ZSPTRF(H,I,Z,I,I)",
    "ZSPTRI(H,I,Z,I,Z,I)",
    "ZSPTRS(H,I,I,Z,I,Z,I,I)",
    "ZSTEDC(H,I,D,D,Z,I,Z,I,D,I,I,I,I)",
    "ZSTEGR(H,H,I,D,D,D,D,I,I,D,I,D,Z,I,I,D,I,I,I,I)",
    "ZSTEIN(I,D,D,I,D,I,I,Z,I,D,I,I,I)",
    "ZSTEMR(H,H,I,D,D,D,D,I,I,I,D,Z,I,I,I,L,D,I,I,I,I)",
    "ZSTEQR(H,I,D,D,Z,I,D,I)",
    "ZSWAP(I,Z,I,Z,I)",
    "ZSYCON(H,I,Z,I,I,D,D,Z,I)",
    "ZSYMM(H,H,I,I,Z,Z,I,Z,I,Z,Z,I)",
    "ZSYMV(H,I,Z,Z,I,Z,I,Z,Z,I)",
    "ZSYR(H,I,Z,Z,I,Z,I)",
    "ZSYR2K(H,H,I,I,Z,Z,I,Z,I,Z,Z,I)",
    "ZSYRFS(H,I,I,Z,I,Z,I,I,Z,I,Z,I,D,D,Z,D,I)",
    "ZSYRK(H,H,I,I,Z,Z,I,Z,Z,I)",
    "ZSYSV(H,I,I,Z,I,I,Z,I,Z,I,I)",
    "ZSYSVX(H,H,I,I,Z,I,Z,I,I,Z,I,Z,I,D,D,D,Z,I,D,I)",
    "ZSYTF2(H,I,Z,I,I,I)",
    "ZSYTRF(H,I,Z,I,I,Z,I,I)",
    "ZSYTRI(H,I,Z,I,I,Z,I)",
    "ZSYTRS(H,I,I,Z,I,I,Z,I,I)",
    "ZTBCON(H,H,H,I,I,Z,I,D,Z,D,I)",
    "ZTBMV(H,H,H,I,I,Z,I,Z,I)",
    "ZTBRFS(H,H,H,I,I,I,Z,I,Z,I,Z,I,D,D,Z,D,I)",
    "ZTBSV(H,H,H,I,I,Z,I,Z,I)",
    "ZTBTRS(H,H,H,I,I,I,Z,I,Z,I,I)",
    "ZTGEVC(H,H,L,I,Z,I,Z,I,Z,I,Z,I,I,I,Z,D,I)",
    "ZTGEX2(L,L,I,Z,I,Z,I,Z,I,Z,I,I,I)",
    "ZTGEXC(L,L,I,Z,I,Z,I,Z,I,Z,I,I,I,I)",
    "ZTGSEN(I,L,L,L,I,Z,I,Z,I,Z,Z,Z,I,Z,I,I,D,D,D,Z,I,I,I,I)",
    "ZTGSJA(H,H,H,I,I,I,I,I,Z,I,Z,I,D,D,D,D,Z,I,Z,I,Z,I,Z,I,I)",
    "ZTGSNA(H,H,L,I,Z,I,Z,I,Z,I,Z,I,D,D,I,I,Z,I,I,I)",
    "ZTGSY2(H,I,I,I,Z,I,Z,I,Z,I,Z,I,Z,I,Z,I,D,D,D,I)",
    "ZTGSYL(H,I,I,I,Z,I,Z,I,Z,I,Z,I,Z,I,Z,I,D,D,Z,I,I,I)",
    "ZTPCON(H,H,H,I,Z,D,Z,D,I)",
    "ZTPMV(H,H,H,I,Z,Z,I)",
    "ZTPRFS(H,H,H,I,I,Z,Z,I,Z,I,D,D,Z,D,I)",
    "ZTPSV(H,H,H,I,Z,Z,I)",
    "ZTPTRI(H,H,I,Z,I)",
    "ZTPTRS(H,H,H,I,I,Z,Z,I,I)",
    "ZTRCON(H,H,H,I,Z,I,D,Z,D,I)",
    "ZTREVC(H,H,L,I,Z,I,Z,I,Z,I,I,I,Z,D,I)",
    "ZTREXC(H,I,Z,I,Z,I,I,I,I)",
    "ZTRMM(H,H,H,H,I,I,Z,Z,I,Z,I)",
    "ZTRMV(H,H,H,I,Z,I,Z,I)",
    "ZTRRFS(H,H,H,I,I,Z,I,Z,I,Z,I,D,D,Z,D,I)",
    "ZTRSEN(H,H,L,I,Z,I,Z,I,Z,I,D,D,Z,I,I)",
    "ZTRSM(H,H,H,H,I,I,Z,Z,I,Z,I)",
    "ZTRSNA(H,H,L,I,Z,I,Z,I,Z,I,D,D,I,I,Z,I,D,I)",
    "ZTRSV(H,H,H,I,Z,I,Z,I)",
    "ZTRSYL(H,H,I,I,I,Z,I,Z,I,Z,I,D,I)",
    "ZTRTI2(H,H,I,Z,I,I)",
    "ZTRTRI(H,H,I,Z,I,I)",
    "ZTRTRS(H,H,H,I,I,Z,I,Z,I,I)",
    "ZTZRQF(I,I,Z,I,Z,I)",
    "ZTZRZF(I,I,Z,I,Z,Z,I,I)",
    "ZUNG2L(I,I,I,Z,I,Z,Z,I)",
    "ZUNG2R(I,I,I,Z,I,Z,Z,I)",
    "ZUNGBR(H,I,I,I,Z,I,Z,Z,I,I)",
    "ZUNGHR(I,I,I,Z,I,Z,Z,I,I)",
    "ZUNGL2(I,I,I,Z,I,Z,Z,I)",
    "ZUNGLQ(I,I,I,Z,I,Z,Z,I,I)",
    "ZUNGQL(I,I,I,Z,I,Z,Z,I,I)",
    "ZUNGQR(I,I,I,Z,I,Z,Z,I,I)",
    "ZUNGR2(I,I,I,Z,I,Z,Z,I)",
    "ZUNGRQ(I,I,I,Z,I,Z,Z,I,I)",
    "ZUNGTR(H,I,Z,I,Z,Z,I,I)",
    "ZUNM2L(H,H,I,I,I,Z,I,Z,Z,I,Z,I)",
    "ZUNM2R(H,H,I,I,I,Z,I,Z,Z,I,Z,I)",
    "ZUNMBR(H,H,H,I,I,I,Z,I,Z,Z,I,Z,I,I)",
    "ZUNMHR(H,H,I,I,I,I,Z,I,Z,Z,I,Z,I,I)",
    "ZUNML2(H,H,I,I,I,Z,I,Z,Z,I,Z,I)",
    "ZUNMLQ(H,H,I,I,I,Z,I,Z,Z,I,Z,I,I)",
    "ZUNMQL(H,H,I,I,I,Z,I,Z,Z,I,Z,I,I)",
    "ZUNMQR(H,H,I,I,I,Z,I,Z,Z,I,Z,I,I)",
    "ZUNMR2(H,H,I,I,I,Z,I,Z,Z,I,Z,I)",
    "ZUNMR3(H,H,I,I,I,I,Z,I,Z,Z,I,Z,I)",
    "ZUNMRQ(H,H,I,I,I,Z,I,Z,Z,I,Z,I,I)",
    "ZUNMRZ(H,H,I,I,I,I,Z,I,Z,Z,I,Z,I,I)",
    "ZUNMTR(H,H,H,I,I,Z,I,Z,Z,I,Z,I,I)",
    "ZUPGTR(H,I,Z,Z,Z,I,Z,I)",
    "ZUPMTR(H,H,H,I,I,Z,Z,Z,I,Z,I)",
  };
  
  char **pszProto = bsearch(&szFuncName, pszProtos,
                            sizeof(pszProtos)/sizeof(char*), 
                            sizeof(char*), Compare);

  return (pszProto) ? *pszProto : 0;
}
