function [varargout] = lapack(varargin);
%LAPACK Call any LAPACK or BLAS function or subroutine.  
%   C = LAPACK('FUNC', A1, A2, ...) 
%   [R1, R2, ...] = LAPACK('FUNC', A1, A2, ...)
%   calls the BLAS/LAPACK function named FUNC, passing it arguments
%   A1, A2, ..., then returns the results in the cell array C or in
%   the arguments R1, R2, ...
%
%   PROTO = LAPACK('FUNC') returns the prototype for the BLAS/LAPACK
%   function named FUNC as a string.
%
%   [R1, R2, ...] = LAPACK(PROTO, A1, A2, ...)
%   C = LAPACK(PROTO, A1, A2, ...)
%   where PROTO is a string similar to that returned by LAPACK('FUNC').
%   This form allows selection of a subset of the arguments to be
%   returned.  It also allows using functions from a newer version of
%   the BLAS/LAPACK libraries where the prototype is not known.
%
%   The PROTO argument is a string of the form 'T=FUNC(T,T,T,...)',
%   where T indicates the FORTRAN type of each argument to the routine
%   as well as the type of any return value. Lowercase values indicate
%   that this is an input value only, while uppercase values indicate
%   that it is an output (and possibly also an input) value.  Valid
%   values for T are
%
%      T    FORTRAN Type       Acceptable Matlab Type
%   --------------------------------------------------------------------------
%     s,S - REAL               single, double*
%     d,D - DOUBLE PRECISION   double*
%     c,C - COMPLEX            single, complex single*, double, complex double 
%     z,Z - COMPLEX*16         double, complex double*
%     i,I - INTEGER            int32*, int64*, double, logical
%     l,L - LOGICAL            int32*, int64, double, logical
%     h,H - CHARACTER          char*
%            * indicates the type that will be returned
%
%   Even if the BLAS/LAPACK routine specifies something as an output,
%   is is okay to specify it as input only (lowercase) if you are not
%   interested in that result.
%
%   The easiest place to find the descriptions of the arguments for
%   BLAS/LAPACK functions is from the FORTRAN source code, which
%   include comments that thoroughly describe each argument.  
%   They are available through the LAPACKHELP command.
%
% Examples:
%
%   Call the BLAS routine xDOT and xDOTC to compute inner product of X and Y
%     dot = lapack('D=DDOT(i,d,i,d,i)', length(X), X, 1, Y, 1)
%     dot = lapack('S=SDOT(i,s,i,s,i)', length(X), X, 1, Y, 1)
%     dot = lapack('Z=zdotc(i,z,i,z,i)', length(X), X, 1, Y, 1)
%   which are all equivalent to the Matlab expression
%     dot = X'*Y
%   the first lapack statement is also equivalent to
%     C = lapack('ddot', length(X), X, 1, Y, 1);
%     dot = C{1}
%
%   Call the LAPACK routine DGESVD to compute the SVD of the m by n matrix X
%     [m,n] = size(X);   % determine the dimensions of X
%     C = lapack('dgesvd', 'A', 'A', m, n, X, m, zeros(n,1), ...
%                zeros(m), m, zeros(n), n, zeros(5*m,1), 5*m, 0);
%     [s,U,VT] = C{[7,8,10]};
%     V = VT';
%   which is equivalent to the Matlab expressions:
%     [U,S,V] = svd(X);
%     s = diag(S);
%   equivalently we could specify a prototype for DGESVD:
%     [s,U,VT] = lapack('DGESVD(h,h,i,i,d,i,D,D,i,D,i,d,i,i)', ...
%                     'A', 'A', m, n, X, m, zeros(n,1), ...
%                      zeros(m), m, zeros(n), n, zeros(5*m,1), 5*m, 0);
%     V = VT';
%   to see a description of what the arguments for DGESVD should be type:
%     lapackhelp('dgesvd');
%
%   Call the LAPACK function LSAME to compare two characters
%     lapack('L=LSAME(h,h)', 'x', 'y')
%   which is equivalent to the Matlab expression
%     'x' == 'y'
%
%   Call the BLAS routine DGEMM
%     C = lapack('DGEMM',TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC);
%   the cell array C will have 13 elements corresponding to the updated
%   arguments passed to the DGEMM function.  To determine what each
%   argument stands for, type lapackhelp('dgemm') in Matlab.
%
% See also LAPACKHELP.

%   Copyright 2007-2009, Timothy Toolan, University of Rhode Island

% If we are in this function, then an appropriate compiled version
% does not exist, so ask the user if they want to build one.
reply = input(['The function ''lapack'' requires a compiled mex version, ' ...
               'create now? (Y/n): '], 's');

if (isempty(reply) || strcmpi(reply(1), 'y'))
  % To ensure the compiled version has a higher precedence than this
  % m-file, it should be placed in the directory where this m-file
  % resides.  We will also assume the C source file resides there.
  [dir, name, ext] = fileparts(which('lapack'));
  if (~length(dir) || ~strcmpi(name, 'lapack'))
    error('Could not determine directory containing this script (lapack.m).');
  end

  % On UNIX, we need to link to the library 'libdl' to allow runtime
  % linking to the internal Matlab LAPACK and BLAS functions.
  libs = '';
  if (isunix) 
    libs = '-ldl'; 
  end
  
  largedims = '';
  [arch, maxsize] = computer;
  if (maxsize > 2^31)
    largedims = '-largeArrayDims';
  end
  
  % Build the executable in the directory that contains this m-file.
  olddir = pwd;
  cd(dir);
  try
    eval(['mex ' largedims ' ' libs  ' lapack.c']); 
  catch exception
    cd(olddir);
    rethrow(exception);
  end
  cd(olddir);
  
  % Verify that we sucessfully created the mex file before we try to call it.
  [dir, name, ext] = fileparts(which('lapack'));
  if (~strcmpi(name, 'lapack') || ~strcmpi(ext, ['.' mexext]))
    error(['Could not create the compiled version of ''lapack''.' ...
           '  Try creating manually.']);
  end
  
  % Call the newly compiled function.
  [varargout{1:nargout}] = lapack(varargin{:});
else
  % If the user pressed 'n' to complation, so throw an error.
  error('Please create a compiled mex version of the function ''lapack''.');
end




