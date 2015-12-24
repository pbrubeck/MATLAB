function lapackhelp(func)
%LAPACKHELP   View online documentation for LAPACK and BLAS.
%   LAPACKHELP(ROUTINE) displays the documented source code for the
%   BLAS/LAPACK routine named ROUTINE.  This contains detailed
%   descriptions about what the routine does and what the arguments
%   should be.
%
%   LAPACKHELP()           display index of all BLAS/LAPACK routines.
%
%   LAPACKHELP('-lapack')  display the LAPACK Home Page.
%   LAPACKHELP('-blas')    display the BLAS Home Page.
%
%   LAPACKHELP('-lug')     display the LAPACK Users' Guide.
%   LAPACKHELP('-lugcomp') descriptions of all Driver and Computational 
%                              Routines from the LAPACK Users' Guide.
%   LAPACKHELP('-lugaux')  descriptions of all Auxilary Routines from 
%                              the LAPACK Users' Guide.
%   LAPACKHELP('-lugblas') Quick Reference Guide to the BLAS from 
%                              the LAPACK Users' Guide.
%
%   Example
%     lapackhelp();
%     lapackhelp('-lug');
%     lapackhelp('DGEMM');
%     lapackhelp('DGESVD');
%
%   See also LAPACK.

%   Copyright 2009, Timothy Toolan, University of Rhode Island

if (nargin < 1)
  webpage = 'http://www.netlib.org/lapack/explore-html/';
elseif (strcmpi(func, '-lapack'))
  webpage = 'http://www.netlib.org/lapack';
elseif (strcmpi(func, '-blas'))
  webpage = 'http://www.netlib.org/blas';
elseif (strcmpi(func, '-lug'))
  webpage = 'http://www.netlib.org/lapack/lug';
elseif (strcmpi(func, '-lugcomp'))
  webpage = 'http://www.netlib.org/lapack/lug/node142.html';
elseif (strcmpi(func, '-lugaux'))
  webpage = 'http://www.netlib.org/lapack/lug/node144.html';
elseif (strcmpi(func, '-lugblas'))
  webpage = 'http://www.netlib.org/lapack/lug/node145.html';
else    
  webpage = ['http://www.netlib.org/lapack/explore-html/' ...
             lower(func) '.f.html'];
end

web(webpage);
