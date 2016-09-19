%DEMO_BIVARIATE   Demo for zeros of a system of bivariate polynomials 
%
% To solve a system of two bivariate polynomials
%
% p(x,y) = 1 + 2x + 3y + 4x^2 + 5xy + 6y^2 + 7x^3 + 8x^2y + 9xy^2 + 10y^3 = 0
% q(x,y) = 10 + 9x + 8y + 7x^2 + 6xy + 5y^2 + 4x^3 + 3x^2y + 2xy^2 +  y^3 = 0
%
% we linearize it as a singular two-parameter eigenvalue problem with 
% matrices of size 5 x 5 using Linearization 1 from the reference. The finite 
% eigenvalues are the zeroes. See the reference for linearizations of the 
% generic bivariate polynomial.
%
% Reference: B. Plestenjak, M. E. Hochstenbach: Roots of bivariate polynomial 
% systems via determinantal representations, arXiv:1506.02291, 2015
%
% See also: TWOPAREIG, DEMO_TRIVARIATE

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision 8.9.2015

A1 = [1     0     0     0     6
      0    -1     0     0     0
      0     0    -1     0     0
      0     0     0    -1     0
      0     0     0     0    -1];

B1 = [2     4     7     0     9
      1     0     0     0     0
      0     1     0     0     0
      0     0     0     0     0
      0     0     0     0     0];

C1 = [3     5     8     0    10
      0     0     0     0     0
      0     0     0     0     0
      1     0     0     0     0
      0     0     0     1     0];

A2 = [10     0     0     0     5
      0    -1     0     0     0
      0     0    -1     0     0
      0     0     0    -1     0
      0     0     0     0    -1];

B2 = [9     7     4     0     2
      1     0     0     0     0
      0     1     0     0     0
      0     0     0     0     0
      0     0     0     0     0];

C2 = [8     6     3     0     1
      0     0     0     0     0
      0     0     0     0     0
      1     0     0     0     0
      0     0     0     1     0];

opts = [];
opts.singular = 1;
[x,y] = twopareig(A1,B1,C1,A2,B2,C2,opts)