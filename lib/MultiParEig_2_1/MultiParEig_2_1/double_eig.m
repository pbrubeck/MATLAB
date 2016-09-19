function lambda = double_eig(A, B, opts)

%DOUBLE_EIG   lambda such that A + lambda*B has a multiple eigenvalue
%
% lambda = DOUBLE_EIG(A, B, opts) return values lambda where 
% A + lambda B has a multiple eigenvalue
% 
% The algorithm translates the problem to a singular quadratic two-parameter eigenvvalue problem
% (A + lambda*B + mu*I)*x = 0
% (A + lambda*B + mu*I)^2*y = 0
% and solves this using a linearization and the staircase method for singular MEP
%
% Matrices should not be larger then 10 x 10, as the staircase algorithm
% is very sensitive.
%
% If you do not get as many results as expected, try to use opts with
% larger rankeps, e.g., opts.rankeps = 1e-8 or smaller

% Reference: A. Muhic, B. Plestenjak:  A method for computing all values 
% lambda such that A + lambda*B has a multiple eigenvalue, Linear Algebra 
% Appl. 440 (2014) 345-359.

% MultiParEig toolbox
% B. Plestenjak and A. Muhic, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision 8.9.2015

n = size(A,1);
I = eye(n);
if strcmp(class(A), 'sym') I = vpa(I); end

opts.singular = 1;

% Linearization to a singular two-parameter eigenvalue problem
[P,Q,R] = linearize_quadtwopar(A^2, A*B+B*A, 2*A, B^2, 2*B, I);
lambda = twopareig(A,-B,-I,P,-Q,-R,opts); 


