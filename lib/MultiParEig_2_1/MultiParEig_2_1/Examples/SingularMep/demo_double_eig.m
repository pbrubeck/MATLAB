%DEMO_DOUBLE_EIG    Demo example for double_eig
%
% Example 1: For matrices
%
% A1 = [1 -2 3; -1 1 2;  1 1 -1]
% B1 = [1 -1 1;  1 1 3; -1 1  2]
%
% we find all values lambda such that A + lambda*B has a multiple eigenvalue
%
% Example 2: For matrices 
%
% A2 = [1 -2 3; -1 1 2;  1 1 -1]
% B2 = diag([2 2 3]) - A2
%
% we find all values lambda such that A + lambda*B has a multiple eigenvalue
%
% See also: DOUBLE_EIG

% Reference: A. Muhic, B. Plestenjak:  A method for computing all values lambda 
% such that A + lambda*B has a multiple eigenvalue, Linear Algebra Appl. 440 (2014) 345-359.

% MultiParEig toolbox
% B. Plestenjak and A. Muhic, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision 8.9.2015

A1 = [1 -2 3; -1 1 2; 1 1 -1]
B1 = [1 -1 1; 1 1 3; -1 1 2]

lambda1 = double_eig(A1,B1)

A2 = A1
B2 = diag([2 2 3]) - A2
 
lambda2 = double_eig(A2,B2)

% confirmation that A + lambda*B has multiple eigenvalues
eigs1 = eig(A1 + lambda1(1)*B1)
eigs2 = eig(A2 + lambda2(1)*B2)

% if you want to see the steps of the staircase algorithm, 
% use opts.showrank = 1

opts=[];
opts.showrank = 1;
lambda1 = double_eig(A1,B1,opts)
lambda2 = double_eig(A2,B2,opts)
