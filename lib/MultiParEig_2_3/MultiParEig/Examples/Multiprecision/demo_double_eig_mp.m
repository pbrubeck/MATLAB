%DEMO_DOUBLE_EIG_MP    Demo example for double_eig in multiprecision
%
% For a pair of 10 x 10 matrices A, B we search for all values lambda such 
% that A + lambda*B has a multiple eigenvalue. We use double_eig that 
% solves the related singular two-parameter eigenvalue problem. Double 
% precision is not enough to solve the problem. The staircase algorithm for
% singular two-parameter eigenvalue problems fails as gaps between singular 
% values become too small.
%
% If we repeat the computation in quadruple precision using mp, we get the
% solution. As a faster alternative in quadruple precision we can also use 
% rank revealing qr instead of SVD.
%
% This example requires Multiprecision Computing Toolbox for MATLAB, see
% http://www.advanpix.com/
%
% See also: DOUBLE_EIG

% Reference: A. Muhic, B. Plestenjak:  A method for computing all values 
% lambda such that A + lambda*B has a multiple eigenvalue, Linear Algebra 
% Appl. 440 (2014) 345-359.

% MultiParEig toolbox
% B. Plestenjak and A. Muhic, University of Ljubljana
% P. Holoborodko, Advanpix LLC.
% FreeBSD License, see LICENSE.txt

% PH 26.11.2016: fixed bug in accuracy check.
% BP 27.11.2016: espace from legacy mode for rand, test if MCT is installed

% Last revision 27.11.2016

% this gives matrices such that the computation fails in double precision
rng('default')
rng(3) 
A1 = rand(10);
B1 = rand(10);

% 1) we try double_eig in double precision
disp('Fast computation in double precision fails, no solutions are found')
disp('------------------------------------------------------------------')
tic
opts = [];
opts.showrank = 1;
lambda1 = double_eig(A1,B1,opts);
ans1 = length(lambda1) % number of solutions found
t1 = toc

% 2) we use double_eig in quadruple precision and SVD
is_numeric_type_supported('mp'); % check if MCT is installed
disp('Computation in quadruple precision works, compare columns gap, d(r), d(r+1) to double precision values')
disp('------------------------------------------------------------------------------------------------------')
tic
A2 = mp(A1);
B2 = mp(B1);
lambda2 = double_eig(A2,B2,opts);
ans2 = length(lambda2)
t2 = toc

% 3) we use double_eig in quadruple precision and rank revealing QR
disp('Computation in quadruple precision with faster rank revealing qr works as well')
disp('------------------------------------------------------------------------------')
tic
opts.rrqr = 1;
lambda3 = double_eig(A2,B2,opts);
ans3 = length(lambda3)
t3 = toc

fprintf('One solution is lambda(1) = %e\n',lambda3(1))
disp('Eigenvalues of A+lambda(1)*B are (there should be a double eigenvalue):')
num2str(eig(A2+lambda3(1)*B2))
