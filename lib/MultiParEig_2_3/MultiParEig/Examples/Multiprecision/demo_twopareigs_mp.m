%DEMO_TWOPAREIGS_MP  demo nonsingular two-parameter eigenvalue problem 
%
% We solve a random complex two-parameter eigenvalue problem
%
% A1 x = lambda B1 x + mu C1 x 
% A2 y = lambda B2 y + mu C2 y,
%
% with random multiprecision matrices of size n x n (n = 100) using twopareigs
%
% First option, which is faster and should work for most problems, is to
% first compute double precision approximations and refine them using TRQI
% in multiprecision
%
% Second option is to apply Krylov-Schur in multiprecision, which is slower, 
% but might be better for certain very ill conditioned problems.

% See also: DEMO_TWOPAREIGS, TWOPAREIGS_KS

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision: 2.1.2017

is_numeric_type_supported('mp'); % check if MCT is installed

n = 100;
neig = 10;
rng('default')
rng(3)
A1 = rand(n,'mp') + 1i*rand(n,'mp'); B1 = rand(n,'mp') + 1i*rand(n,'mp'); C1 = rand(n,'mp') + 1i*rand(n,'mp');
A2 = rand(n,'mp') + 1i*rand(n,'mp'); B2 = rand(n,'mp') + 1i*rand(n,'mp'); C2 = rand(n,'mp') + 1i*rand(n,'mp');

% First option : use double precision and then refine in multiprecision 
% As matrices are already in multiprecision, no options are needed, as this
% is the default behaviour of twopareigs
tic
[lambda1,mu1] = twopareigs(A1,B1,C1,A2,B2,C2,neig);
time1 = toc

eigenvalues1 = [lambda1 mu1]

% Second option : use multiprecision from the beginning
tic
opts = [];
opts.init_double = 0;
[lambda2,mu2] = twopareigs(A1,B1,C1,A2,B2,C2,neig,opts);
time2 = toc

eigenvalues2 = [lambda2 mu2]


