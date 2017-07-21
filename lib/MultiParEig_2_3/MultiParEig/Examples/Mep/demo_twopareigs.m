%DEMO_TWOPAREIGS  demo nonsingular two-parameter eigenvalue problem with small matrices
%
% We solve a two-parameter eigenvalue problem
%
% A1 x = lambda B1 x + mu C1 x 
% A2 y = lambda B2 y + mu C2 y,
%
% with random matrices of size n x n (n = 40) using twopareig and twopareigs
%
% twopareig is a "direct method" that computes all n^2 eigenvalues
% twopareigs_ira and twopareigs_ks are iterative subspace methods that compute
% just neig eigenvalues and can tackle much larger problems than twopareig
%
% A sample output is
%
% time1 =
%     5.2731
% eigenvalues1 =
%   -0.4836 + 0.0000i  -0.0032 + 0.0000i
%    2.8407 + 0.0000i  -0.0184 + 0.0000i
%   -2.9606 + 0.5203i   0.0309 - 0.0253i
%   -2.9606 - 0.5203i   0.0309 + 0.0253i
%    0.6115 + 0.0704i  -0.0081 - 0.0394i
%    0.6115 - 0.0704i  -0.0081 + 0.0394i
%    0.1356 + 0.3016i   0.0234 - 0.0357i
%    0.1356 - 0.3016i   0.0234 + 0.0357i
%   -1.5848 + 0.0000i   0.0449 + 0.0000i
%   -0.1568 + 0.0487i   0.0550 + 0.0581i
% time2 =
%     0.1111
% eigenvalues2 =
%   -0.4836 + 0.0000i  -0.0032 + 0.0000i
%    2.8407 + 0.0000i  -0.0184 + 0.0000i
%    0.6115 + 0.0704i  -0.0081 - 0.0394i
%    0.6115 - 0.0704i  -0.0081 + 0.0394i
%   -2.9606 + 0.5203i   0.0309 - 0.0253i
%   -2.9606 - 0.5203i   0.0309 + 0.0253i
%    0.1356 + 0.3016i   0.0234 - 0.0357i
%    0.1356 - 0.3016i   0.0234 + 0.0357i
%   -1.5848 + 0.0000i   0.0449 + 0.0000i
%   -0.1568 - 0.0487i   0.0550 - 0.0581i
% time3 =
%     0.0788
% eigenvalues3 =
%   -0.4836 + 0.0000i  -0.0032 + 0.0000i
%    2.8407 + 0.0000i  -0.0184 + 0.0000i
%   -2.9606 + 0.5203i   0.0309 - 0.0253i
%   -2.9606 - 0.5203i   0.0309 + 0.0253i
%    0.6115 + 0.0704i  -0.0081 - 0.0394i
%    0.6115 - 0.0704i  -0.0081 + 0.0394i
%    0.1356 + 0.3016i   0.0234 - 0.0357i
%    0.1356 - 0.3016i   0.0234 + 0.0357i
%   -1.5848 + 0.0000i   0.0449 + 0.0000i
%   -0.1937 + 0.0000i   0.1126 + 0.0000i  
%
% See also: TWOPAREIG, TWOPAREIGS_IRA, TWOPAREIGS_KS

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision: 8.9.2015

rng('default')
rng(1);

n = 40;
neig = 10;

A1 = rand(n); B1 = rand(n); C1 = rand(n);
A2 = rand(n); B2 = rand(n); C2 = rand(n);

% we solve the two-parameter eigenvalue problem with "direct" method
tic
[lambda,mu] = twopareig(A1,B1,C1,A2,B2,C2);
time1 = toc

% find the neig eigenvalues with the smallest value |mu|
[tmp,ord] = sort(abs(mu));

eigenvalues1 = [lambda(ord(1:neig)) mu(ord(1:neig))]

% we use iterative method (implicitly restarted Arnoldi) to compute just few eigenvalues with the smallest |mu|
tic
[lambda2,mu2] = twopareigs_ira(A1,B1,C1,A2,B2,C2,neig);
time2 = toc

eigenvalues2 = [lambda2 mu2]

% we use iterative method (Krylov Schur) to compute just few eigenvalues with the smallest |mu|
tic
[lambda3,mu3] = twopareigs_ks(A1,B1,C1,A2,B2,C2,neig);
time3 = toc

eigenvalues3 = [lambda3 mu3]




