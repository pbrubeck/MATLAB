%DEMO_TWOPAREIGS_SI    demo for method twopareigs_si
%
% This example computes first 10 eigenmodes for the Mathieu two-parameter eigenvalue problem
% using subspace Arnoldi method for two-parameter eigenvalue problems
%
% See also: TWOPAREIGS_SI, TWOPAREIGS_IRAm TWOPAREIGS_KS, MATHIEU_MEP

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision: 8.9.2015

n = 200;   % size of the matrices
neig = 10; % number of wanted eigenvalues
[A1,B1,C1,A2,B2,C2] = mathieu_mep(n,n,2,2,1); 

%% First version uses implicitly restarted Arnoldi with full vectors and
% Bartels-Stewart method for the Sylvester equation
tic; [lambda1,mu1,X1,Y1] = twopareigs_ira(A1,B1,C1,A2,B2,C2,neig); t1=toc
mu1

%% Second version uses Krylov-Schur method with full vectors and
% Bartels-Stewart method for the Sylvester equation
% MEPKrylovSchur (10 eigenvalues = 2.2s)
tic; [lambda2,mu2,X2,Y2] = twopareigs_ks(A1,B1,C1,A2,B2,C2,neig); t2=toc
mu2

%% Third version is subspace Arnoldi, where we don't use full vectors
% and solve the projected small problem in each step
opts = [];
opts.lowrank = neig;
opts.window = neig;
opts.maxsteps = 10;
opts.arnsteps = 1;
opts.delta = eps*max(norm(A1),norm(A2));
opts.showinfo = 1;
opts.softlock = 0;
tic; [lambda3,mu3,X3,Y3] = twopareigs_si(A1,B1,C1,A2,B2,C2,neig,opts); t3=toc
mu3
