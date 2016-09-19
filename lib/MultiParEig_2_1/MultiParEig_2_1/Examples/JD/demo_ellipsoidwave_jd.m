%DEMO_ELLIPSOIDWAVE_JD     demo for the Jacobi-Davidson method for 
% three-parameter eigenvalue problems
%
% This example computes first eigenmodes of an ellipsoid with semi-axes 1, 1.5, and 2
% using Jacobi-Davidson method for three-parameter eigenvalue problems
%
% The related Helmholtz equation separated in ellipsoidal coordinates is
% discretized with the Chebyshev collocation and solved as a
% three-parameter eigenvalue problem
%
% See also: ELLIPSOIDWAVE_MEP, THREEPAREIGS_JD, THREEPAREIGS

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision: 8.9.2015

x0 = 1;
y0 = 1.5;
z0 = 2;

% set rho, sigma and tau to 0 or 1 to find eigenmodes of type (rho,sigma,tau) 
rho = 0;
sigma = 0;
tau = 0;

shift = 0; % we are looking for eigenvalues such that [eta-shift| is minimal

%% first solution: we use threepareigs to compute 10 eigenvalues for the 
% discretization with 10 nodes in each direction

N = 15; % number of collocation nodes 

[A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3,t1,t2,t3] = ellipsoidwave_mep(N,N,N,x0,y0,z0,rho,sigma,tau);

% we apply shift and substitution in lambda to make A1, A2 and A3 nonsingular
A1 = A1 + 5*B1 - shift*D1;
A2 = A2 + 5*B2 - shift*D2;
A3 = A3 + 5*B3 - shift*D3;

tic
[lambda1,mu1,eta1,X1,X2,X3] = threepareigs(A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3,20);
toc
[lambda1 mu1 eta1]

%% second solution: we use threepareigs_jd to compute 20 eigenvalues for the 
% discretization with 400 nodes in each direction. The matrices are now so
% large that we can not use threepareigs anymore. As threepareigs_jd calls
% threepareig to solve the smaller projected three-parameter problems, we
% can not too large maxsize with window=0. 

N = 400; % number of collocation nodes 
neig = 20;

[A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3,t1,t2,t3] = ellipsoidwave_mep(N,N,N,x0,y0,z0,rho,sigma,tau);

A1 = A1 + 5*B1 - shift*D1;
A2 = A2 + 5*B2 - shift*D2;
A3 = A3 + 5*B3 - shift*D3;

opts = [];
opts.M1 = inv(A1);
opts.M2 = inv(A2);
opts.M3 = inv(A3);
opts.target = [0 0 0];
opts.harmonic = 1;
opts.extraction = 'mineta';
opts.innersteps = 8;
opts.minsize = 3;
opts.maxsize = 7;
opts.maxsteps = 500;
opts.delta = 5*eps*max(norm(A1),norm(A2));
opts.window = 0;
opts.showinfo = 2;
opts.forcereal = 1

tic
[lambda,mu,eta,X1,X2,X3] = threepareigs_jd(A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3,neig,opts);
toc

[lambda mu eta]
