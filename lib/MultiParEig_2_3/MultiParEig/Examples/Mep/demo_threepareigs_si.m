%DEMO_THREEPAREIGS_SI    demo for threepareigs_si
%
% This example computes first 2 eigenmodes of an ellipsoid with semi-axes 1, 1.5, and 2
% using generalized subspace Arnoldi method for three-parameter eigenvalue problems
% The related Helmholtz equation separated in ellipsoidal coordinates is
% discretized with the Chebyshev collocation and solved as a
% three-parameter eigenvalue problem
%
% See also: THREEPAREIGS_SI, ELLIPSOIDWAVE_MEP

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

[lambda1,mu1,eta1,X1,X2,X3] = threepareigs(A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3,10);
[lambda1 mu1 eta1]

%% second solution: we use threepareigs_si to compute 2 eigenvalues for the 
% discretization with 400 nodes in each direction. The matrices are now so
% large that we can not use threepareigs anymore. As threepareigs_si calls
% threepareigs to solve the smaller projected three-parameter problems, we
% can not use neig>2 or arnsteps>1. If we want to compute more eigenvalues,
% we have to apply different shifts.

N = 400; % number of collocation nodes 
neig = 5;

[A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3,t1,t2,t3] = ellipsoidwave_mep(N,N,N,x0,y0,z0,rho,sigma,tau);

A1 = A1 + 5*B1 - shift*D1;
A2 = A2 + 5*B2 - shift*D2;
A3 = A3 + 5*B3 - shift*D3;

opts = [];
opts.showinfo = 1;
opts.maxsteps = 10;
opts.arnsteps = 1;
opts.lowrank = 2;
opts.window = 10;
opts.refine = 2;
opts.delta = 10*eps*max([norm(A1) norm(A2) norm(A3)]);
opts.usesparse = 0;
opts

tic
[lambda,mu,eta,X1,X2,X3] = threepareigs_si(A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3,neig,opts);
toc
[lambda mu eta]
