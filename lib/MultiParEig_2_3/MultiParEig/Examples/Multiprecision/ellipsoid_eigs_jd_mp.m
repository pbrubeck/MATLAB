function [omega,X1,X2,X3,xi1,xi2,xi3,lambda,mu,eta,omegad] = ellipsoid_eigs_jd_mp(x0,y0,z0,rho,sigma,tau,neig,n1,n2,n3)

%ELLIPSOID_EIGS_JD_MP  First eigenmodes of type (r,s,t) of a tri-axial
% ellipsoid using multiple precision
%
% [omega,X1,X2,X3,xi1,xi2,xi3,lambda,mu,eta,omegad] = ELLIPSOID_EIGS_JD_MP(x0,y0,z0,rho,sigma,tau,neig,n1,n2,n3,opts)
% returns first neig eigenmodes of type (rho,sigma,tau) of an ellipsoid
% with semi-axes x0 < y0 < z0 using Chebyshev collocation and the
% Jacobi-Davidson method. The initial results computed in double precision
% are refined in multiple precision.
%
% The solution of Helmholtz equation with Dirichlet b.c. in ellipsoidal coordinates is
% X1(xi1)*X2(xi2)*X3(xi3), where 0 < xi3 < b < xi2 < a < xi1 < z0,
% where a = sqrt(z0^2-x0^2) and b = sqrt(z0^2-y0^2).
%
% Input:
%   - x0 < y0 < z0: semi-axes of the tri-axial ellipsoid
%   - rho, sigma, tau : eigenmode type (0,0,0), (1,0,0), ..., (1,1,1)
%   - neig: number of wanted eigenvalues
%   - n1, n2, n3: number of Chebyshev collocation nodes for the 1st, 2nd, and 3rd equation
% 
% Output:
%   - omega: eigenfrequencies
%   - X1,X2,X3 : eigenwave functions in nodes xi1, xi2, xi3 
%   - xi1, xi2, xi3: nodes
%   - lambda, eta, mu: eigenvalues of the related three-parameter eigenvalue problem
%   - omegad: double precision eigenfrequencies before refinement
%
%
% This function requires Multiprecision Computing Toolbox for MATLAB, see
% http://www.advanpix.com/
%
% See also: ELLIPSOIDWAVE_MEP, DEMO_ELLIPSOIDWAVE, THREEPAREIGS

% References: 
%  - M. Willatzen and L. C. Lew Yan Voon, Numerical implementation of the ellipsoidal 
%    wave equation and application to ellipsoidal quantum dots, Comput. Phys. Commun. 171 (2005) 1-18.
%  - B. Plestenjak, C. I. Gheorghiu, M. E. Hochstenbach: Spectral
%    collocation for multiparameter eigenvalue problems arising from separable
%    boundary value problems, J. Comp. Phys. 298 (2015) 585-601

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision: 01.12.2016

is_numeric_type_supported('mp'); % check if MCT is installed

% First we solve the problem using double precision, then we refine the
% solution using TRQI and matrices obtained in higher precision
% ------------------------------------------------------------------------
a2 = z0^2-x0^2;
b2 = z0^2-y0^2;
c = a2/b2;

[A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3] = ellipsoidwave_mep(n1,n2,n3,double(x0),double(y0),double(z0),rho,sigma,tau);

shifteta = 0;
shiftlambda = -5;

A1s = A1 - shiftlambda*B1 - shifteta*D1;
A2s = A2 - shiftlambda*B2 - shifteta*D2;
A3s = A3 - shiftlambda*B3 - shifteta*D3;

optsJD.M1 = inv(A1s);
optsJD.M2 = inv(A2s);
optsJD.M3 = inv(A3s);
optsJD.target = [0 0 0];
optsJD.harmonic = 0;
optsJD.extraction = 'mineta';
optsJD.innersteps = 0;
optsJD.minsize = 3;
optsJD.maxsize = 6;
optsJD.maxmaxsize = 9;
optsJD.maxsteps = 1000;
optsJD.delta = eps*max([norm(A1s) norm(A2s) norm(A3s)]);
optsJD.switcheps = 1e6*optsJD.delta;
optsJD.showinfo = 2;
optsJD.forcereal = 1;

[lambda,mu,eta,X1,X2,X3] = threepareigs_jd(A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3,neig,optsJD);

eta = eta + shifteta;
lambda = lambda + shiftlambda;
etaold = eta;
omegad = sqrt(4/double(b2)*eta); 

% Refinement in higher precision using new input data in multiple precision
% Refinement is much faster as using multiple precision in Jacobi-Davidson
% ------------------------------------------------------------------------
% New matrices in mp format
opts.fp_type = 'mp';
[A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3,t1,t2,t3] = ellipsoidwave_mep(n1,n2,n3,x0,y0,z0,rho,sigma,tau,opts);
% convert initial approximations into mp
lambda = mp(lambda);
mu = mp(mu);
eta = mp(eta);
X1 = mp(X1);
X2 = mp(X2);
X3 = mp(X3);
for i = 1:length(lambda)
     [lambda(i),mu(i),eta(i),X1(:,i),X2(:,i),X3(:,i)] = trqi_3p(A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3,X1(:,i),X2(:,i),X3(:,i),5,mp('1e2*eps'));
end

X1 = [zeros(1,neig); X1]; % reconstruction of b.c. in x10

xi1 = sqrt(t1*b2);
xi2 = sqrt(t2*b2);
xi3 = sqrt(t3*b2);

X1 = X1.*(abs( (t1.^(rho/2)) .* ((t1-1).^(sigma/2)) .* ((t1-c).^(tau/2)) )*ones(1,neig));
X2 = X2.*(abs( (t2.^(rho/2)) .* ((t2-1).^(sigma/2)) .* ((t2-c).^(tau/2)) )*ones(1,neig));
X3 = X3.*(abs( (t3.^(rho/2)) .* ((t3-1).^(sigma/2)) .* ((t3-c).^(tau/2)) )*ones(1,neig));

omega = sqrt(4/b2*eta);