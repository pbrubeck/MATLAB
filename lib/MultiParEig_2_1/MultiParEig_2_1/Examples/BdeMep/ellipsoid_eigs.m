function [omega,X1,X2,X3,xi1,xi2,xi3,lambda,mu,eta] = ellipsoid_eigs(x0,y0,z0,rho,sigma,tau,neig,n1,n2,n3)

%ELLIPSOID_EIGS  First eigenmodes of type (r,s,t) of a tri-axial ellipsoid
%
% [omega,lambda,mu,eta,X1,X2,X3,xi1,xi2,xi3] = ELLIPSOID_EIGS(x0,y0,z0,rho,sigma,tau,neig,n1,n2,n3)
% returns first neig eigenmodes of type (rho,sigma,tau) of an ellipsoid
% with semi-axes x0 < y0 < z0
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

% Last revision: 8.9.2015

a2 = z0^2-x0^2;
b2 = z0^2-y0^2;
c = a2/b2;

[A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3,t1,t2,t3] = ellipsoidwave_mep(n1,n2,n3,x0,y0,z0,rho,sigma,tau);
[lambda,mu,eta,X1,X2,X3] = threepareigs(A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3,neig);
X1 = [zeros(1,neig); X1]; % reconstruction of b.c. in x10

xi1 = sqrt(t1*b2);
xi2 = sqrt(t2*b2);
xi3 = sqrt(t3*b2);

X1 = X1.*(abs( (t1.^(rho/2)) .* ((t1-1).^(sigma/2)) .* ((t1-c).^(tau/2)) )*ones(1,neig));
X2 = X2.*(abs( (t2.^(rho/2)) .* ((t2-1).^(sigma/2)) .* ((t2-c).^(tau/2)) )*ones(1,neig));
X3 = X3.*(abs( (t3.^(rho/2)) .* ((t3-1).^(sigma/2)) .* ((t3-c).^(tau/2)) )*ones(1,neig));

omega = sqrt(4/b2*eta);
