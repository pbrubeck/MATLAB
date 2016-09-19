function [omega,XR,YR,z1,z2] = mathieu_modes(n1,n2,mode,alfa,beta,neig)

%MATHIEU_MODES  smallest eigenmodes for Mathieu's system
%
% [omega,XR,YR,z1,z2] = MATHIEU_MODES(n1,n2,mode,alfa,beta,neig) gives
% neig smallest eigenmodes for Mathieu's system of differential equations 
% 
% G''(x) + (lambda - 2mu  cos(2x))G(x) = 0,  0 < x < pi/2,
% F''(y) - (lambda - 2mu cosh(2y))F(y) = 0,  0 < y < xi0,
%
% where xi0 = acosh(alfa/sqrt(alfa^2-beta^2)) by converting in into a 
% two-parameter eigenvalue problem using Chebyshev collocation
%
% Input:
%   - n1, n2: number of points for the first and the second equation
%   - mode : boundary value conditions
%            1 : G'(0) = G'(pi/2)=0, F'(0) = F(xi0)=0
%            2 : G'(0) = G(pi/2)=0,  F'(0) = F(xi0)=0
%            3 : G(0)  = G(pi/2)=0,  F(0) = F(xi0)=0
%            4 : G(0)  = G'(pi/2)=0, F(0) = F(xi0)=0
%   - alfa, beta: big and small radius of the ellipse
%   - neig : number of smallest eigenmodes
%
% Output:
%   - omega : eigenfrequncies
%   - XR, YR: eigenmodes (in separated form)
%   - z1, z2 : collocation points for 1st and 2nd equation (including endpoints)
%
% See also: MATHIEU_MEP, TWOPAREIGS, DEMO_MATHIEU, ELLIPSE_EIGS 

% Reference: C. I. Gheorghiu, M. E. Hochstenbach, B. Plestenjak, J. Rommes: 
% Spectral collocation solutions to multiparameter Mathieu’s system, 
% Appl. Math. Comput. 218 (2012) 11990-12000.

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision: 8.9.2015

[A1,B1,C1,A2,B2,C2,z1,z2,G1,k1,r1,G2,k2,r2] = mathieu_mep(n1,n2,mode,alfa,beta);

if mode==1 % in first mode A1 is singular, we shift lambda to lambda - 5
    A1 = A1 + 5*B1;
    A2 = A2 + 5*B2;
end

h = sqrt(alfa^2-beta^2);
[lambda,mu,X,Y] = twopareigs(A1,B1,C1,A2,B2,C2,neig);
XR = recover_bc(X,G1,k1,r1);
YR = recover_bc(Y,G2,k2,r2);
omega = 2/h*sqrt(mu);

