function [A1,B1,C1,A2,B2,C2,z1,z2,G1,k1,r1,G2,k2,r2] = mathieu_mep(n1,n2,mode,alfa,beta)

%MATHIEU_MEP   Discretizes Mathieu system as a two-parameter eigenvalue problem
%
% [A1,B1,C1,A2,B2,C2,z1,z2,G1,G2] = MATHIEU_MEP(n1,n2,mode,alfa,beta)
% transforms Mathieu system of differential equations 
% 
% G''(x) + (lambda - 2mu  cos(2x))G(x) = 0,  0 < x < pi/2,
% F''(y) - (lambda - 2mu cosh(2y))F(y) = 0,  0 < y < xi0,
%
% where xi0 = acosh(alfa/sqrt(alfa^2-beta^2))into a two-parameter 
% eigenvalue problem using Chebyshev collocation
%
% Input:
%   - n1, n2: number of points for the first and the second equation
%   - mode : boundary value conditions
%            1 : G'(0) = G'(pi/2)=0, F'(0) = F(xi0)=0
%            2 : G'(0) = G(pi/2)=0,  F'(0) = F(xi0)=0
%            3 : G(0)  = G(pi/2)=0,  F(0) = F(xi0)=0
%            4 : G(0)  = G'(pi/2)=0, F(0) = F(xi0)=0
%   - alfa, beta: big and small radius of the ellipse
%
% Output:
%   - A1,B1,C1 : (n1-2)x(n1-2) matrices for the first equation
%   - A2,B2,C2 : (n2-2)x(n2-2) matrices for the second equation
%   - z1, z2 : ni points for 1st and 2nd equation (including endpoints)
%   - G1 : 2 x (n1-2) give-back matrix with b.c. (see Hoepffner) for 1st equation
%   - G2 : 2 x (n2-2) give-back matrix with b.c. (see Hoepffner) for 2nd equation
%   - k1, r1 : kept indices and removed indices for the 1st equation
%   - k2, r2 : kept indices and removed indices for the 2nd equation
%
% See also: BDE2MEP, MATHIEU_MODES, DEMO_MATHIEU

% Reference: C. I. Gheorghiu, M. E. Hochstenbach, B. Plestenjak, J. Rommes: 
% Spectral collocation solutions to multiparameter Mathieu’s system, 
% Appl. Math. Comput. 218 (2012) 11990-12000.

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision: 11.09.2014

h = sqrt(alfa^2-beta^2);
xi0 = acosh(alfa/h);

switch mode
    case 1
        bc1 = [0 1; 0 1]; % G'(0) = G'(pi/2)=0
        bc2 = [0 1; 1 0]; % F'(0) = F(xi0)=0
    case 2
        bc1 = [0 1; 1 0]; % G'(0) = G(pi/2)=0
        bc2 = [0 1; 1 0]; % F'(0) = F(xi0)=0
    case 3
        bc1 = [1 0; 1 0]; % G(0) = G(pi/2)=0
        bc2 = [1 0; 1 0]; % F(0) = F(xi0)=0
    case 4
        bc1 = [1 0; 0 1]; % G(0) = G'(pi/2)=0
        bc2 = [1 0; 0 1]; % F(0) = F(xi0)=0
end

a1 = 0; 
b1 = pi/2;
p1 = 1;
q1 = 0;
r1 = 0;
s1 = -1;
t1 = @(x) 2*cos(2*x);
[z1,A1,B1,C1,G1,k1,r1] = bde2mep(a1,b1,p1,q1,r1,s1,t1,bc1,n1);

a2 = 0; 
b2 = xi0;
p2 = 1;
q2 = 0;
r2 = 0;
s2 = 1;
t2 = @(y) -2*cosh(2*y);
[z2,A2,B2,C2,G2,k2,r2] = bde2mep(a2,b2,p2,q2,r2,s2,t2,bc2,n2);


