function [A1,B1,C1,A2,B2,C2,z1,z2,G1,k1,r1,G2,k2,r2] = lame_mep(xi,n1,n2,opts)

%LAME_MEP  Discretizes Lame system as a two-parameter eigenvalue problem
%
% [A1,B1,C1,A2,B2,C2,z1,z2,G1,k1,r1,G2,k2,r2] = LAME_MEP(xi,n1,n2)
% transforms Lame system of differential equations 
% 
% (1-k^2 cos^2(x)L''(x) + k^2 sin(x)cos(x)L'(x) + ((k^2 mu sin^2(x) + lambda)L(x) = 0,  0 < x < pi,
% (1-kp^2 cos^2(x)N''(x) + kp^2 sin(x)cos(x)N'(x) + ((kp^2 mu sin^2(x) - lambda)L(x) = 0,  0 < x < pi/2,
%
% where k = sin(|pi-xi|/2), k^2 + kp^2 = 1, into a two-parameter 
% eigenvalue problem using Chebyshev collocation
%
% Input:
%   - n1, n2: number of points for the first and the second equation
%   - xi : parameter from (0,2*pi)
%
% Boundary conditions are L(0)=L'(pi)=0 and 
% N'(0)=N'(pi/2)=0 for 0<xi<pi or N(0)=N'(pi/2)=0 for pi<xi<2*pi
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
% See also: BDE2MEP, DEMO_LAME

% Reference: 
%  - Morrison & Lewis, Charge singularity at the corner of a flat 
%    plate, SIAM J. Appl. Math. 31 (1976) 233--250.
%  - B. Plestenjak, C. I. Gheorghiu, M. E. Hochstenbach: Spectral
%    collocation for multiparameter eigenvalue problems arising from separable
%    boundary value problems, J. Comp. Phys. 298 (2015) 585-601

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% BP 03.12.2016: modified to be precision-independent                 
% Last revision: 03.12.2016

narginchk(3, 4);

if nargin < 4, opts = []; end
if isfield(opts,'fp_type') && is_numeric_type_supported(opts.fp_type)  
    class_t = opts.fp_type;   
else
    class_t = class(xi);
end
opts.fp_type = class_t;

k = sin(abs(numeric_t('pi',class_t)-xi)/2);
kp  = sqrt(1-k^2);

bc1 = numeric_t('[1 0; 0 1]',class_t);
a1 = numeric_t('0',class_t); 
b1 = numeric_t('pi',class_t);
p1 = @(x) (1-k^2*cos(x).^2);
q1 = @(x) k^2*sin(x).*cos(x);
r1 = numeric_t('0',class_t); 
s1 = numeric_t('-1',class_t); 
t1 = @(x) -k^2*sin(x).^2;
[z1,A1,B1,C1,G1,k1,r1] = bde2mep(a1,b1,p1,q1,r1,s1,t1,bc1,n1,opts);

if xi<numeric_t('pi',class_t) 
    bc2 = numeric_t('[0 1; 0 1]',class_t);
else
    bc2 = numeric_t('[1 0; 0 1]',class_t);
end

a2 = numeric_t('0',class_t); 
b2 = numeric_t('pi/2',class_t);
p2 = @(y) (1-kp^2*cos(y).^2);
q2 = @(y) kp^2*sin(y).*cos(y);
r2 = numeric_t('0',class_t); 
s2 = numeric_t('1',class_t); 
t2 = @(y) -kp^2*sin(y).^2;
[z2,A2,B2,C2,G2,k2,r2] = bde2mep(a2,b2,p2,q2,r2,s2,t2,bc2,n2,opts);