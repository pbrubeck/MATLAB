function [A1,B1,C1,A2,B2,C2,z1,z2,G1,k1,r1,G2,k2,r2] = besselwave_mep(n1,n2,p,xi0,eta0,BC,opts)

%BESSELWAVE_MEP  Discretizes Bessel wave equations as a two-parameter eigenvalue problem
%
% [A1,B1,C1,A2,B2,C2,z1,z2,G1,k1,r1,G2,k2,r2] = BESSELWAVE_MEP(n1,n2,p,xi0,eta0,BC)
% transforms system of Bessel vawe differential equations 
% 
% x^2 M''(x) + x M'(x) -p^2 M(x) + (lambda x^2 + mu x^4) M(x) = 0,  0 < x < xi0,
% y^2 N''(y) + y N'(y) -p^2 N(y) - (lambda y^2 + mu y^4) N(y) = 0,  0 < y < eta0,
%
% into a two-parameter eigenvalue problem using Chebyshev collocation
%
% Input:
%   - n1, n2: number of points for the first and the second equation
%   - p: parameter from the third equation  (integer >=0) 
%   - xi0, eta0: end points of intervals for the 1st and 2nd equation
%   - BC : boundary value conditions for M at xi0 and N at eta0 
%            [1 1] : M(xi0) = 0  and N(eta0) = 0  (Dirichlet & Dirichlet)
%            [2 2] : M'(xi0) = 0 and N'(eta0) = 0 (Neumann & Neumann)
%            [1 2] : M(xi0) = 0  and N'(eta0) = 0 (Dirichlet & Neumann)
%          in point 0 we always have Neumann b.c. (M'(0)=N'(0)=0)
%   - opts : options
%
% Options in opts:
%   - fp_type: numeric type to use: 'single', 'double', or 'mp' (needs MCT),
%              default is the superior type of input data (x0,y0,z0)
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
% See also: BDE2MEP, DEMO_BESSELWAVE1, DEMO_BESSELWAVE2, DEMO_BESSELWAVE3,
% DEMO_BESSELWAVE_FIGS

% References: 
%  - Lew Yan Voon & Willatzen, Helmholtz equation in parabolic rotational 
%    coordinates: application to wave problems in quantum mechanics and
%    acoustics, Math. Comp. Sim. 65 (2004) 337--349.
%  - B. Plestenjak, C. I. Gheorghiu, M. E. Hochstenbach: Spectral
%    collocation for multiparameter eigenvalue problems arising from separable
%    boundary value problems, J. Comp. Phys. 298 (2015) 585-601

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% BP 03.12.2016: modified to be precision-independent                 
% Last revision: 03.12.2016

narginchk(6, 7);

if nargin < 7, opts = []; end
if isfield(opts,'fp_type') && is_numeric_type_supported(opts.fp_type)  
    class_t = opts.fp_type;   
else
    class_t = superiorfloat(xi0,eta0);
end
opts.fp_type = class_t;

% Make sure all inputs are of the same numeric type.
if ~isa(xi0,class_t),  xi0 = numeric_t(xi0,class_t);  end;
if ~isa(eta0,class_t), eta0 = numeric_t(eta0,class_t);  end;

tmpBC = numeric_t('[0 0]',class_t); 
tmpBC(BC(1)) = 1; 
if p==0
    bc1 = [0 1; tmpBC];
else
    bc1 = [1 0; tmpBC];
end
a1 = numeric_t('0',class_t); 
b1 = xi0;
p1 = @(x) x.^2;
q1 = @(x) x;
r1 = -p^2; 
s1 = @(x) -x.^2; 
t1 = @(x) -x.^4;
[z1,A1,B1,C1,G1,k1,r1] = bde2mep(a1,b1,p1,q1,r1,s1,t1,bc1,n1,opts);

tmpBC = numeric_t('[0 0]',class_t);
tmpBC(BC(2)) = 1;
if p==0
    bc2 = [0 1; tmpBC];
else
    bc2 = [1 0; tmpBC];
end
a2 = numeric_t('0',class_t); 
b2 = eta0;
p2 = @(y)  y.^2;
q2 = @(y)  y;
r2 = -p^2; 
s2 = @(y)  y.^2; 
t2 = @(y) -y.^4;
[z2,A2,B2,C2,G2,k2,r2] = bde2mep(a2,b2,p2,q2,r2,s2,t2,bc2,n2,opts);



