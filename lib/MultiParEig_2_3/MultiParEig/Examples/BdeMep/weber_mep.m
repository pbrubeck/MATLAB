function [A1,B1,C1,A2,B2,C2,z1,z2,G1,k1,r1,G2,k2,r2] = weber_mep(n1,n2,mode,x0,vi0,opts)

%WEBER_MEP  Discretizes 2-parameter Weber system obtained from Helmholtz equation 
% in parabolic cylindrical coordinates as a 2EP
%
% [A1,B1,C1,A2,B2,C2,z1,z2,G1,k1,r1,G2,k2,r2] = WEBER_MEP(n1,n2,mode,x0,vi0)
% discretizes system of differential equations 
% 
% M''(x) - (lambda + mu x^2)M(x) = 0,  0 < x < x0,
% N''(y) + (lambda - mu y^2)N(y) = 0,  -vi0 < y < vi0,
%
% Input:
%   - n1, n2: number of points for the first and the second equation
%   - mode : boundary value conditions (beside M'(x0)=N'(-v0)=N'(v0)=0
%            1 : M(0)=0  2 : M'(0)=0
%   - x0, vi0: end points
%   - opts : options
%
% Options in opts:
%   - fp_type: numeric type to use: 'single', 'double', or 'mp' (needs MCT),
%              default is the superior type of input data (x0,vi0)
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
% In mode 1 we search odd N with N(0)=0, in mode 2 N is even with N'(0)=0
%
% See also: DEMO_WEBER, BDE2MEP, LAME_MEP, MATHIEU_MEP, BESSELWAVE_MEP

% References: 
%  - Willatzen & Lew Yan Voon, Theory of acoustic eigenmodes in parabolic 
%    cylindrical enclosures, Journal of Sound and Vibration, (2005) 251--264.
%  - B. Plestenjak, C. I. Gheorghiu, M. E. Hochstenbach: Spectral
%    collocation for multiparameter eigenvalue problems arising from separable
%    boundary value problems, J. Comp. Phys. 298 (2015) 585-601

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% BP 03.12.2016: modified to be precision-independent                 
% Last revision: 03.12.2016

narginchk(5, 6);

if nargin < 6, opts = []; end
if isfield(opts,'fp_type') && is_numeric_type_supported(opts.fp_type)  
    class_t = opts.fp_type;   
else
    class_t = superiorfloat(x0,vi0);
end
opts.fp_type = class_t;

% Make sure all inputs are of the same numeric type.
if ~isa(x0,class_t),  x0 = numeric_t(x0,class_t);  end;
if ~isa(vi0,class_t), vi0 = numeric_t(vi0,class_t);  end;

if mode == 1
     bc1 = numeric_t('[1 0; 0 1]',class_t); 
     bc2 = numeric_t('[1 0; 0 1]',class_t); % N is odd
else
     bc1 = numeric_t('[0 1; 0 1]',class_t); 
     bc2 = numeric_t('[0 1; 0 1]',class_t); % N is even
end
 
a1 = numeric_t('0',class_t); 
b1 = x0;
p1 = numeric_t('1',class_t);
q1 = numeric_t('0',class_t);
r1 = numeric_t('0',class_t);
s1 = numeric_t('1',class_t);
t1 = @(x) x.^2;
[z1,A1,B1,C1,G1,k1,r1] = bde2mep(a1,b1,p1,q1,r1,s1,t1,bc1,n1,opts);

a2 = numeric_t('0',class_t); 
b2 = vi0;
p2 = numeric_t('1',class_t);
q2 = numeric_t('0',class_t);
r2 = numeric_t('0',class_t);
s2 = numeric_t('-1',class_t);
t2 = @(y) y.^2;
[z2,A2,B2,C2,G2,k2,r2] = bde2mep(a2,b2,p2,q2,r2,s2,t2,bc2,n2,opts);


