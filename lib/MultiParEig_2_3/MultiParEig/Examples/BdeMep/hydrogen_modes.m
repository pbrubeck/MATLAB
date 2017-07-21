function [S,E] = hydrogen_modes(R,N1,N2,b,opts)

%HYDROGEN_MODES  first four modes for hydrogen molecular ion H2+ in 2D
%
% [S,E] = HYDROGEN_MODES(R,N1,N2) returns separations constants S and 
% electronic energy E for the first four modes 1sg+, 2pu+, 2sg+, and 2pu- 
% from Table 1 and Table 2 (page 2201) in Patil,  Hydrogen molecular ion 
% and molecule in two dimensions, Journal of Chemical Physics 118, (2003) 2197-2205.
%
% Input:
%   - R : internuclear separation (should be >= 0.1)
%   - N1, N2 : number of points for the Laguerre collocation (N1) and Chebyshev collocation (N2) (default 60)
%   - b : scaling parameter for Laguerre collocation (default 1)
%   - opts : options
%
% Options in opts:
%   - fp_type: numeric type to use: 'single', 'double', or 'mp' (needs MCT),
%              default is the superior type of input data (x0,y0,z0)
%   - all options of twopareigs (init_double, refine_mp, ...)
%
% See also: HYDROGEN_MEP, DEMO_HYDROGEN, TWOPAREIGS

% References: 
%  - Patil, Hydrogen molecular ion and molecule in two dimensions, 
%    Journal of Chemical Physics 118, (2003) 2197-2205.
%  - B. Plestenjak, C. I. Gheorghiu, M. E. Hochstenbach: Spectral
%    collocation for multiparameter eigenvalue problems arising from separable
%    boundary value problems, J. Comp. Phys. 298 (2015) 585-601

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% BP 03.12.2016: modified to be precision-independent                 
% Last revision: 03.12.2016

narginchk(1, 5);

if nargin < 2, N1 = 60; end
if nargin < 3, N2 = 60; end
if nargin < 4, b = 1; end
if nargin < 5, opts = []; end
if isfield(opts,'fp_type') && is_numeric_type_supported(opts.fp_type)  
    class_t = opts.fp_type;   
else
    class_t = class(b);
end
opts.fp_type = class_t;

% Make sure all inputs are of the same numeric type.
if ~isa(b,class_t), b = numeric_t(b,class_t);  end;

S = zeros(1,4,class_t);
E = zeros(1,4,class_t);

% first we compute modes 1sg+, 2pu+, and 2sg+, where k = 0
[A1,B1,C1,A2,B2,C2,z1,z2] = hydrogen_mep(N1,N2,b,R,0,opts);
% we shift so that A1 and A2 are nonsingular and that the smallest absolute
% mu corresponds to smallest mu as a real number
A1 = A1+B1+10*C1; 
A2 = A2+B2+10*C2;
[lambda,mu,X1,X2] = twopareigs(A1,B1,C1,A2,B2,C2,6,opts); 
for j=1:6
    na = count_sign_changes(X1(:,j));
    nb = count_sign_changes(X2(:,j));
    if (na==0) && (nb==0) % mode 1sg+
        S(1,1) = lambda(j) - 1;
        E(1,1) = mu(j) - 10;
    end
    if (na==0) && (nb==1) % mode 2pu+
        S(1,2) = lambda(j) - 1;
        E(1,2) = mu(j) - 10;
    end
    if (na==1) && (nb==0) % mode 2sg+
        S(1,3) = lambda(j) - 1;
        E(1,3) = mu(j) - 10;
    end
end

% we compute mode 2pu-, where k = 1
[A1,B1,C1,A2,B2,C2,z1,z2] = hydrogen_mep(N1,N2,b,R,1,opts);
A1 = A1+B1+10*C1;
A2 = A2+B2+10*C2;
[lambda,mu,X1,X2] = twopareigs(A1,B1,C1,A2,B2,C2,3,opts);
for j=1:3
    na = count_sign_changes(X1(:,j));
    nb = count_sign_changes(X2(:,j));
    if (na==0) && (nb==0) % mode 2pu-
        S(1,4) = lambda(j) - 1;
        E(1,4) = mu(j) - 10;
    end
end