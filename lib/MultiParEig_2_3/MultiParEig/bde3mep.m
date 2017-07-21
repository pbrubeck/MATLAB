function [z,A,B,C,D,G,kd,rd] = bde3mep(a,b,p,q,r,s,t,u,bc,N,opts)

%BDE3MEP Discretizes 3-parameter BDE using Chebyshev collocation
%
% [z,A,B,C,D,G,kd,rd] = BDE3MEP(a,b,p,q,r,s,t,u,bc,N,opts) discretizes a 3-parameter BDE
% 
%    p(x)y''(x) + q(x)y'(x) + r(x)y(x) = lambda s(x)y(x) + mu t(x)y(x) + eta u(x)y(x)
%
% with boundary conditons 
% 
%    alpha1 y(a) + beta1 y'(a) = 0
%    alpha2 y(b) + beta2 y'(b) = 0
%
%    if alpha1=beta1=0 then b.c. is q(a)y'(a) + r(a)y(a) = lambda s(a)y(a) + mu t(a)y(a) + eta u(a)y(a)
%    if alpha2=beta2=0 then b.c. is q(b)y'(b) + r(b)y(b) = lambda s(b)y(b) + mu t(b)y(b) + eta u(b)y(b)
%
% into matrices A,B,C,D such that Ax = lambda Bx + mu Cx + eta Dx
%
% Input:
%   - a,b : interval [a,b]
%   - p,q,r,s,t,u : coeficient functions  - function handles or constants 
%   - bc : boundary conditions : bc = [alpha1 beta1; alpha2 beta2]
%   - N : number of Chebyshev points 
%   - opts : options
%
% Options in opts:
%   - fp_type: numeric type to use: 'single', 'double', or 'mp' (needs MCT),
%              default is the superior type of input data
%
%  Output:
%   - z : vector of Chebyshev nodes scaled to [a,b] and ordered from b to a
%   - A,B,C,D : matrices of size from (N-2)x(N-2) to NxN (depends on b.c.)
%   - G : give-back matrix with boundary conditions (see Hoepffner)
%   - kd : kept degrees of freedom from 1:N
%   - rd : removed degrees of freedom from 1:N
%
% See also BDE2MEP, RECOVER_BC.

% References: 
% 1) J. Hoepffner, Implementation of boundary conditions, 
%    www.fukagata.mech.keio.ac.jp/~jerome/ (2007)
% 2) B. Plestenjak, C. I. Gheorghiu, M. E. Hochstenbach: Spectral
%    collocation for multiparameter eigenvalue problems arising from separable
%    boundary value problems, J. Comp. Phys. 298 (2015) 585-601

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% BP 29.11.2016: modified to be precision-independent                 
% Last revision: 29.11.2016

narginchk(10, 11);

if nargin < 11, opts = []; end
if isfield(opts,'fp_type') && is_numeric_type_supported(opts.fp_type)  
    class_t = opts.fp_type;   
else
    class_t = superiorfloat(a,b);
end

% Make sure all inputs are of the same numeric type.
if ~isa(a,class_t),  a = numeric_t(a,class_t);   end;
if ~isa(b,class_t),  b = numeric_t(b,class_t);   end;
if ~isa(bc,class_t), bc = numeric_t(bc,class_t); end;

[x,tmpD] = chebdif_mp(N,2,class_t); % attention - nodes in chebdif go from 1 to -1 
D0 = eye(N,class_t);
D1 = tmpD(:,:,1);
D2 = tmpD(:,:,2);

z = x/2*(b-a)+1/2*(a+b); % points (substituted from [-1,1]) 
sf = 2/(b-a); % substitution factor

if isa(p,'function_handle'), diagP = diag(p(z)); else diagP = numeric_t(p,class_t)*eye(N); end
if isa(q,'function_handle'), diagQ = diag(q(z)); else diagQ = numeric_t(q,class_t)*eye(N); end
if isa(r,'function_handle'), diagR = diag(r(z)); else diagR = numeric_t(r,class_t)*eye(N); end
if isa(s,'function_handle'), diagS = diag(s(z)); else diagS = numeric_t(s,class_t)*eye(N); end
if isa(t,'function_handle'), diagT = diag(t(z)); else diagT = numeric_t(t,class_t)*eye(N); end
if isa(u,'function_handle'), diagU = diag(u(z)); else diagU = numeric_t(u,class_t)*eye(N); end

A = sf^2*diagP*D2 + sf*diagQ*D1 + diagR;
B = diagS;
C = diagT;
D = diagU;

% boundary condition in a
if norm(bc(1,:))==0
    % we take boundary condition from the DE (note that this works only when p(a)=0)
    rd = numeric_t([],class_t);  
    kd = 1:N;
    MBC2 = numeric_t([],class_t);  
else
    MBC2 = bc(1,:)*[D0(N,:); sf*D1(N,:)]; % b.c. in point a 
    rd = N;
    kd = 1:N-1;
end

% boundary condition in b
if norm(bc(2,:))==0
    % we take boundary condition from the DE (note that this works only when p(b)=0)
    MBC1 = numeric_t([],class_t);  
else
    MBC1 = bc(2,:)*[D0(1,:); sf*D1(1,:)]; % b.c. in point b
    rd = [1 rd];
    kd = kd(2:end);
end

if isempty(rd)
    G = numeric_t('0',class_t);  
else
    CM = [MBC1; MBC2]; % constraint matrix
    G = -CM(:,rd)\CM(:,kd); % give-back matrix
    A = A(kd,kd) + A(kd,rd)*G;
    B = B(kd,kd);
    C = C(kd,kd);
    D = D(kd,kd);
end

