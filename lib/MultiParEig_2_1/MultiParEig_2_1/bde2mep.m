function [z,A,B,C,G,kd,rd] = bde2mep(a,b,p,q,r,s,t,bc,N)

%BDE2MEP Discretizes 2-parameter BDE using Chebyshev collocation
%
% [x,zt,A,B,C,G,kd,rd] = BDE2MEP(a,b,p,q,r,s,t,bc,N) discretizes a two-parameter DE
%  
%    p(x)y''(x) + q(x)y'(x) + r(x)y(x) = lambda s(x)y(x) + mu t(x)y(x)
%
% with boundary conditons 
% 
%    alpha1 y(a) + beta1 y'(a) = 0 for (alpha1,beta1)<>(0,0)
%    alpha2 y(b) + beta2 y'(b) = 0 for (alpha2,beta2)<>(0,0)
%
%    if alpha1=beta1=0 then b.c. is q(a)y'(a) + r(a)y(a) = lambda s(a)y(a) + mu t(a)y(a)
%    if alpha2=beta2=0 then b.c. is q(b)y'(b) + r(b)y(b) = lambda s(b)y(b) + mu t(b)y(b)
%
% into matrices A,B,C such that Ax = lambda Bx + mu Cx
%
% Input:
%   - a,b : interval [a,b]
%   - p,q,r,s,t : coeficient functions  - function handles or constants 
%   - bc : boundary conditions : bc = [alpha1 beta1; alpha2 beta2]
%   - N : number of Chebyshev points 
%
%  Output:
%   - z : vector of Chebyshev nodes scaled to [a,b] and ordered from b to a
%   - A,B,C : matrices of size from (N-2)x(N-2) to NxN (depends on b.c.)
%   - G : give-back matrix with boundary conditions (see Hoepffner)
%   - kd : kept degrees of freedom from 1:N
%   - rd : removed degrees of freedom from 1:N
%
% Package dmsuite by J.A.C. Weideman is required (or other chebdif method)
%
% See also: BDE3MEP, RECOVER_BC.

% References: 
% 1) J. Hoepffner, Implementation of boundary conditions, 
%    www.fukagata.mech.keio.ac.jp/~jerome/ (2007)
% 2) B. Plestenjak, C. I. Gheorghiu, M. E. Hochstenbach: Spectral
%    collocation for multiparameter eigenvalue problems arising from separable
%    boundary value problems, J. Comp. Phys. 298 (2015) 585-601

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision: 8.9.2015

[x,tmpD] = chebdif(N,2); % note that nodes in chebdif go from 1 to -1 
D0 = eye(N);
D1 = tmpD(:,:,1);
D2 = tmpD(:,:,2);

z = x/2*(b-a)+1/2*(a+b); % points from b to a (substituted from [-1,1]) 
sf = 2/(b-a); % substitution factor

if isa(p,'function_handle'), diagP = diag(p(z)); else diagP = p*eye(N); end
if isa(q,'function_handle'), diagQ = diag(q(z)); else diagQ = q*eye(N); end
if isa(r,'function_handle'), diagR = diag(r(z)); else diagR = r*eye(N); end
if isa(s,'function_handle'), diagS = diag(s(z)); else diagS = s*eye(N); end
if isa(t,'function_handle'), diagT = diag(t(z)); else diagT = t*eye(N); end

A = sf^2*diagP*D2 + sf*diagQ*D1 + diagR;
B = diagS;
C = diagT;

% boundary condition in a
if norm(bc(1,:))==0
    % we take boundary condition from the DE (note that this works only when p(a)=0)
    MBC2 = [];
    rd = [];
    kd = 1:N;
else
    MBC2 = bc(1,:)*[D0(N,:); sf*D1(N,:)]; % b.c. in point a 
    rd = N;
    kd = 1:N-1;
end

% boundary condition in b
if norm(bc(2,:))==0
    % we take boundary condition from the DE (note that this works only when p(b)=0)
    MBC1 = [];
else
    MBC1 = bc(2,:)*[D0(1,:); sf*D1(1,:)]; % b.c. in point b
    rd = [1 rd];
    kd = kd(2:end);
end

if isempty(rd)
    G = 0;
else
    CM = [MBC1; MBC2]; % constraint matrix
    G = -CM(:,rd)\CM(:,kd); % give-back matrix
    A = A(kd,kd) + A(kd,rd)*G;
    B = B(kd,kd);
    C = C(kd,kd);
end

