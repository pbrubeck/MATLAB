function [lambda,mu,X1,X2,flag] = twopareigs_ks(A1,B1,C1,A2,B2,C2,k,opts)

%TWOPAREIGS_KS   Krylov-Schur for two-parameter eigenvalue problem
%
% [lambda,mu,X1,X2,flag] = TWOPAREIGS_KS((A1,B1,C1,A2,B2,C2,k,OPTS) 
% returns the eigenvales and eigenvectors with the smallest k eigenvalues 
% mu of the two-parameter eigenvalue problem with real matrices
%
% A1*x = lambda*B1*x + mu*C1*x
% A2*y = lambda*B2*y + mu*C2*y
%
% using Krylov-Schur algorithm on the generalized eigenvalue
% problem Delta2*z = lambda*Delta0*z, where z = kron(x,y) and
% Delta0 = kron(B1,C2) - kron(C1,B2)
% Delta2 = kron(B1,A2) - kron(A1,B2)
% 
% When building Krylov subspace, Bartels-Stewart method is used to solve 
% the Sylvester equation related to linear system Delta2*w = Delta0*z
%
% Input:
%   - A1, B1, C1, A2, B2, C2 : real matrices
%   - k : number of eigenvalues(6)
%   - opts: options (see below)
%
% Output: 
%   - lambda, mu: eigenvalue parts (eigenvalues are (lambda(j),mu(j))
%   - X1, X2: components of decomposable right eigenvectors (eigenvector is kron(X1(:,j),X2(:,j)) such that 
%       (A1-lambda(j)*B1-mu(j)*C1)X1(:,j)=0
%       (A2-lambda(j)*B2-mu(j)*C1)X2(:,j)=0
%   - flag: convergence (0), no convergence (1)
% 
% Possible options in opts:
%   - divA : (default 1) : divide by A1 and A2, or by B1 and B2 (0)
%   - v0 : initial vector (randn(n,1))
%   - tol : absolute tolerance (1e-12)
%   - maxit : maximum number of iterations (100)
%   - p : minimal size of search space (2k), maximal size is 2p
%   - disp : display progress - 1 or not (0)
%   - refine : TRQ refinement steps (1)
%   - epscluster : distance between eigenvalues from the same cluster (1e-4)
%
% For faster evaluation package lapack from MatlabCentral is required
%
% See also: TWOPAREIG, TWOPAREIGS, TWOPAREIGS_IRA, TWOPAREIGS_JD,
% TWOPAREIGS_SI, THREEPAREIGS.

% Reference: K. Meerbergen, B. Plestenjak: A Sylvester-Arnoldi type method 
% for the generalized eigenvalue problem with two-by-two operator determinants, 
% Numer. Linear Algebra Appl. 22 (2015) 1131-1146

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% BP 26.08.2015 : change option w0 to v0, use sylvester in Matlab 2014a
% Last revision: 8.9.2015

n1 = size(A1,1);
n2 = size(A2,1);
n = n1*n2;

if nargin<7, k=6; end      % if not specified, 6 eigenvalues are returned 
if nargin<8, opts = []; end
if isfield(opts,'divA'),       divA = opts.divA;              else divA = 1;          end
if isfield(opts,'v0'),         v0 = opts.v0;                  else v0 = randn(n,1);   end
if isfield(opts,'tol'),        tol = opts.tol;                else tol = 1e-12;       end
if isfield(opts,'maxit'),      maxit = opts.maxit;            else maxit = 100;       end
if isfield(opts,'p'),          p = opts.p;                    else p = 2*k;           end
if isfield(opts,'disp'),       display = opts.disp;           else display = 0;       end
if isfield(opts,'refine'),     refine = opts.refine;          else refine = 2;        end
if isfield(opts,'epscluster'), epscluster = opts.epscluster;  else epscluster = 1e-4; end

if p<k+1, p = k+1; end % search subspace must be at least p+1

if display
    fprintf('\n MEPKrylovSchur parameters:\n')
    fprintf('  - number of wanted eigenvalues: %d\n',k)
    fprintf('  - subspace size from %d to %d\n',p,2*p)
    fprintf('  - absolute tolerance: %e\n',tol)
    fprintf('  - maximum number of iterations: %d\n',maxit)
    fprintf('  - TRQ refinement steps: %d\n',refine)
    fprintf('  - epscluster: %e\n',epscluster)
end

if exist('lapack','file') || ~verLessThan('matlab', '8.3')
    uselapack = 1;
else
    uselapack = 0;
end

if uselapack
    % Schur decomposition (real) for Bartels-Stewart
    if divA
        [U1,R1] = schur(A2\B2);
        [U2,R2] = schur(-transpose(B1)/transpose(A1));
    else
        [U1,R1] = schur(B2\A2);
        [U2,R2] = schur(-transpose(A1)/transpose(B1));
    end
else
    if divA
        [U1,R1] = schur(A2\B2,'complex');
        [U2,R2] = schur(-transpose(B1)/transpose(A1),'complex');
    else
        [U1,R1] = schur(B2\A2,'complex');
        [U2,R2] = schur(-transpose(A1)/transpose(B1),'complex');
    end
end

% y = multGamma2(x) solves linear system Delta2*y = Delta0*x 
% using Sylvester equation and method SylvBSUTReal
function y = multGamma2Real(x)
    MX = reshape(x,n2,n1); % MX = mat(x)
    FF = C2*MX*transpose(B1) - B2*MX*transpose(C1); % f = Delta0*w
    if divA
        y = reshape(sylvreal(U1,R1,U2,R2, -A2\FF/transpose(A1)),n1*n2,1); % wt = Delta2\Delta0*w
    else
        y = reshape(sylvreal(U1,R1,U2,R2, B2\FF/transpose(B1)),n1*n2,1); % wt = Delta2\Delta0*w
    end
end

% y = multGamma2(x) solves linear system Delta2*y = Delta0*x 
% using Sylvester equation and method SylvBSUT and takes real part
function y = multGamma2(x)
    MX = reshape(x,n2,n1); % MX = mat(x)
    FF = C2*MX*transpose(B1) - B2*MX*transpose(C1); % f = Delta0*w
    if divA
        y = reshape(sylv(U1,R1,U2,R2, -A2\FF/transpose(A1),0),n1*n2,1); % wt = Delta2\Delta0*w
    else
        y = reshape(sylv(U1,R1,U2,R2, B2\FF/transpose(B1),0),n1*n2,1); % wt = Delta2\Delta0*w
    end
end

% we start the Arnoldi subspace with initial vector v0 of size n
W = v0/norm(v0);
H = zeros(1,0);
iter = 0;
err = Inf;
while (err>tol) && (iter < maxit)
    [W,H] = ArnoldiExpand(W,H,2*p,uselapack); % we expand the Arnoldi subspace to maximum size 2*p
    WM = W(:,1:end-1); % all columns of the active subspace except the last one
    HM = H(1:end-1,:); % square Hessenberg matrix
    Hend = H(end,:); % last row in H
    Wend = W(:,end); % last column in W
    
    [Q1,HM1] = schur(HM); % Schur form
    ritz = ordeig(HM1); % Ritz values
    [tilda,ord]=sort(abs(ritz));
    ordritz = ritz(ord); % ordered Ritz values
    invord = []; 
    for i=1:length(ritz); 
        invord(ord(i)) = i; 
    end
    for j=2:length(ritz)
        if abs(ordritz(j)-conj(ordritz(j-1)))/abs(ordritz(j))<1e-10 % conjugate pair
            invord(ord(j)) = j-1;
        end
    end
    [Q2,T] = ordschur(Q1,HM1,invord); % reordered Schur form
        
    WM = WM*Q2;
    HM = T;

    if abs(ordritz(p)-conj(ordritz(p+1)))/abs(ordritz(p))<1e-10 
        meja = p+1; % keep conjugate complex pairs together
    else
        meja = p;
    end
    
    H = HM(1:meja,1:meja);
    W = [WM(:,1:meja) Wend];
    ztilda = Hend*Q2(:,1:meja);
    H(meja + 1,:) = ztilda;
    err = norm(H(meja+1,1:k));
    iter = iter + 1;
    if display
        disp(sprintf('Iteration %3d, error: %e',[iter err]))
    end
end

if err<tol
    
    % we extract eigenvalues and eigenvector of MEP
    X1 = zeros(n1,k);
    X2 = zeros(n2,k);
    [X,D] = eig(H(1:k,1:k));
    % clustering of eigenvalues
    [order,clstart,clsize,newmu] = clusters(diag(D),epscluster);
    Z = W(:,1:k)*X(:,order); % eigenvectors sorted in clusters
    for c = 1:length(clstart)
        % for each cluster we compute matrices Gi = W'*Deltai*W
        % that give as lambda and correct eigenvectors
        % in cluster of size 1 this is just Rayleigh quotient
        j1 = clstart(c);
        j2 = clstart(c)+clsize(c)-1;
        [G0,G1,tilda] = mult_delta(A1,B1,C1,A2,B2,C2,Z(:,j1:j2),Z(:,j1:j2));
        [tmpX,tmpD] = eig(G1,G0);
        Z(:,j1:j2) = Z(:,j1:j2)*tmpX;
        lambda(j1:j2,1) = diag(tmpD);
    end
    mu = 1./newmu;
    for j = 1:k
       mz = reshape(Z(:,j),n2,n1);
       if refine
           [tmpU,tilda,tmpV] = svd(mz);
           refx1 = conj(tmpV(:,1));
           refx2 = tmpU(:,1);
           [refl,refu,refx1,refx2] = trqi(A1,B1,C1,A2,B2,C2,refx1,refx2,refine,eps);
           mu(j,1) = refu;
           lambda(j,1) = refl;
           X1(:,j) = refx1;
           X2(:,j) = refx2;
       else    
           % extraction of eigenvectors 
           X1(:,j) = min_sing_vec(A1-lambda(j)*B1-mu(j)*C1);
           X2(:,j) = min_sing_vec(A2-lambda(j)*B2-mu(j)*C2);
       end
    end
    flag = 0;
else
    % no convergence
    warning('twopareigs_ks has not converged')
    lambda = [];
    mu = [];
    X1 = [];
    X2 = [];
    flag = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [W,H] = ArnoldiExpand(W,H,m,uselapack)

% [W,H] = ArnoldiExpand(W,H,m) expands the active Krylov subspace 

k0 = size(W,2);
for j = k0:m;
    % we multiply by matrix using Sylvester equation and B-S algorithm
    % z = Inv(Delta2)*Delta0*W(:,j)
    if uselapack
        z = multGamma2Real(W(:,j)); 
    else
        z = real(multGamma2(W(:,j))); 
    end        
    t0 = norm(z);
    for i = 1:j
        H(i,j) = W(:,i)'*z;
        z = z - H(i,j)*W(:,i);
    end
    t = norm(z);
    if t < 0.7*t0 % reorthogonalization 
        for i = 1:j
            h2 = W(:,i)'*z;
            H(i,j) = H(i,j) + h2;
            z = z - h2*W(:,i);
        end
        t = norm(z);
    end 
    if t == 0
        return
    end
    H(j+1,j) = t;
    W(:,j+1) = z/H(j+1,j);
end

end % ArnoldiExpand

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end % twopareigs_ks