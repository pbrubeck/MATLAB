function [lambda,mu,X1,X2,flag] = twopareigs_ks(A1,B1,C1,A2,B2,C2,neig,opts)

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
%   - A1, B1, C1, A2, B2, C2 : matrices
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
%   - refine : TRQ refinement steps (2)
%   - epscluster : distance between eigenvalues from the same cluster (1e-4)
%   - fp_type: numeric type to use ('single', 'double', or 'mp' (needs MCT),
%     use only if you want to change the default - the superior type of input data
%   - init_double: (default 1) if the numerical type used is not double,
%     problem is solved first in double precision and then refined, if you
%     set to 0, then iterative methods will use multiprecision, which is
%     much slower
%   - refine_mp : TRQ refinement steps in multiprecision (10)
%
% For faster evaluation package lapack from MatlabCentral is required in
% Matlab below 2014a
%
% See also: TWOPAREIG, TWOPAREIGS, TWOPAREIGS_IRA, TWOPAREIGS_JD,
% TWOPAREIGS_SI, THREEPAREIGS.

% Reference: K. Meerbergen, B. Plestenjak: A Sylvester-Arnoldi type method 
% for the generalized eigenvalue problem with two-by-two operator determinants, 
% Numer. Linear Algebra Appl. 22 (2015) 1131-1146

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% BP 07.12.2016 : minor code improvements
% BP 03.12.2016 : fixed bug if last eigenvalues is from conjugate pair
% BP 02.12.2016 : modified to be precision-independent, support for complex matrices
% BP 26.08.2015 : change option w0 to v0, use sylvester in Matlab 2014a
% Last revision: 7.12.2016

narginchk(6, 8);

n1 = size(A1,1);
n2 = size(A2,1);
n = n1*n2;

if nargin<7, k = 6; else k = neig; end      % if not specified, 6 eigenvalues are returned 

% Analyse user supplied options, if any.
if nargin<8, opts = []; end
if isfield(opts,'fp_type') && is_numeric_type_supported(opts.fp_type)  
    class_t = opts.fp_type;   
else
    class_t = superiorfloat(A1,B1,C1,A2,B2,C2);
end
if isfield(opts,'init_double'), init_double = opts.init_double;  else init_double = 1; end
% class_ks is the class that is used in Krylov-Schur
if init_double
    class_ks = 'double';
else
    class_ks = class_t;
end
if isfield(opts,'divA'),        divA = opts.divA;                else divA = 1;          end
if isfield(opts,'v0'),          v0 = opts.v0;                    else v0 = randn(n,1,class_ks); end
if isfield(opts,'tol'),         tol = opts.tol;                  else tol = numeric_t('1e4*eps',class_ks); end
if isfield(opts,'maxit'),       maxit = opts.maxit;              else maxit = 100;       end
if isfield(opts,'p'),           p = opts.p;                      else p = 2*k;           end
if isfield(opts,'disp'),        display = opts.disp;             else display = 0;       end
if isfield(opts,'refine'),      refine = opts.refine;            else refine = 2;        end
if isfield(opts,'refine_mp'),   refine_mp = opts.refine_mp;      else refine_mp = 10;    end
if isfield(opts,'epscluster'),  epscluster = opts.epscluster;    else epscluster = numeric_t('1e-4',class_ks); end

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

if isreal(A1) && isreal(B1) && isreal(C1) && isreal(A2) && isreal(B2) && isreal(C2)
    realMEP = 1; 
else
    realMEP = 0;
end

if init_double
    MA1 = double(A1); MB1 = double(B1); MC1 = double(C1);
    MA2 = double(A2); MB2 = double(B2); MC2 = double(C2);
else
    MA1 = A1; MB1 = B1; MC1 = C1; 
    MA2 = A2; MB2 = B2; MC2 = C2;
end

if ~verLessThan('matlab', '8.3') || strcmpi(class_ks,'mp')
    uselapack = 0; % Matlab at least 2014a or MCT
elseif exist('lapack','file') 
    uselapack = 1; % old Matlab with lapack package
else
    uselapack = 2; % old Matlab without lapack package
end

if uselapack<2
    % Schur decomposition (real or complex) for Bartels-Stewart
    if divA
        [U1,R1] = schur(MA2\MB2);
        [U2,R2] = schur(-transpose(MB1)/transpose(MA1));
    else
        [U1,R1] = schur(MB2\MA2);
        [U2,R2] = schur(-transpose(MA1)/transpose(MB1));
    end
else
    % in old Matlab we use complex Schur even for real matrices to solve the Sylvester equaiton
    if divA
        [U1,R1] = schur(MA2\MB2,'complex');
        [U2,R2] = schur(-transpose(MB1)/transpose(MA1),'complex');
    else
        [U1,R1] = schur(MB2\MA2,'complex');
        [U2,R2] = schur(-transpose(MA1)/transpose(MB1),'complex');
    end
end

% y = multGamma2(x) solves linear system Delta2*y = Delta0*x 
% using Sylvester equation and method sylv
function y = multGamma2(x)
    MX = reshape(x,n2,n1); % MX = mat(x)
    FF = MC2*MX*transpose(MB1) - MB2*MX*transpose(MC1); % f = Delta0*w
    if divA
        y = reshape(sylv(U1,R1,U2,R2, -MA2\FF/transpose(MA1),uselapack),n1*n2,1); % wt = Delta2\Delta0*w
    else
        y = reshape(sylv(U1,R1,U2,R2, MB2\FF/transpose(MB1),uselapack),n1*n2,1); % wt = Delta2\Delta0*w
    end
end

% we start the Arnoldi subspace with initial vector v0 of size n
W = v0/norm(v0);
H = zeros(1,0,class_ks);
iter = 0;
err = Inf;
while (err>tol) && (iter < maxit)
    [W,H] = ArnoldiExpand(W,H,2*p); % we expand the Arnoldi subspace to maximum size 2*p
    WM = W(:,1:end-1); % all columns of the active subspace except the last one
    HM = H(1:end-1,:); % square Hessenberg matrix
    Hend = H(end,:); % last row in H
    Wend = W(:,end); % last column in W   
    
    [Q1,HM1] = schur(HM); % Schur form
    ritz = ordeig(HM1); % Ritz values
    [tilda,ord] = sort(abs(ritz));
    ordritz = ritz(ord); % ordered Ritz values
    invord(ord) = 1:length(ritz);
    for j = 2:length(ritz)
        if HM1(j,j-1)
            invord(ord(j)) = j-1;
        end
    end
    [Q2,T] = ordschur(Q1,HM1,invord); % reordered Schur form
        
    WM = WM*Q2;
    HM = T;

    if HM(p+1,p); 
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
        fprintf('Iteration %3d, error: %e\n',[iter err])
    end
end

if err<tol
    
    % we extract eigenvalues and eigenvector of MEP
    X1 = zeros(n1,k,class_ks);
    X2 = zeros(n2,k,class_ks);
    % last block in H might be 2x2
    if H(k+1,k)
        k = k + 1;
    end
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
        [G0,G1,tilda] = mult_delta(MA1,MB1,MC1,MA2,MB2,MC2,Z(:,j1:j2),Z(:,j1:j2));
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
           [refl,refu,refx1,refx2] = trqi(MA1,MB1,MC1,MA2,MB2,MC2,refx1,refx2,refine,eps);
           mu(j,1) = refu;
           lambda(j,1) = refl;
           X1(:,j) = refx1;
           X2(:,j) = refx2;
       else    
           % extraction of eigenvectors 
           X1(:,j) = min_sing_vec(MA1-lambda(j)*MB1-mu(j)*MC1);
           X2(:,j) = min_sing_vec(MA2-lambda(j)*MB2-mu(j)*MC2);
       end
    end
    flag = 0;
    if k>neig
        X1 = X1(:,1:neig);
        X2 = X2(:,1:neig);
        lambda = lambda(1:neig);
        mu = mu(1:neig);
    end
else
    % no convergence
    warning('twopareigs_ks has not converged')
    lambda = [];
    mu = [];
    X1 = [];
    X2 = [];
    flag = 1;
end

if init_double && ~strcmpi(class_t,'double')
    % TRQI refinement in higher precision
    RX1 = zeros(n1,k,class_t);
    RX2 = zeros(n2,k,class_t);
    Rlambda = zeros(k,1,class_t);
    Rmu = zeros(k,1,class_t);
    for j = 1:length(lambda)
        [Rlambda(j,1),Rmu(j,1),RX1(:,j),RX2(:,j)] = trqi(A1,B1,C1,A2,B2,C2,numeric_t(X1(:,j),class_t),numeric_t(X2(:,j),class_t),refine_mp,numeric_t('eps',class_t));
    end
    X1 = RX1;
    X2 = RX2;
    lambda = Rlambda;
    mu = Rmu;
end
        
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [W,H] = ArnoldiExpand(W,H,m)

% [W,H] = ArnoldiExpand(W,H,m) expands the active Krylov subspace 

k0 = size(W,2);
for j = k0:m;
    % we multiply by matrix using Sylvester equation and B-S algorithm
    % z = Inv(Delta2)*Delta0*W(:,j)
    z  = multGamma2(W(:,j));
    if (uselapack==2) && realSE % we take only the real part if a real system was solved using complex Schur 
        z = real(z);
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