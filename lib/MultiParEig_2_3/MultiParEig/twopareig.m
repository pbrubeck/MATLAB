function [lambda,mu,XR,YR,XL,YL,report] = twopareig(A1,B1,C1,A2,B2,C2,opts)

%TWOPAREIG   Solve a two-parameter eigenvalue problem
%
% [lambda,mu,XR,YR,XL,YL] = TWOPAREIG(A1,B1,C1,A2,B2,C2,opts) returns
% eigenvalues and eigenvectors of the two-parameter eigenvalue problem
%
% A1 x = lambda B1 x + mu C1 x 
% A2 y = lambda B2 y + mu C2 y
%
% Input:
%   - A1, B1, C1, A2, B2, C2: matrices
%   - opts: options (see below)
%
% Output: 
%   - lambda, mu: eigenvalue parts (eigenvalues are (lambda(j),mu(j))
%   - XR, YR: components of decomposable right eigenvectors
%     (eigenvector is kron(XR(:,j),YR(:,j)) such that 
%       (A1-lambda(j)*B1-mu(j)*C1)*XR(:,j)=0
%       (A2-lambda(j)*B2-mu(j)*C2)*YR(:,j)=0
%   - XL, YL: components of decomposable left eigenvectors 
%     (eigenvector is kron(XL(:,j),YL(:,j)) such that 
%       (A1-lambda(j)*B1-mu(j)*C1)'*XL(:,j)=0
%       (A2-lambda(j)*B2-mu(j)*C2)'*YL(:,j)=0
%   - report: details of compression method in case of a singular problem 
%       rows are [mode m n r s(1) s(r) s(r+1) choice)], 
%       where mode is 1: CR step 1, 2: CR step 2, 3: RC step 1, 4. RC step 2
%       m,n: size of Delta matrices, r: rank
%       s(1), s(r), s(r+1): singular values
%       choice: how was rank determined (see NUMRANK for details)
%
% Operator determinants Delta0, Delta1, and Delta2 are used, where
% Delta0 = kron(C2, B1) - kron(B2, C1)
% Delta1 = kron(C2, A1) - kron(A2, C1)
% Delta2 = kron(A2, B1) - kron(B2, A1)
%
% Options in opts:
%   - singular (0): set to 1 for a singular problem, i.e., det(Delta0)=0
%   - epscluster (1e-6): relative distance between eigenvalues in a cluster
%   - fast (1): use fast algorithm (can fail for multiple eigenvalues) 
%     or slow algorithm (0) with clustering
%   - inviter (1): use inverse iteration for eigenvectors or slow svd (0)
%   - all options of auxiliary functions
%   - novectors (0): set to 1 when report is important and vectors are not
%   - maxgensize (0): maximal size of the regular part (0 means no bound)
%   - fp_type: numeric type to use ('single', 'double', or 'mp' (needs MCT),
%     use only if you want to change the default - the superior type of input data
%
% See also: TWOPAREIGS, TWOPAREIGS_JD, TWOPAREIGS_SI, THREEPAREIG, MULTIPAREIG

% Reference: M. E. Hochstenbach, T. Kosir, B. Plestenjak: A Jacobi-Davidson 
% type method for the two-parameter eigenvalue problem, SIAM J. Matrix Anal. 
% Appl. 26 (2005) 477-497

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% P. Holoborodko, Advanpix LLC.
% FreeBSD License, see LICENSE.txt

% BP 06.09.2015 : use extract_regular_part_np
% BP 24.08.2016 : double clustering in QZ (lambda first and then mu in blocks)
% BP 03.11.2016 : small speedup in computation od mu (inspired by Pavel Holoborodko's changes in multipareig)
% PH 22.11.2016 : modified to be precision-independent, i.e. be able to work with 
%                 any numeric type (whether it is built-in 'double'/'single' or custom like 'mp').
% BP 26.11.2016 : option fp_type
% PH 26.11.2016 : code simplifications, small fixes and clean-ups.

% Last revision: 26.11.2016

narginchk(6, 7);

% Analyse user supplied options, if any.
if nargin < 7, opts = []; end
if isfield(opts,'fp_type') && is_numeric_type_supported(opts.fp_type)  
    class_t = opts.fp_type;   
else
    class_t = superiorfloat(A1,B1,C1,A2,B2,C2);
end

if isfield(opts,'epscluster'), epscluster = opts.epscluster;   else, epscluster = numeric_t('1e-6',class_t);   end
if isfield(opts,'fast'),       fast       = opts.fast;         else, fast       = 1;                           end
if isfield(opts,'inviter'),    inviter    = opts.inviter;      else, inviter    = 1;                           end
if isfield(opts,'singular'),   singular   = opts.singular;     else, singular   = 0;                           end
if isfield(opts,'novectors'),  novectors  = opts.novectors;    else, novectors  = 0;                           end
if isfield(opts,'maxgensize'), maxgensize = opts.maxgensize;   else, maxgensize = 0;                           end

% Make sure all inputs are of the same numeric type.
if ~isa(A1,class_t), A1 = numeric_t(A1,class_t); end;
if ~isa(B1,class_t), B1 = numeric_t(B1,class_t); end;
if ~isa(C1,class_t), C1 = numeric_t(C1,class_t); end;
if ~isa(A2,class_t), A2 = numeric_t(A2,class_t); end;
if ~isa(B2,class_t), B2 = numeric_t(B2,class_t); end;
if ~isa(C2,class_t), C2 = numeric_t(C2,class_t); end;

% Default outputs
lambda = numeric_t([],class_t); 
mu     = numeric_t([],class_t); 
XR     = numeric_t([],class_t);  
YR     = numeric_t([],class_t);  
XL     = numeric_t([],class_t);  
YL     = numeric_t([],class_t); 
report = numeric_t([],class_t);

% Compute delta matrices 
[Delta0, Delta1, Delta2] = twopar_delta(A1,B1,C1,A2,B2,C2);

if singular
    [DeltaCell, report] = extract_regular_part_np({Delta0,Delta1,Delta2}, opts);
    Delta0 = DeltaCell{1}; Delta1 = DeltaCell{2}; Delta2 = DeltaCell{3}; 
end

% Quick return in degenerate cases
if ((size(Delta0,1) == 0) || (size(Delta0,1) ~= size(Delta0,2))) || ...  % no regular part was found
   ((maxgensize      > 0) && (size(Delta0,1)  > maxgensize))             % regular part is too large
     return
end

n = size(Delta0,1);
if fast
    tmp = Delta0\[Delta1 Delta2];
    Gamma1 = tmp(:,1:n);
    Gamma2 = tmp(:,n+1:end); 
    [Q1,D1] = schur(Gamma1,'complex');
    lambda = diag(D1);
    GQ = Gamma2*Q1;
    mu = zeros(n,1,class_t);
    for i = 1:n
        mu(i) = Q1(:,i)'*GQ(:,i);
    end
else
    [S0,S1,Q,Z,order,start,csize,tlambda] = clustered_qz(Delta0,Delta1,epscluster); %#ok<*ASGLU>
    if max(csize)==1
        lambda = tlambda;
        DZ = Delta2*Z;
        mu = zeros(n,1,class_t);
        for i = 1:n
            mu(i) = Q(i,:)*DZ(:,i);
        end
        mu = mu./diag(S0);
    else
        S2 = Q*Delta2*Z;
        for k = 1:length(start)
            partS0 = S0(start(k):(start(k)+csize(k)-1),start(k):(start(k)+csize(k)-1));
            partS1 = S1(start(k):(start(k)+csize(k)-1),start(k):(start(k)+csize(k)-1));
            partS2 = S2(start(k):(start(k)+csize(k)-1),start(k):(start(k)+csize(k)-1));
            % additional clustering in block partS2
            [partS0,partS2,Q2,Z2,order2,start2,csize2,partmu] = clustered_qz(partS0,partS2,epscluster);
            mu = [mu; partmu]; %#ok<*AGROW>
            partS1 = Q2*partS1*Z2;
            for j=1:length(start2)
                % computing lambda
                blockS0 = partS0(start2(j):(start2(j)+csize2(j)-1),start2(j):(start2(j)+csize2(j)-1));
                blockS1 = partS1(start2(j):(start2(j)+csize2(j)-1),start2(j):(start2(j)+csize2(j)-1));
                partlambda = eig(blockS1,blockS0);
                lambda = [lambda; partlambda];
            end
        end
    end
end

if (~novectors) && (nargout > 2) && all(isfinite(lambda))
    % extraction of eigenvectors (individually using inverse iteration or SVD)    
    
    % Generate initial vectors (in case of inverse iteration) only
    % once. This gives us a bit of speed-up, especially in extended precision
    % case, where multiple calls to randn might be noticeable.
    if inviter
        x10 = randn(size(A1,1),1,class_t);
        x20 = randn(size(A2,1),1,class_t);    
    else
        x10 = numeric_t([],class_t);
        x20 = numeric_t([],class_t);
    end;
   
    for k = 1:length(lambda)
        [xr,xl] = min_sing_vec(A1-lambda(k)*B1-mu(k)*C1,inviter,x10,x10);
        XR(:,k) = xr; XL(:,k) = xl;
        [yr,yl] = min_sing_vec(A2-lambda(k)*B2-mu(k)*C2,inviter,x20,x20);
        YR(:,k) = yr; YL(:,k) = yl;
    end   
end
end % twopareig