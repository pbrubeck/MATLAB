function [lambda,mu,X1,X2,flag] = twopareigs_si(A1,B1,C1,A2,B2,C2,neig,opts)

%TWOPAREIGS_SI   Subspace iteration for two-parameter eigenvalue problem
%
% [lambda,mu,X1,X2,flag] = TWOPAREIGS_SI(A1,B1,C1,A2,B2,C2,neig,opts)
% returns neig eigenvalues (lambda,mu) with the smallest |mu| of the 
% two-parameter eigenvalue problem
%
% A1*x = lambda*B1*x + mu*C1*x
% A2*y = lambda*B2*y + mu*C2*y
%
% using subspace iteration with Arnoldi expansion and restart based on 
% selected Ritz vectors on the generalized eigenvalue problem 
% Delta2*z = lambda*Delta0*z, where z = kron(x,y) and
% Delta0 = kron(B1,C2) - kron(C1,B2)
% Delta2 = kron(B1,A2) - kron(A1,B2)
%
% Input:
%   - A1, B1, C1, A2, B2, C2 : matrices, A1 and A2 have to be nonsingular
%   - neig : number of eigenvalues (6)
%   - opts : options (see below)
% 
% Output: 
%   - lambda, mu : eigenvalue parts (eigenvalues are (lambda(j),mu(j))
%   - X1, X2 : components of decomposable right eigenvectors (eigenvector is kron(X1(:,j),X2(:,j)), such that
%       (A1-lambda(j)*B1-mu(j)*C1)X1(:,j)=0
%       (A2-lambda(j)*B2-mu(j)*C2)X2(:,j)=0
%   - flag : convergence (0), no convergence (1)
%
% Options are (default values in parenthesis):
%   - delta : absolute tolerance (1e-8)
%   - deltas : tolerance for soft locking  - if both residuals are below the tolerance, then the vector is soft locked (1e-4)
%   - arnsteps : how many steps of Arnoldi we do in Arnoldi expansion (1)
%   - lowrank : after each iteration we take the leading lowrank ritz vectors (neig+1)
%   - softlock : 0 or (1) take ritz vectors from place lowrank+1 to window if their residual is small
%   - window : how many Ritz values of the projected problem do we compute
%   - maxsteps : maximum number of outer steps
%   - showinfo : display full progress (2), just eigenvales found (1), nothing (0), default is 1
%   - harmonic : set to 1 to use harmonic instead of Ritz values (0) - use this for interior eigenvalue
%   - fp_type: numeric type to use ('single', 'double', or 'mp' (needs MCT),
%     use only if you want to change the default - the superior type of input data
%   - all options for twopareigs
%
% See also: THREEPAREIGS_SI, TWOPAREIG, TWOPAREIGS, TWOPAREIGS_IRA, TWOPAREIGS_KS,
% TWOPAREIGS_JD

% Reference: K. Meerbergen, B. Plestenjak: A Sylvester-Arnoldi type method 
% for the generalized eigenvalue problem with two-by-two operator determinants, 
% Report TW 653, Department of Computer Science, KU Leuven, 2014.

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% BP 02.12.2016 : modified to be precision-independent, support for complex matrices
% BP 06.09.2015 : no more seed to set random generator
% Last revision: 3.12.2016

narginchk(6, 8);

if nargin<7, neig = 6; end
% Analyse user supplied options, if any.
if nargin<8, opts = []; end
if isfield(opts,'fp_type') && is_numeric_type_supported(opts.fp_type)  
    class_t = opts.fp_type;   
else
    class_t = superiorfloat(A1,B1,C1,A2,B2,C2);
end

% options for the outer loop
if isfield(opts,'delta'),       delta = opts.delta;              else delta = numeric_t('1e8*eps',class_t);           end
if isfield(opts,'deltas'),      deltas = opts.deltas;            else deltas = numeric_t('1e4',class_t)*delta;        end
if isfield(opts,'arnsteps'),    arnsteps = opts.arnsteps;        else arnsteps = 1;           end
if isfield(opts,'maxsteps'),    maxsteps = opts.maxsteps;        else maxsteps = 20;          end
if isfield(opts,'showinfo'),    showinfo = opts.showinfo;        else showinfo = 1;           end
if isfield(opts,'harmonic'),    harmonic = opts.harmonic;        else harmonic = 0;           end
if isfield(opts,'lowrank'),     lowrank = opts.lowrank;          else lowrank = neig+1;       end
if isfield(opts,'window'),      window = opts.window;            else window = 2*neig;        end
if isfield(opts,'softlock'),    softlock = opts.softlock;        else softlock = 1;           end
% options for the inner solver twopareigs
if isfield(opts,'refine'),      refine = opts.refine;            else refine = 2;             end
if isfield(opts,'tol'),         tol = opts.tol;                  else tol = numeric_t('1e4*eps',class_t);            end

sparseA = issparse(A1);

lambda = numeric_t([],class_t); 
mu = numeric_t([],class_t); 
X1 = numeric_t([],class_t); 
X2 = numeric_t([],class_t);  

% Preparation for Arnoldi expansion and other computations
% If matrices are sparse we solve system in each step, otherwise we precompute the inverse.
if ~sparseA % we explicitly divide by A1 and A2 
   MC1 = A1\C1;
   MC2 = A2\C2;
   MB2 = A2\B2;
   MB1 =  A1\B1;
end

n1 = size(A1,1); 
n2 = size(A2,1);

% initial random subspaces
V = randn(n1,1); V = V/norm(V);
U = randn(n2,1); U = U/norm(U);

step = 0;
converged = 0;
flag = 0;
while (step < maxsteps) && (converged < neig)
   step = step + 1;

   % Block Arnoldi expansion
   if ~sparseA
        F = [MC2*U MB2*U];
        G = [MB1*V MC1*V];
        Un = block_krylov(MB2,F,arnsteps);
        Vn = block_krylov(MB1,G,arnsteps);
   else
        F = [A2\(C2*U) A2\(B2*U)];
        G = [A1\(B1*V) A1\(C1*V)];
        Un = block_krylov(B2,F,arnsteps,A2);
        Vn = block_krylov(B1,G,arnsteps,A1);
   end

   % Projections of initial matrices on the subspace 
   PA1 = A1*Vn;   PB1 = B1*Vn;   PC1 = C1*Vn;
   PA2 = A2*Un;   PB2 = B2*Un;   PC2 = C2*Un;
   if ~harmonic
      projA1 = Vn'*PA1;   projB1 = Vn'*PB1;   projC1 = Vn'*PC1;
      projA2 = Un'*PA2;   projB2 = Un'*PB2;   projC2 = Un'*PC2;
   else 
      Vleft = PA1;
      Uleft = PA2;
      projA1 = Vleft'*PA1;   projB1 = Vleft'*PB1;   projC1 = Vleft'*PC1;
      projA2 = Uleft'*PA2;   projB2 = Uleft'*PB2;   projC2 = Uleft'*PC2;
   end
   
   % Ritz (harmonic) values check
   tmpsize1 = size(projA1,1);
   tmpsize2 = size(projA2,1);
   newsize = min([tmpsize1*tmpsize2 window]);
   [sigma,tau,pXr,pYr] = twopareigs(projA1,projB1,projC1,projA2,projB2,projC2,newsize,opts);
   tmpvel = length(sigma);

   % compute the residuals and check if they are small enough
   converged = 0;
   for j = 1:tmpvel
       ost(j,1) = norm(PA1*pXr(:,j)-sigma(j)*(PB1*pXr(:,j))-tau(j)*(PC1*pXr(:,j))); 
       ost(j,2) = norm(PA2*pYr(:,j)-sigma(j)*(PB2*pYr(:,j))-tau(j)*(PC2*pYr(:,j)));
       if (ost(j,1)<delta) && (ost(j,2)<delta)
           converged = converged + 1;
           mark(j) = 1;
       else
           mark(j) = 0;
       end
    end
    if showinfo == 2
        sweep = [sigma(1:tmpvel) tau(1:tmpvel) ost] 
    end
    if converged < neig
       % we form new subspace out of current approx.
       Vnext = numeric_t([],class_t); 
       Unext = numeric_t([],class_t); 
       tmpsize = min(lowrank,tmpvel);
       locked = 0;
       if softlock
           % Soft locking (we keep candidates that would be missed otherwise)
           for j = tmpsize+1:tmpvel
               if (ost(j,1)<deltas) && (ost(j,2)<deltas)
                   Vnext = [Vnext real(pXr(:,j)) imag(pXr(:,j))];
                   Unext = [Unext real(pYr(:,j)) imag(pYr(:,j))];
                   locked = locked + 1;
               end
           end
       end
       if showinfo
          fprintf('Step %3d, subspace size %3d x %3d, converged %2d, soft locked %2d\n',...
             step,tmpsize1,tmpsize2,converged,locked)
       end
       Vnext = [Vnext real(pXr(:,1:tmpsize)) imag(pXr(:,1:tmpsize))];
       Unext = [Unext real(pYr(:,1:tmpsize)) imag(pYr(:,1:tmpsize))];
       V = Vn*orth(Vnext);
       U = Un*orth(Unext);
    end       
end

if converged
    if showinfo
       fprintf('Step %3d, subspace size %3d x %3d, converged %2d, soft locked %2d\n',...
           step,tmpsize1,tmpsize2,converged,0)
    end
    ind = 1;
    for j=1:tmpvel
        if mark(j)
            lambda(ind,1) = sigma(j);
            mu(ind,1)= tau(j);
            X1(:,ind) = Vn*pXr(:,j);
            X2(:,ind) = Un*pYr(:,j);
            ind = ind + 1;
            if ind > neig % we have enough eigenvalues
                break;
            end
        end
    end
end

if converged < neig
    fprintf('Method has not converged, returning %d eigenvalues from last step\n',converged);
    flag = 1;
end


