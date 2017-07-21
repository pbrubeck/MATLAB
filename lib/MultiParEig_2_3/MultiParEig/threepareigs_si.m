function [lambda,mu,eta,X1,X2,X3,flag] = threepareigs_si(A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3,neig,opts)

%THREEPAREIGS_SI   Subspace iteration for a three-parameter eigenvalue problem
% [lambda,mu,X1,X2,X3,flag] = THREEPAREIGS_SI(A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3,neig,opts)
% returns neig eigenvalues (lambda,mu,eta) with the smallest |eta| of the 
% two-parameter eigenvalue problem
%
% A1*x = lambda*B1*x + mu*C1*x + eta*D1*x
% A2*y = lambda*B2*y + mu*C2*y + eta*D2*y
% A3*z = lambda*B3*z + mu*C3*z + eta*D3*z
%
% using subspace iteration with Arnoldi expansion and restart based on 
% selected Ritz vectors on the generalized eigenvalue problem 
% Delta3*w = eta*Delta0*w, 
% where Delta0 and Delta3 are corresponding operator determinants
% Delta0 =   | B1 C1 D1; B2 C2 D2; B3 C3 D3 |
% Delta3 = - | B1 C1 A1; B2 C2 A2; B3 C3 A3 |
%
% Input:
%   - A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3 : matrices, A1, A2, and A3 have to be nonsingular
%   - neig : number of eigenvalues (6)
%   - opts : options (see below)
% 
% Output: 
%   - lambda, mu, eta : eigenvalue parts (eigenvalues are (lambda(j),mu(j),eta(j))
%   - X1, X2, X3 : components of decomposable right eigenvectors 
%     (eigenvector is kron(X1(:,j),kron(X2(:,j),X3(:,j))), such that
%       (A1-lambda(j)*B1-mu(j)*C1-eta(j)*D1)X1(:,j)=0
%       (A2-lambda(j)*B2-mu(j)*C2-eta(j)*D2)X2(:,j)=0
%       (A3-lambda(j)*B3-mu(j)*C3-eta(j)*D3)X3(:,j)=0
%   - flag : convergence (0), no convergence (1)
%
% Options are (default values in parenthesis):
%   - delta : absolute tolerance (1e-8)
%   - arnsteps : how many steps of Arnoldi we do in Arnoldi expansion (1)
%   - lowrank : after each iteration we take the leading lowrank ritz vectors (neig+1)
%   - window : how many Ritz values we consider for soft locking (2*neig), set to lowrank for no soft locking
%   - maxsteps : maximum number of outer steps
%   - showinfo : display full progress (2), just eigenvales found (1), nothing (0), default is 1
%   - harmonic : set to 1 to use harmonic instead of Ritz values (0) - use this for interior eigenvalue
%   - all options for twopareigs
%
% Please note that the method works only for really small number of wanted
% eigenvalues as the projected three-parameter problem can not be solved
% when its dimension is too large. 
% Safe settings are lowrank = window = 2, arnsteps = 1 
%
% See also: THREEPAREIG, THREEPAREIGS, THREEPAREIGS_JD, TWOPAREIGS_SI,
% DEMO_THREEPAREIGS_SI.

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% BP 06.09.2015 : no more seed to set random generator
% Last revision: 8.9.2015

if nargin<13, neig = 6; end
if nargin<14, opts = []; end

% options for the outer loop
if isfield(opts,'delta'),       delta = opts.delta;              else delta = 1e-8;           end
if isfield(opts,'arnsteps'),    arnsteps = opts.arnsteps;        else arnsteps = 1;           end
if isfield(opts,'maxsteps'),    maxsteps = opts.maxsteps;        else maxsteps = 20;          end
if isfield(opts,'showinfo'),    showinfo = opts.showinfo;        else showinfo = 1;           end
if isfield(opts,'harmonic'),    harmonic = opts.harmonic;        else harmonic = 0;           end
if isfield(opts,'lowrank'),     lowrank = opts.lowrank;          else lowrank = neig+1;       end
if isfield(opts,'window'),      window = opts.window;            else window = 2*neig;        end
% options for the inner solver twopareigs
if isfield(opts,'refine'),      refine = opts.refine;            else refine = 2;             end
if isfield(opts,'tol'),         tol = opts.tol;                  else tol = 1e-12;            end

sparseA = issparse(A1);

lambda = []; mu = []; eta = []; X1 = []; X2 = []; X3 = []; 

% Preparation for Arnoldi expansion and other computations
% If matrices are sparse we solve system in each step, otherwise we precompute the inverse.
if ~sparseA % we explicitly divide by A1 and A2 
   MB1 = A1\B1; MC1 = A1\C1; MD1 = A1\D1;
   MB2 = A2\B2; MC2 = A2\C2; MD2 = A2\D2;
   MB3 = A3\B3; MC3 = A3\C3; MD3 = A3\D3;
end

n1 = size(A1,1); 
n2 = size(A2,1);
n3 = size(A3,1);

% initial random subspaces
U1 = randn(n1,2); U1 = U1/norm(U1);
U2 = randn(n2,2); U2 = U2/norm(U2);
U3 = randn(n3,2); U3 = U3/norm(U3);

step = 0;
converged = 0;
flag = 0;
while (step < maxsteps) && (converged < neig)
   step = step + 1;

   % Block Arnoldi expansion
   if ~sparseA
        F1 = [MB1*U1 MC1*U1 MD1*U1];
        F2 = [MB2*U2 MC2*U2 MD2*U2];
        F3 = [MB3*U3 MC3*U3 MD3*U3];
        U1n = block_krylov_3p(MB1,MC1,F1,arnsteps);
        U2n = block_krylov_3p(MB2,MC2,F2,arnsteps);
        U3n = block_krylov_3p(MB3,MC3,F3,arnsteps);
   else
        F1 = A1\[B1*U1 C1*U1 D1*U1];
        F2 = A2\[B2*U2 C2*U2 D2*U2];
        F3 = A3\[B3*U3 C3*U3 D3*U3];
        U1n = block_krylov_3p(B1,C1,F1,arnsteps,A1);
        U2n = block_krylov_3p(B2,C2,F2,arnsteps,A2);
        U3n = block_krylov_3p(B3,C3,F3,arnsteps,A3);
   end
   % Projections of initial matrices on the subspace 
   PA1 = A1*U1n;  PB1 = B1*U1n;  PC1 = C1*U1n;  PD1 = D1*U1n;
   PA2 = A2*U2n;  PB2 = B2*U2n;  PC2 = C2*U2n;  PD2 = D2*U2n;
   PA3 = A3*U3n;  PB3 = B3*U3n;  PC3 = C3*U3n;  PD3 = D3*U3n;
   if ~harmonic
      prA1 = U1n'*PA1; prB1 = U1n'*PB1; prC1 = U1n'*PC1; prD1 = U1n'*PD1;
      prA2 = U2n'*PA2; prB2 = U2n'*PB2; prC2 = U2n'*PC2; prD2 = U2n'*PD2;
      prA3 = U3n'*PA3; prB3 = U3n'*PB3; prC3 = U3n'*PC3; prD3 = U3n'*PD3;
   else 
      U1left = PA1;
      U2left = PA2;
      U3left = PA3;
      prA1 = U1left'*PA1; prB1 = U1left'*PB1; prC1 = U1left'*PC1; prD1 = U1left'*PD1;
      prA2 = U2left'*PA2; prB2 = U2left'*PB2; prC2 = U2left'*PC2; prD2 = U2left'*PD2;
      prA3 = U3left'*PA3; prB3 = U3left'*PB3; prC3 = U3left'*PC3; prD3 = U3left'*PD3;
   end
   
   % Ritz (harmonic) values check
   tmpsize1 = size(prA1,1);
   tmpsize2 = size(prA2,1);
   tmpsize3 = size(prA3,1);
   newsize = min([tmpsize1*tmpsize2*tmpsize3 window]);
   [sigma,tau,rho,pXr,pYr,pZr] = threepareigs(prA1,prB1,prC1,prD1,prA2,prB2,prC2,prD2,prA3,prB3,prC3,prD3,window,opts);
   tmpvel = length(sigma);

   % compute the residuals and check if they are small enough
   converged = 0;
   for j = 1:tmpvel
       ost(j,1) = norm(PA1*pXr(:,j)-sigma(j)*PB1*pXr(:,j)-tau(j)*PC1*pXr(:,j)-rho(j)*PD1*pXr(:,j)); 
       ost(j,2) = norm(PA2*pYr(:,j)-sigma(j)*PB2*pYr(:,j)-tau(j)*PC2*pYr(:,j)-rho(j)*PD2*pYr(:,j));
       ost(j,3) = norm(PA3*pZr(:,j)-sigma(j)*PB3*pZr(:,j)-tau(j)*PC3*pZr(:,j)-rho(j)*PD3*pZr(:,j));
       if (ost(j,1)<delta) && (ost(j,2)<delta) && (ost(j,3)<delta)
           converged = converged + 1;
           mark(j) = 1;
       else
           mark(j) = 0;
       end
    end
    if showinfo == 2
        sweep = [sigma(1:tmpvel) tau(1:tmpvel) rho(1:tmpvel) ost] 
    end
    if converged < neig
       % we form new subspace out of current approx.
       tmpsize = min(lowrank,length(sigma));
       if showinfo
           fprintf('Step %3d, subspace size %3d x %3d x %3d, converged %2d\n',step,tmpsize1,tmpsize2,tmpsize3,converged)
       end
       V1next = [real(pXr(:,1:tmpsize)) imag(pXr(:,1:tmpsize))];
       V2next = [real(pYr(:,1:tmpsize)) imag(pYr(:,1:tmpsize))];
       V3next = [real(pZr(:,1:tmpsize)) imag(pZr(:,1:tmpsize))];
       U1 = U1n*orth(V1next);
       U2 = U2n*orth(V2next);
       U3 = U3n*orth(V3next);
    end       
end

if converged
    if showinfo
        fprintf('Step %3d, subspace size %3d x %3d x %3d, converged %2d\n',step,tmpsize1,tmpsize2,tmpsize3,converged)
    end
    ind = 1;
    for j=1:tmpvel
        if mark(j)
            lambda(ind,1) = sigma(j);
            mu(ind,1) = tau(j);
            eta(ind,1) = rho(j);
            X1(:,ind) = U1n*pXr(:,j);
            X2(:,ind) = U2n*pYr(:,j);
            X3(:,ind) = U3n*pZr(:,j);
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


