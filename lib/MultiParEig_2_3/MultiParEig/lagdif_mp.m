function [x, DM] = lagdif_mp(N, M, b, class_t)

%  The function [x, DM] = lagdif(N, M, b, class_t) computes the
%  differentiation matrices D1, D2, ..., DM on Laguerre points of class_t.
%
%  Input:
%  N:    Number of points, i.e., order of differentiation matrices (integer).
%  M:    Number of derivatives required (integer).
%  b:    Scaling parameter (real, positive).
%  class_t:  numeric type to use (optional), default is 'double'
%
%  Note:  0 < M < N-1.
%
%  Output:
%  x:    Vector of nodes (zeros of Laguerre polynomial of degree N-1,
%        plus x = 0), all scaled by the parameter b.
%  DM:   DM(1:N,1:N,l) contains l-th derivative matrix, l=1..M.
%
% This function was taken from Matlab package DMSUITE by J.A.C Weideman
% and modified to be precision-independent, i.e. be able to work with any 
% numeric type (e.g. built-in 'double'/'single' or custom like 'mp').

%  J.A.C. Weideman, S.C. Reddy 1998.

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% BP 04.12.2016: modified to be precision-independent
% Last revision: 04.12.2016
 
if nargin<4
    class_t = class(b);
end

x = lagroots(N-1,class_t);            % Compute Laguerre roots
x = [0; x];                           % Add a node at x=0 to facilitate
                                      % the implementation of BCs.

alpha = exp(-x/2);                    % Compute weights.

for ell = 1:M                             % Set up beta matrix s.t.
 beta(ell,:) = numeric_t('-0.5',class_t)^ell*ones(size(x')); % beta(ell,j) is (l^th derivative
end                                       % alpha(x))/alpha(x)
                                          % evaluated at x = x(j).

DM = poldif(x, alpha, beta,class_t);          % Compute differentiation matrix (b=1).

x = x/b;                              % Scale nodes by the factor b.

for ell = 1:M                           % Adjust for b not equal to 1.
 DM(:,:,ell) = (b^ell)*DM(:,:,ell);
end

%------------------------------------------------------------------------
function r = lagroots(N,class_t);

%  The function r = lagroots(N) computes the roots of the 
%  Laguerre polynomial of degree N.

% This function was taken from Matlab package DMSUITE by J.A.C Weideman
% and modified to be precision-independent, i.e. be able to work with any 
% numeric type (e.g. built-in 'double'/'single' or custom like 'mp').

%  J.A.C. Weideman, S.C. Reddy 1998.

J = diag([1:2:2*N-1])-diag([1:N-1],1)-diag([1:N-1],-1);  % Jacobi matrix
J = numeric_t(J,class_t); % we can do this as J has integer entries
r = sort(eig(J));                                % Compute eigenvalues

%------------------------------------------------------------------------
function DM = poldif(x, malpha, B, class_t)

%  The function DM =  poldif(x, maplha, B) computes the
%  differentiation matrices D1, D2, ..., DM on arbitrary nodes.
%
%  The function is called with either two or three input arguments.
%  If two input arguments are supplied, the weight function is assumed 
%  to be constant.   If three arguments are supplied, the weights should 
%  be defined as the second and third arguments.
%
%  Input (constant weight):
%
%  x:        Vector of N distinct nodes.
%  malpha:   M, the number of derivatives required (integer).
%  B:        Omitted.
%
%  Note:     0 < M < N-1.
%
%  Input (non-constant weight):
%
%  x:        Vector of N distinct nodes.
%  malpha:   Vector of weight values alpha(x), evaluated at x = x(k).
%  B:        Matrix of size M x N,  where M is the highest 
%            derivative required.  It should contain the quantities 
%            B(ell,j) = beta(ell,j) = (ell-th derivative
%            of alpha(x))/alpha(x),   evaluated at x = x(j).
%
%  Output:
%  DM:       DM(1:N,1:N,ell) contains ell-th derivative matrix, ell=1..M.

% This function was taken from Matlab package DMSUITE by J.A.C Weideman
% and modified to be precision-independent, i.e. be able to work with any 
% numeric type (e.g. built-in 'double'/'single' or custom like 'mp').

%  J.A.C. Weideman, S.C. Reddy 1998

       N = length(x);                      
       x = x(:);                     % Make sure x is a column vector 

   alpha = malpha(:);                % Make sure alpha is a column vector
       M = length(B(:,1));           % First dimension of B is the number 
 
        I = eye(N,class_t);          % Identity matrix.
        L = logical(I);              % Logical identity matrix.

       XX = x(:,ones(1,N,class_t));
       DX = XX-XX';                  % DX contains entries x(k)-x(j).

    DX(L) = ones(N,1,class_t);               % Put 1's one the main diagonal.

         c = alpha.*prod(DX,2);       % Quantities c(j).

        C = c(:,ones(1,N,class_t)); 
        C = C./C';                   % Matrix with entries c(k)/c(j).
   
        Z = 1./DX;                   % Z contains entries 1/(x(k)-x(j))
     Z(L) = zeros(N,1,class_t);              % with zeros on the diagonal.

        X = Z';                      % X is same as Z', but with 
     X(L) = [];                      % diagonal entries removed.
        X = reshape(X,N-1,N);

        Y = ones(N-1,N,class_t);             % Initialize Y and D matrices.
        D = eye(N,class_t);                  % Y is matrix of cumulative sums,
                                     % D differentiation matrices.
for ell = 1:M
        Y   = cumsum([B(ell,:); ell*Y(1:N-1,:).*X]); % Diagonals
        D   = ell*Z.*(C.*repmat(diag(D),1,N) - D);   % Off-diagonals
     D(L)   = Y(N,:);                                % Correct the diagonal
DM(:,:,ell) = D;                                     % Store the current D
end
