function [x, DM] = chebdif_mp(N, M, class_t)

%CHEBDIF_MP Differentiation matrices on Chebyshev nodes
%
% [x, DM] =  chebdif(N,M,class_t) computes the differentiation 
% matrices D1, D2, ..., DM on Chebyshev nodes of class_t. 
% 
% Input:
% N:        Size of differentiation matrix.        
% M:        Number of derivatives required (integer),  0 < M <= N-1.
% class_t:  numeric type to use (optional), default is 'double'
%
% Output:
% x:        Chebyshev nodes
% DM:       DM(1:N,1:N,ell) contains ell-th derivative matrix, ell=1..M.
%
% The code implements two strategies for enhanced accuracy suggested by 
% W. Don and S. Solomonoff in SIAM J. Sci. Comp. Vol. 6, pp. 1253--1268 
% (1994). The two strategies are (a) the use of trigonometric identities 
% to avoid the computation of differences x(k)-x(j) and (b) the use of the
% "flipping trick" which is necessary since sin t can be computed to high
% relative precision when t is small whereas sin (pi-t) cannot.
% Note added May 2003:  It may, in fact, be slightly better not to
% implement the strategies (a) and (b).   Please consult the following
% paper for details:   "Spectral Differencing with a Twist", by
% R. Baltensperger and M.R. Trummer, to appear in SIAM J. Sci. Comp. 
%
% This function was taken from Matlab package DMSUITE by J.A.C Weideman
% and modified to be precision-independent, i.e. be able to work with any 
% numeric type (e.g. built-in 'double'/'single' or custom like 'mp').

% J.A.C. Weideman, S.C. Reddy 1998. Help notes modified by JACW, May 2003.

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% BP 28.11.2016: modified to be precision-independent
% Last revision: 28.11.2016

narginchk(2, 3);

if nargin < 3, 
    class_t = 'double'; 
end

     I = eye(N);                          % Identity matrix.     
     L = logical(I);                      % Logical identity matrix.

    n1 = floor(N/2); n2  = ceil(N/2);     % Indices used for flipping trick.

     k = [0:N-1]';                        % Compute theta vector.
    mp_pi = numeric_t('pi',class_t); 
    th = k*mp_pi/(N-1);

     x = sin(mp_pi*[N-1:-2:1-N]'/(2*(N-1))); % Compute Chebyshev points.

     T = repmat(th/2,1,N);                
    DX = 2*sin(T'+T).*sin(T'-T);          % Trigonometric identity. 
    DX = [DX(1:n1,:); -flipud(fliplr(DX(1:n2,:)))];   % Flipping trick. 
 DX(L) = ones(N,1,class_t);               % Put 1's on the main diagonal of DX.

     C = toeplitz(numeric_t('-1',class_t).^k);  % C is the matrix with 
C(1,:) = C(1,:)*2; C(N,:) = C(N,:)*2;     % entries c(k)/c(j)
C(:,1) = C(:,1)/2; C(:,N) = C(:,N)/2;

     Z = 1./DX;                           % Z contains entries 1/(x(k)-x(j))  
  Z(L) = zeros(N,1,class_t);              % with zeros on the diagonal.

     D = eye(N,class_t);                  % D contains diff. matrices.
                                          
for ell = 1:M
          D = ell*Z.*(C.*repmat(diag(D),1,N) - D); % Off-diagonals
       D(L) = -sum(D');                            % Correct main diagonal of D
DM(:,:,ell) = D;                                   % Store current D in DM
end
