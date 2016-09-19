function Delta = multipar_delta(A)

%MULTIPAR_DELTA   Delta matrices for a multiparameter eigenvalue problem
%
% Delta = MULTIPAR_DELTA(A) returns set of k+1 Delta matrices, which are 
% operator determinants related to the multiparameter eigenvalue problem
%
% A{1,1} x1 = lambda(1) A{1,2} x1 +  ... + lambda(k) A{1,k+1} x1 
% A{2,1} x2 = lambda(1) A{2,2} x2 +  ... + lambda(k) A{2,k+1} x2 
% ...
% A{k,1} xk = lambda(1) A{k,2} xk +  ... + lambda(k) A{k,k+1} xk 
%
% See also: TWOPAR_DELTA, THREEPAR_DELTA

% MultiParEig toolbox
% B. Plestenjak and A. Muhic, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision: 8.9.2015

k = length(A); % number of parameters + 1 

% computation of Delta matrices
Delta = cell(1,length(A));
for r = 0:k-1
    Delta(r+1) = {(-1)^(r)*krondet(A, r)};
end;

end

%------------------------------------------------------------------------

function krondelta =  krondet(A,index)
% recursive computation of Delta matrices - operator determinants

if nargin == 2
    n = length(A);
    ind = [1:index, (index+2):n];
    A = A(:, ind);
end;
n = length(A);
if n == 1
    krondelta = A{1,1}; 
    return;
end;
deter = 0; s = 1;
for k = 1 : n
    indj = [1:(k-1), (k+1):n];
    indi = 2: n;   
    deter =  deter + s*kron( A{1, k},  krondet(A(indi, indj))  );
    s = -s;
end;

krondelta = deter;
end