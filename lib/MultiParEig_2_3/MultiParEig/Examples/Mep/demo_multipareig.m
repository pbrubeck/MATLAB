%DEMO_MULTIPAREIG   demo five-parameter eigenvalue problem with 4 x 4 matrices
% 
% We solve a five-parameter eigenvalue problem
%
% A11 x1 = lambda1 A12 x1 + lambda2 A13 x1 + ... + lambda5 A16 x2
% A21 x2 = lambda1 A22 x2 + lambda2 A23 x2 + ... + lambda5 A26 x2
% ...
% A51 x5 = lambda1 A52 x5 + lambda2 A53 x5 + ... + lambda5 A56 x5
%
% where Aij is a random 4x4 matrix
%
% The problem has 1024 eigenvalues
%
% See also: MULTIPAREIG

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision: 8.9.2015

% we generate the matrices of the problem
n = 5;   % 5

A = cell(n,n+1);
for i=1:n
    for j=1:n+1
        A{i,j} = rand(4);
    end
end

% solution of the problem
[lambda,X,Y] = multipareig(A);

% number of eigenvalues
neig = length(lambda)

% test if we really have eigenvalues (we test only the first eigenvalue)
for i=1:n
     r = A{i,1};
     for k=2:n+1
       r = r - lambda(i,k-1)*A{i,k};
     end;
     fprintf('%.7e\n',min(svd(r)));    
end;
