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
A = cell(5,6);
for i=1:5
    for j=1:6
        A{i,j} = rand(4);
    end
end

% solution of the problem
[lambda,X,Y] = multipareig(A);

% number of eigenvalues
neig = length(lambda)

% test if we really have eigenvalues (we test the first eigenvalue)
minsing = min(svd(A{1,1} - lambda(1,1)*A{1,2} - lambda(1,2)*A{1,3} - lambda(1,3)*A{1,4} - lambda(1,4)*A{1,5} - lambda(1,5)*A{1,6}))

