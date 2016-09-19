function [A1,B1,C1,A2,B2,C2,lambda,mu] = random_sparse_2ep(n,dens)

%RANDOM_SPARSE_2EP   random sparse two-parameter eigenvalue problem
%
%[A1,B1,C1,A2,B2,C2,lambda,mu] = RANDOM_SPARSE_2EP(n,dens) builds
% a random two-parameter eigenvalue problem
%   A1 x = l B1 x + u C1 x
%   A2 y = l B2 y + u C2 y
% with sparse matrices of size n x n and calculates eigenvalues (lambda,mu)
%
% Input:
%   n : size of matrices
%   dens : density parameter used in the construction
%
% The obtaines two-parameter eigenvalue problem is similar to a 
% two-parameter eigenvalue problem with diagonal matrices and this 
% enables us to compute all the eigenvalues aand use them for tests
% in numerical methods

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision: 8.9.2015

% we want the same example
if verLessThan('matlab', '7.7')
    rand('state',1);
else
    rng('default')
    rng(0);
end

U1 = sprand(n,n,dens)+speye(n); 
U2 = sprand(n,n,dens)+speye(n);
V1 = sprand(n,n,dens)+speye(n); 
V2 = sprand(n,n,dens)+speye(n);
a1 = rand(n,1)+1i*rand(n,1)-0.5-0.5i;
b1 = rand(n,1)+1i*rand(n,1)-0.5-0.5i;
c1 = rand(n,1)+1i*rand(n,1)-0.5-0.5i;
a2 = rand(n,1)+1i*rand(n,1)-0.5-0.5i;
b2 = rand(n,1)+1i*rand(n,1)-0.5-0.5i;
c2 = rand(n,1)+1i*rand(n,1)-0.5-0.5i;
lambda = zeros(n*n,1);
mu = zeros(n*n,1);

% in order to obtain all eigenvalues we need to solve two-parameter problem with diagonal matrices
% eigenvalues are intersections of lines a1+x b1+y c2=0 and a2+x b2+y c3=0

for k=1:n
   for j=1:n
      tmp = [b1(k),c1(k); b2(j), c2(j)]\[a1(k);a2(j)];
      lambda((k-1)*n+j) = tmp(1);
      mu((k-1)*n+j) = tmp(2);
   end
end

% we multiply matrices by nonsingular matrices so that they are not diagonal

A1 = V1*spdiag(a1)*U1;
B1 = V1*spdiag(b1)*U1;
C1 = V1*spdiag(c1)*U1;
A2 = V2*spdiag(a2)*U2;
B2 = V2*spdiag(b2)*U2;
C2 = V2*spdiag(c2)*U2;

end

function A = spdiag(d);

n = length(d);
A = spconvert([1:n; 1:n; d']');

end
