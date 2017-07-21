%DEMO_MULTIPAREIG_MP   demo four-parameter eigenvalue problem with 4 x 4 matrices
% 
% We solve a four-parameter eigenvalue problem with random 4x4 matrices,
% that has 256 eigenvalues
%
% This example requires Multiprecision Computing Toolbox for MATLAB, see
% http://www.advanpix.com/
%
% See also: DEMO_MULTIPAREIG, MULTIPAREIG

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% P. Holoborodko, Advanpix LLC.
% FreeBSD License, see LICENSE.txt

% PH 24.11.2016: extented accuracy check for all equations.
% PH 26.11.2016: fixed bug in second loop boundary.
% BP 27.11.2016: report error if MCT is not installed

% Last revision: 27.11.2016

is_numeric_type_supported('mp'); % check if MCT is installed

% we generate the matrices of the problem in mp format
n = 4;   
class_t = 'mp'; 

A = cell(n,n+1);
for i=1:n
    for j=1:n+1
        A{i,j} = rand(4,class_t);
    end
end

% solution of the problem
[lambda,X,Y] = multipareig(A);

% number of eigenvalues
neig = length(lambda)

% test if we really have eigenvalues (just the first eigenvalue)
disp('Minimal singular values of equations for the first eigenvalue (should be zero):')
for i=1:n
     r = A{i,1};
     for k=2:n+1
       r = r - lambda(i,k-1)*A{i,k};
     end;
     fprintf('%.7e\n',min(svd(r)));    
end;
