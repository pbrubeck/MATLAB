%DEMO_SINGULAR_TWOPAREIG  demo singular two-parameter eigenvalue problems with 2 x 2 matrices
%
% We solve a two-parameter eigenvalue problem
%
% A1 x = lambda B1 x + mu C1 x 
% A2 y = lambda B2 y + mu C2 y,
%
% where 
%
% A1 = [1  2; 3  4]; B1 = [1  1; -1 1]; C1 = [2  1; 5 1],
% A2 = [1 -2; 3 -5]; B2 = [1 -1; -2 3]; C2 = [1 -1; 3 1].
%
% In this case Delta0 = kron(B1,C2) - kron(C1,B2) is singular and we solve
% this as a singular two-parameter eigenvalue problem
% 
% The output should include:
%
% eigenvalues =
%
%   0.2662 + 0.4284i   0.9071 - 0.9634i
%   0.2662 - 0.4284i   0.9071 + 0.9634i
%  -0.4244            -0.0130          
%
% See also: TWOPAREIG

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision: 17.9.2016

% matrices of the two-parameter eigenvalue problem
A1 = [1  2;  3  4]; 
B1 = [1  1; -1  1]; 
C1 = [2  1;  5  1];
A2 = [1 -2;  3 -5]; 
B2 = [1 -1; -2  3]; 
C2 = [1 -1;  3  1];

% Delta0 operator determinant
Delta0 = kron(B1,C2) - kron(C1,B2)

% rank of matrix Delta0 is 3 -> Delta0 is singular
ran = rank(Delta0)

% we solve this as a singular two-parameter eigenvalue problem
opts.singular = 1;
[lambda,mu,Xr,Yr,Xl,Yl] = twopareig(A1,B1,C1,A2,B2,C2,opts);

eigenvalues = [lambda mu]

% check that eigenvalues and eigenvectors are correct
for k = 1:size(eigenvalues,1)
    minsing1 = min(svd((A1-lambda(k)*B1-mu(k)*C1)));
    minsing2 = min(svd((A2-lambda(k)*B2-mu(k)*C2)));
    normres1 = norm((A1-lambda(k)*B1-mu(k)*C1)*Xr(:,k));
    normres2 = norm((A2-lambda(k)*B2-mu(k)*C2)*Yr(:,k));
    normres3 = norm(Xl(:,k)'*(A1-lambda(k)*B1-mu(k)*C1));
    normres4 = norm(Yl(:,k)'*(A2-lambda(k)*B2-mu(k)*C2));
    fprintf('Minimal singular values : (%7.1e, %7.1e), right residuals: (%7.1e, %7.1e), left residuals: (%7.1e,%7.1e)\n',...
        minsing1,minsing2,normres1,normres2,normres3,normres4)
end
    

