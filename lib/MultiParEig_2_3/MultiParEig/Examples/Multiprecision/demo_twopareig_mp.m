%DEMO_TWOPAREIG_MP   demo nonsingular two-parameter eigenvalue problems
%with 2 x 2 matrices using multiple precision
%
% We solve a two-parameter eigenvalue problem
%
% A1 x = lambda B1 x + mu C1 x 
% A2 y = lambda B2 y + mu C2 y,
%
% where 
%
% A1 = [1  2; 3  4]; B1 = [3  1; -1 1]; C1 = [2  1; 5 1],
% A2 = [1 -2; 3 -5]; B2 = [1 -1; -2 3]; C2 = [2 -1; 3 1].
%
% The output should include:
%
% eigenvalues =
%
%   -3.5718             5.6063          
%    3.9014            -1.0824          
%   -0.1364 + 0.0800i   0.0259 + 0.2820i
%   -0.1364 - 0.0800i   0.0259 - 0.2820i
%
% This example requires Multiprecision Computing Toolbox for MATLAB, see
% http://www.advanpix.com/
%
% The residuals obtained in this method are much smaller than residuals
% obtained using the standard double precision (see DEMO_TWOPAREIG).

% See also: DEMO_TWOPAREIG_MP TWOPAREIG

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% BP 27.11.2016: report error if MCT is not installed

% Last revision: 27.11.2016

is_numeric_type_supported('mp'); % check if MCT is installed

% matrices of the two-parameter eigenvalue problem
A1 = mp('[1  2;  3  4]'); 
B1 = mp('[3  1; -1  1]'); 
C1 = mp('[2  1;  5  1]');
A2 = mp('[1 -2;  3 -5]'); 
B2 = mp('[1 -1; -2  3]'); 
C2 = mp('[2 -1;  3  1]');

% Delta0 operator determinant
Delta0 = kron(B1,C2) - kron(C1,B2)

% rank of matrix Delta0 is 4 -> Delta0 is nonsingular
ran = rank(Delta0)

% We solve the two-parameter eigenvalue problem. The computation is done in
% quadruple precision because this is the numeric format for A1,...,C2
[lambda,mu,Xr,Yr,Xl,Yl] = twopareig(A1,B1,C1,A2,B2,C2);

eigenvalues = [lambda mu]

% check that eigenvalues and eigenvectors are correct
for k = 1:length(eigenvalues)
    minsing1 = min(svd((A1-lambda(k)*B1-mu(k)*C1)));
    minsing2 = min(svd((A2-lambda(k)*B2-mu(k)*C2)));
    normres1 = norm((A1-lambda(k)*B1-mu(k)*C1)*Xr(:,k));
    normres2 = norm((A2-lambda(k)*B2-mu(k)*C2)*Yr(:,k));
    normres3 = norm(Xl(:,k)'*(A1-lambda(k)*B1-mu(k)*C1));
    normres4 = norm(Yl(:,k)'*(A2-lambda(k)*B2-mu(k)*C2));
    fprintf('Minimal singular values : (%7.1e, %7.1e), right residuals: (%7.1e, %7.1e), left residuals: (%7.1e,%7.1e)\n',...
        minsing1,minsing2,normres1,normres2,normres3,normres4)
end
    

