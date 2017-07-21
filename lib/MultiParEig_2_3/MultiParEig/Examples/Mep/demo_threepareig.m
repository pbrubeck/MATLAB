%DEMO_THREEPAR_EIG   demo nonsingular three-parameter eigenvalue problems with 2 x 2 matrices
%
% We solve a three-parameter eigenvalue problem
%
% A1 x1 = lambda B1 x1 + mu C1 x1 + eta D1 x1 
% A2 x2 = lambda B2 x2 + mu C2 x2 + eta D2 x2 
% A3 x3 = lambda B3 x3 + mu C3 x3 + eta D3 x3 
%
% where 
%
% A1 = [1  2; 3  4]; B1 = [3  1; -1 1]; C1 = [2  1; 5 1], D1 = [ 2 3; 4 -5];
% A2 = [1 -2; 3 -5]; B2 = [1 -1; -2 3]; C2 = [2 -1; 3 1], D2 = [ 1 2; 5 -1];
% A3 = [-1 1; 3  1]; B3 = [2 -1; -1 3]; C3 = [4  1; 1 1], D3 = [-1 2; 2 -1];
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
% The output should include
%
% eigenvalues =
%
%    3.7403            -0.9296            -0.0595          
%   -1.8708             3.9018            -2.9162          
%   -1.1879 + 0.6087i   0.1447 - 0.3354i  -0.5016 + 0.2853i
%   -1.1879 - 0.6087i   0.1447 + 0.3354i  -0.5016 - 0.2853i
%    0.0303 + 0.3501i   0.4042 - 0.8872i   0.2823 - 0.0022i
%    0.0303 - 0.3501i   0.4042 + 0.8872i   0.2823 + 0.0022i
%    0.0794            -0.2437             0.5393          
%   -0.3718             1.5605            -1.0989          
%
% See also: THREEPAREIG

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision: 8.9.2015

% matrices of the two-parameter eigenvalue problem
A1 = [1  2; 3  4]; 
B1 = [3  1; -1 1]; 
C1 = [2  1; 5 1]; 
D1 = [ 2 3; 4 -5];
A2 = [1 -2; 3 -5]; 
B2 = [1 -1; -2 3]; 
C2 = [2 -1; 3 1];
D2 = [ 1 2; 5 -1];
A3 = [-1 1; 3  1]; 
B3 = [2 -1; -1 3]; 
C3 = [4  1; 1 1];
D3 = [-1 2; 2 -1];

% we solve the two-parameter eigenvalue problem
[lambda,mu,eta,Xr,Yr,Zr,Xl,Yl,Zl] = threepareig(A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3);

eigenvalues = [lambda mu eta]

% check that eigenvalues and eigenvectors are correct
for k = 1:length(eigenvalues)
    minsing1 = min(svd((A1-lambda(k)*B1-mu(k)*C1-eta(k)*D1)));
    minsing2 = min(svd((A2-lambda(k)*B2-mu(k)*C2-eta(k)*D2)));
    minsing3 = min(svd((A3-lambda(k)*B3-mu(k)*C3-eta(k)*D3)));
    normres1 = norm((A1-lambda(k)*B1-mu(k)*C1-eta(k)*D1)*Xr(:,k));
    normres2 = norm((A2-lambda(k)*B2-mu(k)*C2-eta(k)*D2)*Yr(:,k));
    normres3 = norm((A3-lambda(k)*B3-mu(k)*C3-eta(k)*D3)*Zr(:,k));
    normres4 = norm(Xl(:,k)'*(A1-lambda(k)*B1-mu(k)*C1-eta(k)*D1));
    normres5 = norm(Yl(:,k)'*(A2-lambda(k)*B2-mu(k)*C2-eta(k)*D2));
    normres6 = norm(Zl(:,k)'*(A3-lambda(k)*B3-mu(k)*C3-eta(k)*D3));
    fprintf('Minimal singular values : (%7.1e, %7.1e, %7.1e), right residuals: (%7.1e, %7.1e, %7.1e), left residuals: (%7.1e,%7.1e, %7.1e)\n',...
        minsing1,minsing2,minsing3,normres1,normres2,normres3,normres4,normres5,normres6)
end
    

