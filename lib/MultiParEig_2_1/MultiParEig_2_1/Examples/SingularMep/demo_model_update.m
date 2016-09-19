%DEMO_MODEL_UPDATE   demo model updating problem with two parameters from Cottin
%
% For a given matrices A, B, and C of the same size and values s1 and s2
% we are looking for lambda and mu such that A - lambda*B - mu*C has 
% eigenvalues s1 and s2 (we do not care about the remaining eigenvalues)
% 
% We can write this as a singular two-parameter eigenvalue problem
% 
% (A - s1 I) x = lambda B x + mu C x 
% (A - s2 I) y = lambda B y + mu C y 
%
% In the generic case there are n(n-1)/2 solutions.
%
% References:
%  - N. Cottin. Dynamic model updating — a multiparameter eigenvalue problem, 
%    Mech. Syst. Signal Pr. 15 (2001) 649–665.
%  - A. Muhic, B. Plestenjak, On the singular two-parameter eigenvalue problem, 
%    Electron. J. Linear Algebra 18 (2009) 420–437.

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision: 8.9.2015

if verLessThan('matlab', '7.7')
    rand('state',1);
else
    rng('default')
    rng(1);
end

n = 10; % computation is feasible for approx. n < 50

A = rand(n);
B = rand(n);
C = rand(n);
I = eye(n);

% prescribed eigenvalues are 2 and 3
s1 = 2;
s2 = 3;

% we solve a singular two-parameter eigenvalue problem
opts.singular = 1;
[lambda,mu] = twopareig(A-s1*I,B,C,A-s2*I,B,C,opts);

% all possible solutions
solution = [lambda mu];

% number of solutions
solutions_found = length(solution)

% we write first 10 solutions
first10 = solution(1:min([10 length(solution)]),:)

% check that s1 and s2 are really eigenvalues of A-lambda*B-mu*C 
% (for first 10 solutions)
for k = 1:min([10 length(solution)])
    minsing1 = min(svd((A-lambda(k)*B - mu(k)*C - s1*I)));
    minsing2 = min(svd((A-lambda(k)*B - mu(k)*C - s2*I)));
    fprintf('Minimal singular values for A - si*I - lambda*B - mu C: (%7.1e, %7.1e), \n',minsing1,minsing2)
end

% another check
eig1 = eig(A-lambda(1)*B-mu(1)*C)
