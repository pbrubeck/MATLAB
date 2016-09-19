function X = sylv(Q,R,U,S,C,uselapack)

%SYLV   Solves the complex Sylvester equation
%
% X = SYLV(Q,R,U,S,C,uselapack) solves the complex Sylvester equation
% A*X + X*B = C, where A = Q*R*Q' and B = U*S*U', Q,U are unitary and 
% R,S are upper triangular complex matrices 
%
% uselapack = 1: uses faster package lapack from MatlabCentral (default 0)
% (relevant only for Matlab below 2014a)
% 
% See also: SYLVREAL.

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision: 8.9.2015

if nargin<5  || isempty(uselapack), uselapack = 0; end

% If the version of Matlab is 2014a or above, lapack is not needed
% and the builtin function sylvester is used
if ~verLessThan('matlab', '8.3'); uselapack = 1; end

m = size(R,1);
n = size(S,1);
X = zeros(m,n);
F = Q'*C*U;

if uselapack
    % Faster version that uses Lapack
    if verLessThan('matlab', '8.3'); 
        M = lapack('ztrsyl','N','N',1,m,n,R,m,S,n,F,m,1,0);
        X = M{10};
    else
        M = sylvester(R,S,F);
        X = M;
    end
else
    % Alternative version without calling Lapack for Matlab below 2014a
    X(:,1) = (R + S(1,1)*eye(m))\F(:,1);
    for k = 2:n
        X(:,k) = (R + S(k,k)*eye(m))\(F(:,k) - X(:,1:k-1)*S(1:k-1,k));
    end
end

X = Q*X*U';
