function X = sylvreal(Q,R,U,S,C)

%SYLVREAL  Solves the real Sylvester equation
%
% X = SYLVREAL(Q,R,U,S,C) solves the real Sylvester equation
% A*X + X*B = C, where A = Q*R*Q' and B = U*S*U', Q,U are orthogonal
% and R,S are real quasi upper triangular matrices 
% 
% In Matlab below 2014a package lapack from MatlabCentral is required
%
% See also: SYLV.

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision: 8.9.2015

m = size(R,1);
n = size(S,1);
F = Q'*C*U;

if verLessThan('matlab', '8.3'); 
    % Routine DTRSYL from Lapack is used
    M = lapack('dtrsyl','N','N',1,m,n,R,m,S,n,F,m,1,0);
    X = Q*M{10}*U';
else
    M = sylvester(R,S,F);
    X = Q*M*U';
end