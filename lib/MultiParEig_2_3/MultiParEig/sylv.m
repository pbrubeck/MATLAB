function X = sylv(Q,R,U,S,C,uselapack)

%SYLV   Solves the complex Sylvester equation
%
% X = SYLV(Q,R,U,S,C,uselapack) solves the Sylvester equation
% A*X + X*B = C, where A = Q*R*Q' and B = U*S*U', and 
% a) Q,U unitary and R,S upper triangular matrices, or
% b) Q,U real orthogonal and R,S real quasi-upper triagular matrices
%
% uselapack = 0: use function sylvester from Matlab (default), 1: use
%    lapack (package is required), 2: use alterative Bartels-Stewart, in 
%    this case matrices R,S have to be triagular. Options 1 and 2 are
%    possible only if Matlab below 2014a is used
% 
% Function assumes that all matrics have the appropriate form

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% BP 02.12.2016: modified to be precision-independent, merge with sylvreal  
% Last revision: 2.12.2016

% Validate number of input parameters
narginchk(5,6);

class_t = superiorfloat(Q,R,U,S,C);

% Make sure all inputs are of the same numeric type.
if ~isa(Q,class_t), Q = numeric_t(Q,class_t); end
if ~isa(R,class_t), R = numeric_t(R,class_t); end
if ~isa(U,class_t), U = numeric_t(U,class_t); end
if ~isa(S,class_t), S = numeric_t(S,class_t); end
if ~isa(C,class_t), C = numeric_t(C,class_t); end

if isreal(Q) && isreal(R) && isreal(U) && isreal(S) && isreal(C)
    realSys = 1; 
else
    realSys = 0;
end

if nargin<6  || isempty(uselapack), uselapack = 0; end

% If the version of Matlab is 2014a or above, lapack is not needed
% and the builtin function sylvester is used
% If mp is used, sylvester ia available in all versions of Matlab
if ~verLessThan('matlab', '8.3') || strcmpi(class_t,'mp'), uselapack = 0; end

m = size(R,1);
n = size(S,1);
F = Q'*C*U;

if uselapack == 0
   X = sylvester(R,S,F);
elseif uselapack == 1 
    if realSys
        % Routine DTRSYL from Lapack is used for a real problem
        M = lapack('dtrsyl','N','N',1,m,n,R,m,S,n,F,m,1,0);
    else
        % Routine ZTRSYL from Lapack is used for a complex problem
        M = lapack('ztrsyl','N','N',1,m,n,R,m,S,n,F,m,1,0);
    end
    X = M{10}; 
else
    % Alternative version without calling Lapack for Matlab below 2014a
    X = zeros(m,n);
    X(:,1) = (R + S(1,1)*eye(m))\F(:,1);
    for k = 2:n
        X(:,k) = (R + S(k,k)*eye(m))\(F(:,k) - X(:,1:k-1)*S(1:k-1,k));
    end
end

X = Q*X*U';
