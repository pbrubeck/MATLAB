function [ran,U,S,V] = extended_svd(A, opts)

%EXTENDED_SVD  SVD and numerical rank estimation
%
% [ran,U,S,V] = EXTENDED_SVD(A, opts)
% returns SVD of a matrix A and estimates its numerical rank
%
% Output:
%   ran: rank of matrix A
%   U,S,V: factors from [U,S,V] = svd(A);
%
% See also: NUMRANK.

% MultiParEig toolbox
% B. Plestenjak and A. Muhic, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% BP 21.07.2015 : Sometimes svd does not converge, a possible solution is
% to try svd on the transposed matrix and transpose back the results

% Last revision: 06.09.2015

if nargin<2, opts=[]; end

[m1, m2] = size(A);
opts.m1 = m1;
opts.m2 = m2;

if nargout <= 1
    S = svd(A);
    d = S;
else
    [U,S,V] = SVDZero(A);
    if min(size(S))==1
       d = S(1,1);
    else
       d = diag(S);
    end
end
    
if isempty(d) 
    ran = 0;
else
    ran = numrank(d, opts);
end

end

%%-------------------------------------------------------------------------
%% Auxiliary function SVDZero
%%-------------------------------------------------------------------------

function [BigU,BigS,BigV] = SVDZero(A)
        
%SVDZero Optimized SVD for matrix with many zero columns of rows
% [U,S,V] = SVDZero(A) returns singular value decomposition of matrix A. 
% It is more efficient than SVD if matrix has many zero columns or rows
        
% Bor Plestenjak 2014
[n1,n2] = size(A);
        
% we check for zero columns and rows only if matrix is large enough
% otherwise we use standard svd
if min(n1,n2) < 100
    [BigU, BigS, BigV] = SVDCatch(A);
    return
end
        
col = sum(abs(A));
row = sum(abs(A'));
        
rind1 = find(row>0);
rind2 = find(row==0);
r1 = length(rind1);
        
cind1 = find(col>0);
cind2 = find(col==0);
c1 = length(cind1);
        
% we do reduction only if we save enough
if min(r1,c1)>0.7*min(n1,n2)
    [BigU, BigS, BigV] = SVDCatch(A);
    return
end

revr([rind1 rind2])=1:n1;
revc([cind1 cind2])=1:n2;
        
[U,S,V] = SVDCatch(A(rind1,cind1));
        
BigU = [U zeros(r1,n1-r1); zeros(n1-r1,r1) eye(n1-r1)];
BigU = BigU(revr,:);
        
BigV = [V zeros(c1,n2-c1); zeros(n2-c1,c1) eye(n2-c1,n2-c1)];
BigV = BigV(revc,:);
        
BigS = [S zeros(r1,n2-c1); zeros(n1-r1,n2)];
        
end

%%-------------------------------------------------------------------------
%% Auxiliary function SVDCatch
%%-------------------------------------------------------------------------

function [U,S,V] = SVDCatch(A)

%SVDCatch SVD with error control for cases when Matlab fails to return SVD 

% Bor Plestenjak 2015
k = 0; 
success = 0;

while (k<5) && (~success)
    if k>0
        fprintf('SVD failed, try %d\n',k+1)
    end
    try
        [U,S,V] = svd(A);
        success = 1;
    catch err1
        try 
            [V,S,U] = svd(A');
            S = S';
            success = 1;
        end
        k = k + 1;
        A = 1.1234 * A;
    end
end
end


