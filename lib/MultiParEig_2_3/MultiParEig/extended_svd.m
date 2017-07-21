function [ran,U,S,V,report] = extended_svd(A, opts)
%EXTENDED_SVD  SVD and numerical rank estimation
%
% [ran,U,S,V,report] = EXTENDED_SVD(A, opts)
% returns SVD of a matrix A and estimates its numerical rank
%
% Output:
%   ran: rank of matrix A
%   U,S,V: factors from [U,S,V] = svd(A);
%   report: additional information: [rank s(1) s(rank) s(rank+1) choice], 
%      where s(i) is i-th singular value
%
% Options in opts:
%    - rrqr (0) : set to 1 or 2 to use rank-revealing qr with column pivoting instead of SVD, where
%         - 1: only orthogonal basis U and ran are computed from qrcp
%         - 2: only orthogonal basis V and ran are computed from qrcp
%
% See also: NUMRANK.

% MultiParEig toolbox
% B. Plestenjak and A. Muhic, University of Ljubljana
% P. Holoborodko, Advanpix LLC.
% FreeBSD License, see LICENSE.txt

% BP 21.07.2015 : Sometimes svd does not converge, a possible solution is to try svd on the transposed matrix and transpose back the results
% BP 15.09.2016 : Additional result extra with choice details, option of fixed rank
% BP 23.09.2016 : option of using qr with column pivoting instead of SVD
% BP 21.11.2016 : faster QRZero for rank revealing qr
% PH 22.11.2016 : precision-independent version
% BP 26.11.2016 : bug fixed
% PH 26.11.2016 : code simplification and cleaning-up.

% Last revision: 26.11.2016

class_t = class(A);

if nargin<2, opts=[]; end
if isfield(opts,'rrqr'), rrqr = opts.rrqr;  else rrqr = 1;  end;

[m1, m2] = size(A);
opts.m1 = m1;
opts.m2 = m2;

if nargout <= 1
    if rrqr
        [Q,R] = QRZero(A); %#ok<*ASGLU>
        d = abs(diag(R));
    else
        S = svd(A);
        d = S;
    end
else
    if rrqr>0
        if rrqr == 1
            [U,S] = QRZero(A);
            V = numeric_t([],class_t);
        else
            [V,S] = QRZero(A');
            U = numeric_t([],class_t);
            S = S';
        end
    else    
        [U,S,V] = SVDZero(A);
    end
    if min(size(S))==1
       d = abs(S(1,1));
    else
       d = abs(diag(S));
    end
end

if isempty(d)
    ran = 0;
    report = numeric_t([0 0 0 0 0],class_t);
else
    [ran,choice] = numrank(d, opts);
    if ran<length(d)
        if ran==0
            report = [ran d(1) 0 0 0];
        else
            report = [ran d(1) d(ran) d(ran+1) choice];
        end
    else
        report = [ran d(1) d(ran) 0 choice];
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
if min(r1,c1)>0.9*min(n1,n2)
    [BigU, BigS, BigV] = SVDCatch(A);
    return
end

revr([rind1 rind2])=1:n1;
revc([cind1 cind2])=1:n2;
        
[U,S,V] = SVDCatch(A(rind1,cind1));
        
BigU = [U zeros(r1,n1-r1,class_t); zeros(n1-r1,r1,class_t) eye(n1-r1,class_t)];
BigU = BigU(revr,:);
        
BigV = [V zeros(c1,n2-c1,class_t); zeros(n2-c1,c1,class_t) eye(n2-c1,n2-c1,class_t)];
BigV = BigV(revc,:);
        
BigS = [S zeros(r1,n2-c1,class_t); zeros(n1-r1,n2,class_t)];
        
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
        A = numeric_t('1.1234',class_t) * A;
    end
end
end

%%-------------------------------------------------------------------------
%% Auxiliary function QRZero
%%-------------------------------------------------------------------------

function [BigU,BigS] = QRZero(A)
        
%QRZero Optimized rank revealing QR for matrix with many zero columns of rows
% [U,S] = QRZero(A) returns QR decomposition of matrix A
% It is more efficient than QR if matrix has many zero columns or rows
        
% Bor Plestenjak 2016
[n1,n2] = size(A);
        
% we check for zero columns and rows only if matrix is large enough
% otherwise we use standard qr
if min(n1,n2) < 100
    [BigU, BigS, tilda] = qr(A);
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
if min(r1,c1)>0.9*min(n1,n2)
    [BigU, BigS, tilda] = qr(A);
    return
end

revr([rind1 rind2])=1:n1;
revc([cind1 cind2])=1:n2;
        
[U,S,tilda] = qr(A(rind1,cind1));
        
BigU = [U zeros(r1,n1-r1,class_t); zeros(n1-r1,r1,class_t) eye(n1-r1,class_t)];
BigU = BigU(revr,:);
               
BigS = [S zeros(r1,n2-c1,class_t); zeros(n1-r1,n2,class_t)];
        
end

end % extended_svd