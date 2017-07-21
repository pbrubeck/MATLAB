function V = block_krylov(B,F,k,A)

%BLOCK_KRYLOV   Orthogonal basis for block Krylov subspace
%
% V = BLOCK_KRYLOV(B,F,k) returns orthogonal basis for the
% block Krylov subspace (F,B*F,B^2*F,...,B^(k-1)*F)
%
% V = BLOCK_KRYLOV(B,F,k,A) returns orthogonal basis for the
% block Krylov subspace (F,inv(A)*B*F,(inv(A)*B)^2*F,...,(inv(A)*B)^(k-1)*F)
%
% See also: BLOCK_KRYLOV_3P.

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% BP 03.12.2016 : modified to be precision-independent
% Last revision: 3.12.2016

narginchk(3, 4);
matrixA = (nargin == 4);

if matrixA
    class_t = superiorfloat(B,F,A);
else
    class_t = superiorfloat(B,F);
end

% Make sure all inputs are of the same numeric type.
if ~isa(B,class_t), B = numeric_t(B,class_t); end;
if ~isa(F,class_t), F = numeric_t(F,class_t); end;
if matrixA && (~isa(A,class_t)), A = numeric_t(A,class_t); end;

opts.rankeps = numeric_t('100*eps',class_t);
opts.heuristic = 0;
   
[rk,U,tildaS,tildaV] = extended_svd(F,opts);
Q = U(:,1:rk);
V = Q;

zac(1) = 1;
kon(1) = rk;

for j = 1:k
    W = B*Q;
    if matrixA
        W = A\W;
    end
    zacn = norm(W,'fro');
    for i = 1:j
        tmp = V(:,zac(i):kon(i))'*W;
        W = W - V(:,zac(i):kon(i))*tmp;
    end
    % Repeated orhtogonalization
    for i=1:j
        tmp = V(:,zac(i):kon(i))'*W;
        W = W - V(:,zac(i):kon(i))*tmp;
    end
    [rk,U,tildaS,tildaV] = extended_svd(W,opts);
    Q = U(:,1:rk);
    if rk == 0
        return
    end
    V = [V Q];
    zac(j+1) = kon(j)+1;
    kon(j+1) = kon(j)+rk;
end
    