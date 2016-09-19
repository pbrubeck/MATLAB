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

% Last revision: 8.9.2015

opts.rankeps = 100*eps;
opts.heuristic = 0;
matrixA = (nargin == 4);
   
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
    