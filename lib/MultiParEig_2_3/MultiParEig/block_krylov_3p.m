function V = block_krylov_3p(B,C,F,k,A)

%BLOCK_KRYLOV_3P   Orthogonal basis for generalized block Krylov subspace
% 
% V = BLOCK_KRYLOV_3P(B,C,F,k,A) returns orthogonal basis for the
% generalized block subspace, where 
% 
% k=1 : V = Lin(F, B1*F, C1*F)
% k=2 : V = Lin(F, B1*F, C1*F, B1^2*F, B1*C1*F, C1*B1*F, C1^2*F),
% k=3 : V = Lin(F, B1*F, C1*F, ..., B1^3*F, ... , C1^3*F), ..., 
%
% where B1 = inv(A)*B and C1 = inv(A)*C (if A is not supplied, B1=B, C1=C)
%
% See also BLOCK_KRYLOV.

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision: 8.9.2015

opts.rankeps = 100*eps;
opts.heuristic = 0;
matrixA = (nargin == 5);
    
[rk,U,tildaS,tildaV] = extended_svd(F,opts);
Q = U(:,1:rk);
V = Q;

zac(1) = 1;
kon(1) = rk;

for j = 1:k
    W = [B*Q C*Q];
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
    