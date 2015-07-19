function [Q, R, x]=GramSchmidt(A, b)
% Returns an orthonormal basis for C(A), the QR decomposition, and
% optionally solves Ax=b.
n=size(A,2);
Q=A;
for i=1:n
    for j=i+1:n
        Q(:,j)=Q(:,j)-proj(Q(:,j), Q(:,i));
    end
end
for i=1:n
    n=norm(Q(:,i));
    if(n>1E-15)
        Q(:,i)=Q(:,i)/n;
    end
end
R=Q'*A;
if nargin > 1
    x=bwdSubst(R, Q'*b);
end
end