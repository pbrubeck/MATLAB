function [L, U, x] = Crout(A, b)
% Returns the LU decomposition of A and optionally solves Ax=b.
n=size(A,1);
L=zeros(n,n);
U=eye(n);
for i=1:n
    for j=1:i
        L(i,j)=A(i,j)-L(i,1:j-1)*U(1:j-1,j);
    end
    for j=i+1:n
        U(i,j)=(A(i,j)-L(i,1:i-1)*U(1:i-1,j))/L(i,i);
    end
end
if nargin > 1
    x=bwdSubst(U, fwdSubst(L, b));
end
end