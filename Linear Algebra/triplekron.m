function [ X ] = triplekron(A, F)
n=size(A,1); I=eye(n);
[Q, U]=schur(A);
FF=ge3kv(Q, F);
Z=zeros(n, n, n);
for i=n:-1:1
    for j=n:-1:1
        Z(:,j,i)=(U+(U(i,i)+U(j,j))*I)\(FF(:,j,i)-Z(:,j+1:end,i)*U(j,j+1:end)');
    end
    FF=FF-reshape(kron(U(:,i), Z(:,:,i)), size(FF));
end
X=ge3kv(Q', Z);
end

function [x] = ge3kv(A, x)
n=size(A,1); n2=n*n;
x=reshape(A*reshape(x, [n n2]), size(x));
x=permute(x, [2 3 1]);
x=reshape(A*reshape(x, [n n2]), size(x));
x=permute(x, [2 3 1]);
x=reshape(A*reshape(x, [n n2]), size(x));
x=permute(x, [2 3 1]);
end