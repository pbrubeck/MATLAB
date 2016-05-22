function [X] = kronsolve(A, B, C)
% solves the Sylvester equation A*X+X*B=C
Z=zeros(size(C));
I=eye(size(A));
[Q, U]=schur(A);
[P, V]=schur(B');
G=Q'*C*P;
n=size(B,1);
for i=n:-1:1
    Z(:,i)=(U+V(i,i)*I)\(G(:,i)-Z(:,i+1:end)*V(i,i+1:end)');
end
X=Q*Z*P';
end