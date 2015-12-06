function [X] = kronsolve(A, B, C)
% solves the Sylvester equation A*X+X*B=C
I=eye(size(A));
X=zeros(size(C));
[Q1, U1]=schur(A);
[Q2, U2]=schur(B');
CC=Q1'*C*Q2;
n=size(B, 1);
for i=n:-1:1
    X(:, i)=(U1+U2(i,i)*I)\(CC(:,i)-X(:,i+1:end)*U2(i,i+1:end)');
end
X=Q1*X*Q2';
disp(norm(A*X+X*B-C, 'fro'))
end