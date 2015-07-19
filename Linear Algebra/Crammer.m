function x=Crammer(A, b)
% Solves Ax=b by Crammer's Rule.
n=size(A,1);
x=zeros(n,1);
for i=1:n
    B=A;
    B(:,i)=b;
    x(i)=det(B);
end
x=x/det(A);
end