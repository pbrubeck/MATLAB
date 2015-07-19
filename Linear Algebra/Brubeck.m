function x=Brubeck(A, b)
n=size(A,1);
D=diag(diag(A), 0);
C=zeros(n+1, n+1);
C(1:n, 1:n+1)=[eye(n)-D\A, D\b];
C(n+1, n+1)=1;
for i=1:7   
    C=C*C;
end
x=C(1:n, n+1);
end