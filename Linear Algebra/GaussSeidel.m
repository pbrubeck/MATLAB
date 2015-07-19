function x=GaussSeidel(A, b, x)
% Solves Ax=b by clever iteration.
n=size(A,1);
D=diag(diag(A),0);
A=eye(n)-D\A;
b=D\b;
ii=0;
err=1;
while(err>=1E-10 && ii<30)
    x0=x;
    for j=1:n
        x(j)=A(j,:)*x+b(j);
    end
    err=max(abs((x-x0)./x));
    ii=ii+1;
end