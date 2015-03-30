function x=GaussSeidel(A, b, x)
n=size(A,1);
D=diag(diag(A),0);
A=eye(n)-D\A;
b=D\b;
ii=0;
err=1;
while(err>=0.02 && ii<20)
    x0=x;
    for j=1:n
        x(j)=A(j,:)*x+b(j);
    end
    err=max(abs((x-x0)./x));
    ii=ii+1;
    fprintf('i=%d\t err=%f\n', ii, err);
end