function x=Jacobi(A, b, x0)
n=size(A,1);
D=diag(diag(A),0);
A=eye(n)-D\A;
b=D\b;
ii=0;
err=1;
while(err>=1E-10 && ii<20)
    x=A*x0+b;
    err=max(abs((x-x0)./x));
    x0=x;
    ii=ii+1;
    fprintf('i=%d\t err=%f\n', ii, err);
end