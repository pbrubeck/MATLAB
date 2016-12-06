function [ c ] = hermSeries(f, n)
sig=1/sqrt(2);
[x,w]=gauherm(n,0,sig);
Psi=zeros(n);
a=zeros(n,1); a(end)=1;
for i=1:n
    Psi(i,:)=HermitePsi(a(end-i+1:end),x);
end
w=w.*exp(x.^2/(2*sig^2));
c=Psi*diag(w)*f(x(:));
end