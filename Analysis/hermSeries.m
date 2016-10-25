function [ c ] = hermSeries(f, n)
[x,w]=GaussHermite(n,0,1/sqrt(2));
Psi=zeros(n);
for i=1:n
    a=[zeros(i-1,1); 1];
    Psi(i,:)=HermitePsi(a,x);
end
w=w.*exp(x.^2);
c=Psi*diag(w)*f(x(:));
end