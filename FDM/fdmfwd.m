function [D, x] = fdmfwd(a, b, n)
% Finite difference matrix, foward difference scheme O(h)
h=(b-a)/(n-1);
x=linspace(a,b,n); x=x(:);
r=zeros(n,1); 
r(1:2)=[-1,1]/h;
c=zeros(n,1);
c(1)=r(1);
D=toeplitz(c,r);
D(n,n-1:n)=D(n-1,n-1:n);
end