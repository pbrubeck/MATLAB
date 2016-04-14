function [D, x] = fdmcntr(a, b, n)
% Finite difference matrix, central difference scheme O(h^2)
h=(b-a)/(n-1);
x=linspace(a,b,n); x=x(:);
r=zeros(n,1);
r(2)=1/(2*h);
c=zeros(n,1);
c(2)=-1/(2*h);
D=toeplitz(c,r);
D(1,1:2)=[-1,1]/h;
D(n,n-1:n)=[-1,1]/h;
end