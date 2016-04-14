function [D, x] = fdmbwd(a, b, n)
% Finite difference matrix, backward difference scheme O(h)
h=(b-a)/(n-1);
x=linspace(a,b,n); x=x(:);
c=zeros(n,1);
c(1:2)=[1,-1]/h;
r=zeros(n,1); 
r(1)=c(1);
D=toeplitz(c,r);
D(1,1:2)=D(2,1:2);
end