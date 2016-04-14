function [D, x] = fdmbwd2(a, b, n)
% Finite difference matrix, backward difference scheme O(h^2)
h=(b-a)/(n-1);
x=linspace(a,b,n); x=x(:);
c=zeros(n,1); 
c(1:3)=[3,-4,1]/(2*h);
r=zeros(n,1);
r(1)=c(1);
D=toeplitz(c,r);
D(1:2,1:4)=[-3,4,-1,0;-1,0,1,0]/(2*h);
end