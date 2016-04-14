function [D, x] = fdmcntr2(a, b, n)
% Finite difference matrix, central difference scheme O(h^4)
h=(b-a)/(n-1);
x=linspace(a,b,n); x=x(:);
r=zeros(n,1);
r(2:3)=[8,-1]/(12*h);
c=zeros(n,1);
c(2:3)=[-8,1]/(12*h);
D=toeplitz(c,r);
D(1:2,1:4)=[-3,4,-1,0;-1,0,1,0]/(2*h);
D(n-1:n,n-3:n)=[0,-1,0,1;0,1,-4,3]/(2*h);
end