function [D, x] = fdmfwd2(a, b, n)
% Finite difference matrix, foward difference scheme O(h^2)
h=(b-a)/(n-1);
x=linspace(a,b,n); x=x(:);
r=zeros(n,1); 
r(1:3)=[-3,4,-1]/(2*h);
c=zeros(n,1);
c(1)=r(1);
D=toeplitz(c,r);
D(n-1:n,n-3:n)=[0,-1,0,1;0,1,-4,3]/(2*h);
end