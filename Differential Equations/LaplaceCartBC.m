function [ uu ] = LaplaceCartBC(N)
%LAPLACECARTBC Summary of this function goes here
%   Detailed explanation goes here
[D,x]=chebD(N);
[xx, yy]=meshgrid(x); xx=xx(:); yy=yy(:);
D2=D^2; I=eye(N);
L=kron(D2, I)+kron(I, D2);

b=find(abs(xx)==1 | abs(yy)==1);
L(b,:)=zeros(4*(N-1), N^2);
L(b,b)=eye(4*(N-1));
rhs=zeros(N^2, 1);
rhs(b)=(yy(b)==1 & xx(b)<0).*sin(pi*xx(b)).^4+0.2*(xx(b)==1).*sin(3*pi*yy(b));
u=L\rhs; uu=reshape(u, N, N);
surf(x, x, uu);
end

