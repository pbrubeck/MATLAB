function [Dxx, Dyy, xx, yy] = chebLapCart(N)
% Matrices to calculate the laplacian on the square [-1,1]^2
% usage: Lap(U)=Dxx*U+U*Dyy'
% domain x:[-1,1], y:[-1,1]
[D, x]=chebD(N);
D2=D^2;
Dxx=D2;
Dyy=D2;
[xx,yy]=meshgrid(x);
end