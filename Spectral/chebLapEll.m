function [Dxx, Jx, Dyy, Jy, x, y] = chebLapEll(a, b, N, M)
% Matrices to calculate the laplacian on the ellipse
% usage: Jx*Lap(U)+Lap(U)*Jy=Dxx*U+U*Dyy'
% domain x:[-xi0,xi0], y:(0,2*pi]
f=sqrt(a^2-b^2);
xi0=acosh(a/f);
[Dyy,y]=fourD2(M);
[Dx,x]=chebD(N);
x=xi0*x; Dx=Dx/xi0; Dxx=Dx*Dx;
Jx=diag(f^2/2*cosh(2*x));
Jy=diag(-f^2/2*cos(2*y));
end
