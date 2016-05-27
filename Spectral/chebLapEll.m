function [F, G, Duu, Dvv, u, v] = chebLapEll(N,M)
% Matrices to calculate the laplacian on the ellipse
% usage: F*Lap(U)+Lap(U)*G=Duu*U+U*Dvv'
% domain r:[-1,1], phi:(0,2*pi]
[Dvv,v]=periodicD2(M);
[Du,u]=chebD(N);
Duu=Du^2;
F=diag(cosh(2*u)/2);
G=diag(-cos(2*v)/2);
end
